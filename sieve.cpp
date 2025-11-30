#include "sieve.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <execution>
#include <immintrin.h>
#include <mutex>
#include <vector>

#ifdef __linux__
    #include <unistd.h>

ulong get_usable_l1_bytes() { return (ulong)sysconf(_SC_LEVEL1_DCACHE_SIZE); }
#elif defined(_WIN32)
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>

using ulong = unsigned long;

ulong get_usable_l1_bytes() {
    DWORD len = 0;
    BOOL ok   = GetLogicalProcessorInformation(nullptr, &len);
    if (!ok) {
        DWORD err = GetLastError();
        if (err != ERROR_INSUFFICIENT_BUFFER) return 0;
    }

    std::vector<uint8_t> buffer(len);
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION info =
        reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION>(buffer.data());
    ok = GetLogicalProcessorInformation(info, &len);
    if (!ok) return 0;

    DWORD count = len / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
    for (DWORD i = 0; i < count; ++i)
        if (info[i].Relationship == RelationCache && info[i].Cache.Level == 1 &&
            info[i].Cache.Type == CacheData && info[i].Cache.Size > 0)
            return (ulong)info[i].Cache.Size;

    return 0;
}
#else
using ulong = unsigned long;

ulong get_usable_l1_bytes() { return 0; }
#endif

static const ulong WHEEL                    = 30UL;
static const std::array<ulong, 8> RESIDUES  = {1UL,  7UL,  11UL, 13UL,
                                               17UL, 19UL, 23UL, 29UL};
static const size_t NUM_RESIDUES            = 8;
static const ulong FALLBACK_USABLE_L1_BYTES = 32'768UL;
static const ulong SEGMENT_SIZE =
    WHEEL *
    (get_usable_l1_bytes() ? get_usable_l1_bytes() : FALLBACK_USABLE_L1_BYTES);
static const ulong LCM_7_11_13_30 = 30UL * 7UL * 11UL * 13UL; // = 30030;

static uint8_t modinv30(uint8_t a) {
    a %= 30;
    for (uint8_t x = 1; x < 30; ++x)
        if ((a * x) % 30 == 1) return x;
    return 0;
}

static void make_offset_table(std::array<std::array<uint8_t, 30>, 30> &table) {
    for (int d = 0; d < 30; ++d)
        for (int inv = 0; inv < 30; ++inv)
            table[d][inv] = (uint8_t)((d * inv) % 30);
}

static std::vector<ulong> simple_sieve(ulong limit) {
    if (limit < 2) return {};
    std::vector<char> is_prime(limit + 1, 1);
    is_prime[0] = is_prime[1] = 0;
    for (ulong p = 2; p * p <= limit; ++p) {
        if (!is_prime[p]) continue;
        for (ulong q = p * p; q <= limit; q += p) is_prime[q] = 0;
    }
    std::vector<ulong> out;
    for (ulong i = 2; i <= limit; ++i)
        if (is_prime[i]) out.push_back(i);
    return out;
}

static std::array<int8_t, 30> make_residue_index() {
    std::array<int8_t, 30> idx;
    idx.fill(-1);
    for (size_t i = 0; i < RESIDUES.size(); ++i) idx[RESIDUES[i]] = (int8_t)i;
    return idx;
}

struct BasePrimeInfo {
    ulong p;
    uint8_t inv30;
};

std::vector<ulong> segmented_wheel_sieve(ulong n) {
    const ulong limit = 2 * n;
    if (limit < 2) return {};

    static const auto residue_index = make_residue_index();
    const ulong pattern_bytes       = (LCM_7_11_13_30 + WHEEL - 1) / WHEEL;

    std::array<std::array<uint8_t, 30>, 30> OFFSET_TABLE;
    make_offset_table(OFFSET_TABLE);

    const ulong L = (ulong)std::floor(std::sqrt((long double)limit));
    std::vector<ulong> base_primes;
    if (L <= SEGMENT_SIZE) base_primes = simple_sieve(L);
    else {
        ulong nprime = (L + 1) / 2;
        auto rec     = segmented_wheel_sieve(nprime);
        base_primes.reserve(rec.size());
        for (ulong p : rec)
            if (p <= L) base_primes.push_back(p);
    }

    std::vector<BasePrimeInfo> base_infos;
    base_infos.reserve(base_primes.size());
    for (ulong p : base_primes)
        if (p >= 17 && (p % 2 != 0) && (p % 3 != 0) && (p % 5 != 0)) {
            uint8_t inv = modinv30((uint8_t)(p % 30));
            base_infos.push_back(BasePrimeInfo{p, inv});
        }

    std::vector<uint8_t> pattern(pattern_bytes, 0x00);
    auto set_composite_in_pattern = [&](ulong number) {
        if (number == 0) return;
        ulong block = number / WHEEL;
        if (block >= pattern_bytes) return;
        ulong r   = number % WHEEL;
        int8_t bi = residue_index[r];
        if (bi >= 0) pattern[block] |= (uint8_t)(1u << bi);
    };

    for (ulong p : {7UL, 11UL, 13UL}) {
        ulong start = p * p;
        if (start >= LCM_7_11_13_30) start = p;
        for (ulong m = start; m < LCM_7_11_13_30; m += p)
            set_composite_in_pattern(m);
    }

    const ulong aligned_low_start = 0UL;

    std::vector<ulong> final_primes;
    final_primes.reserve((size_t)(limit / (std::log(limit ? limit : 2.0L)) +
                                  10));
    if (2 <= limit) final_primes.push_back(2);
    if (3 <= limit) final_primes.push_back(3);
    if (5 <= limit) final_primes.push_back(5);

    std::mutex append_mutex;

    auto process_segment = [&](size_t seg_i) {
        ulong low = aligned_low_start + (ulong)seg_i * SEGMENT_SIZE;
        if (low > limit) return;
        ulong high            = std::min(low + SEGMENT_SIZE, limit + 1),
              segment_numbers = (high > low) ? (high - low) : 0,
              bytes_needed    = (segment_numbers + WHEEL - 1) / WHEEL;
        if (bytes_needed == 0) return;

        std::vector<uint8_t> bits(bytes_needed);
        for (ulong pos = 0; pos < bytes_needed;) {
            ulong copy_len = std::min<ulong>(pattern_bytes, bytes_needed - pos);
            std::memcpy(bits.data() + pos, pattern.data(), (size_t)copy_len);
            pos += copy_len;
        }

        std::for_each(
            std::execution::par_unseq, base_infos.begin(), base_infos.end(),
            [&](const BasePrimeInfo &info) {
                ulong p       = info.p;
                uint8_t inv30 = info.inv30;

                ulong m = (low + p - 1) / p * p;
                if (m < p * p) m = p * p;

                uint8_t m_mod30 = (uint8_t)(m % 30);

                for (size_t ridx = 0; ridx < NUM_RESIDUES; ++ridx) {
                    uint8_t r = (uint8_t)RESIDUES[ridx];
                    int d     = (int)r - (int)m_mod30;
                    if (d < 0) d += 30;
                    uint8_t s = OFFSET_TABLE[d][inv30];

                    ulong first = m + (ulong)s * p;
                    if (first < p * p) {
                        ulong step = 30 * p,
                              k    = (p * p - first + step - 1) / step;
                        first += k * step;
                    }

                    if (first >= high) continue;

                    ulong start_byte  = (first - low) / WHEEL;
                    uint8_t bit_mask  = (uint8_t)(1u << ridx);
                    ulong byte_stride = p;

                    for (ulong idx = start_byte; idx < bytes_needed;
                         idx += byte_stride) {
                        bits[idx] |= bit_mask;
                    }
                }
            }
        );

        std::vector<ulong> local_primes;
        local_primes.reserve(1'024);

        for (ulong byte_i = 0; byte_i < bytes_needed; ++byte_i) {
            uint8_t b = bits[(size_t)byte_i];
            if (b == 0xFFu) continue; // all residues composite
            for (int bit = 0; bit < (int)NUM_RESIDUES; ++bit) {
                if ((b >> bit) & 1u) continue; // composite
                ulong p = low + byte_i * WHEEL + RESIDUES[bit];
                if (p <= 1) continue;
                if (p > limit) continue;
                local_primes.push_back(p);
            }
        }

        {
            std::lock_guard<std::mutex> lk(append_mutex);
            final_primes.insert(
                final_primes.end(), local_primes.begin(), local_primes.end()
            );
        }
    };

    std::vector<size_t> segment_indices;
    {
        for (size_t idx = 0;; ++idx) {
            ulong low = (ulong)idx * SEGMENT_SIZE;
            if (low > limit) break;
            segment_indices.push_back(idx);
        }
    }

    std::for_each(
        std::execution::par, segment_indices.begin(), segment_indices.end(),
        process_segment
    );

    std::sort(final_primes.begin(), final_primes.end());
    final_primes.erase(
        std::unique(final_primes.begin(), final_primes.end()),
        final_primes.end()
    );

    std::vector<ulong> result;
    result.reserve(final_primes.size());
    for (ulong p : final_primes)
        if (p >= 2 && p <= limit) result.push_back(p);

    return result;
}
