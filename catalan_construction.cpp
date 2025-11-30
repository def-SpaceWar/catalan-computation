#include "catalan_construction.hpp"
#include "catalan.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <gmp.h>
#include <map>
#include <set>
#include <string>
#include <sys/types.h>
#include <vector>

#ifdef _WIN32
    #include <windows.h>
#else
    #include <sys/sysinfo.h>
    #include <unistd.h>
#endif

typedef uint64_t ulong;

// Parameters for memory budgeting (tweakable)
static constexpr double MEM_FRACTION         = 0.85; // f
static constexpr unsigned CONCURRENCY_FACTOR = 4,    // c
    MEMORY_MULTIPLICATION_OVERHEAD           = 4,    // F
    MIN_TARGET_BITS = 1'024; // don't allow target_bits too small

static ulong get_free_ram_bytes() {
#ifdef _WIN32
    MEMORYSTATUSEX st;
    st.dwLength = sizeof(st);
    if (!GlobalMemoryStatusEx(&st)) return 512ULL * 1'024ULL * 1'024ULL;
    return (ulong)st.ullAvailPhys;
#else
    struct sysinfo si;
    if (sysinfo(&si) == 0) {
        ulong free_bytes = (ulong)si.freeram * (ulong)si.mem_unit;
        return free_bytes;
    } else {
        return 512ULL * 1'024ULL * 1'024ULL;
    }
#endif
}

static mpz_t* mpz_new_one() {
    mpz_t* z = (mpz_t*)std::malloc(sizeof(mpz_t));
    if (!z) std::abort();
    mpz_init_set_ui(*z, 1);
    return z;
}

static mpz_t* mpz_new_from_ull(ulong v) {
    mpz_t* z = (mpz_t*)std::malloc(sizeof(mpz_t));
    if (!z) std::abort();
    mpz_init(*z);
    mpz_import(*z, 1, -1, sizeof(ulong), 0, 0, &v);
    return z;
}

static void mpz_free(mpz_t* z) {
    if (!z) return;
    mpz_clear(*z);
    std::free(z);
}

static inline size_t mpz_bitlen(const mpz_t x) {
    if (mpz_cmp_ui(x, 0) == 0) return 0;
    return mpz_sizeinbase(x, 2);
}

static mpz_t* product_of_ulongs_balanced(const std::vector<ulong> &primes) {
    if (primes.empty()) return mpz_new_one();

    struct Node {
        size_t bits;
        mpz_t* v;
    };

    struct Cmp {
        bool operator()(Node const &a, Node const &b) const {
            if (a.bits != b.bits) return a.bits < b.bits;
            return a.v < b.v;
        }
    };

    std::multiset<Node, Cmp> pool;
    pool.clear();

    for (ulong p : primes) {
        mpz_t* z = mpz_new_from_ull(p);
        pool.insert(Node{mpz_bitlen(*z), z});
    }

    while (pool.size() > 1) {
        auto it1 = pool.begin();
        Node n1  = *it1;
        pool.erase(it1);
        auto it2 = pool.begin();
        Node n2  = *it2;
        pool.erase(it2);

        mpz_t* prod = (mpz_t*)std::malloc(sizeof(mpz_t));
        if (!prod) std::abort();
        mpz_init(*prod);
        mpz_mul(*prod, *(n1.v), *(n2.v));

        mpz_free(n1.v);
        mpz_free(n2.v);

        pool.insert(Node{mpz_bitlen(*prod), prod});
    }

    if (pool.empty()) return mpz_new_one();
    mpz_t* result = pool.begin()->v;
    pool.erase(pool.begin());
    return result;
}

static mpz_t*
balanced_product_tree_components(std::vector<mpz_t*> &components) {
    struct Node {
        size_t bits;
        mpz_t* v;
    };

    struct Cmp {
        bool operator()(Node const &a, Node const &b) const {
            if (a.bits != b.bits) return a.bits < b.bits;
            return a.v < b.v;
        }
    };

    std::multiset<Node, Cmp> pool;
    for (mpz_t* p : components) pool.insert(Node{mpz_bitlen(*p), p});
    components.clear();

    while (pool.size() > 1) {
        auto it1 = pool.begin();
        Node n1  = *it1;
        pool.erase(it1);
        auto it2 = pool.begin();
        Node n2  = *it2;
        pool.erase(it2);

        mpz_t* prod = (mpz_t*)std::malloc(sizeof(mpz_t));
        if (!prod) std::abort();
        mpz_init(*prod);
        mpz_mul(*prod, *(n1.v), *(n2.v));

        mpz_free(n1.v);
        mpz_free(n2.v);

        pool.insert(Node{mpz_bitlen(*prod), prod});
    }

    if (pool.empty()) return mpz_new_one();
    mpz_t* result = pool.begin()->v;
    pool.erase(pool.begin());
    return result;
}

static inline ulong ceil_div_u64(ulong a, ulong b) { return (a + b - 1) / b; }

static std::vector<mpz_t*> build_components_from_factorization(
    const std::vector<std::pair<ulong, ulong>> &factors, ulong target_bits
) {
    std::map<ulong, std::vector<ulong>> groups;
    for (auto &pe : factors) {
        ulong p = pe.first, e = pe.second;
        if (p == 2) continue;
        groups[e].push_back(p);
    }
    std::vector<mpz_t*> components;

    for (auto &kv : groups) {
        ulong e                    = kv.first;
        std::vector<ulong> &primes = kv.second;
        if (primes.empty()) continue;

        ulong bits_per_prime = 0;
        for (ulong p : primes) {
            bits_per_prime += (ulong)std::floor(std::log2((long double)p)) + 1;
        }
        bits_per_prime = std::max<ulong>(1, bits_per_prime / primes.size());

        ulong K          = (ulong)primes.size();
        ulong total_bits = (ulong)e * K * bits_per_prime;
        ulong N_e        = ceil_div_u64(total_bits, target_bits);
        if (N_e == 0) N_e = 1;
        ulong chunk_size = ceil_div_u64(K, N_e);

        for (size_t pos = 0; pos < primes.size(); pos += chunk_size) {
            size_t end = std::min(primes.size(), pos + chunk_size);
            std::vector<ulong> chunk_primes;
            chunk_primes.reserve(end - pos);
            for (size_t i = pos; i < end; ++i)
                chunk_primes.push_back(primes[i]);

            mpz_t* chunk_prod = product_of_ulongs_balanced(chunk_primes);
            if (e > 1) {
                mpz_t* powed = (mpz_t*)std::malloc(sizeof(mpz_t));
                if (!powed) std::abort();
                mpz_init(*powed);
                mpz_pow_ui(*powed, *chunk_prod, (unsigned long)e);
                mpz_free(chunk_prod);
                chunk_prod = powed;
            }
            components.push_back(chunk_prod);
        }
    }

    if (!components.size()) { components.push_back(mpz_new_one()); }
    return components;
}

static std::vector<ulong> catalan_construct_kbit(ulong n, size_t k_bits) {
    auto factors = catalan_prime_factorization(n);

    ulong free_ram = get_free_ram_bytes();
    double A = (double)free_ram, f = MEM_FRACTION,
           c          = (double)CONCURRENCY_FACTOR,
           F          = (double)MEMORY_MULTIPLICATION_OVERHEAD;
    ulong target_bits = (ulong)std::floor((A * f * 8.0) / (c * F));
    if (target_bits < MIN_TARGET_BITS) target_bits = MIN_TARGET_BITS;

    ulong two_exp = 0;
    for (auto &pe : factors)
        if (pe.first == 2) two_exp = pe.second;
    std::vector<std::pair<ulong, ulong>> factors_no2;
    factors_no2.reserve(factors.size());
    for (auto &pe : factors)
        if (pe.first != 2) factors_no2.push_back(pe);

    auto components =
        build_components_from_factorization(factors_no2, target_bits);
    mpz_t* final_prod = balanced_product_tree_components(components);

    if (two_exp > 0)
        mpz_mul_2exp(*final_prod, *final_prod, (mp_bitcnt_t)two_exp);

    mpz_t r;
    mpz_init(r);
    if (k_bits > 0) {
        mpz_tdiv_r_2exp(r, *final_prod, (mp_bitcnt_t)k_bits);
    } else {
        mpz_set(r, *final_prod);
    }

    size_t limbs = (k_bits + 63) / 64;
    if (limbs == 0) limbs = 1;
    std::vector<ulong> out(limbs, 0);
    size_t written = 0;
    mpz_export(out.data(), &written, -1, sizeof(ulong), 0, 0, r);

    mpz_clear(r);
    mpz_free(final_prod);
    return out;
}

std::array<ulong, 4> catalan_construct_256(ulong n) {
    auto v = catalan_construct_kbit(n, 256);
    std::array<ulong, 4> out{};
    for (size_t i = 0; i < out.size(); ++i)
        out[i] = (i < v.size() ? v[i] : 0ULL);
    return out;
}

std::array<ulong, 8> catalan_construct_512(ulong n) {
    auto v = catalan_construct_kbit(n, 512);
    std::array<ulong, 8> out{};
    for (size_t i = 0; i < out.size(); ++i)
        out[i] = (i < v.size() ? v[i] : 0ULL);
    return out;
}

std::array<ulong, 16> catalan_construct_1024(ulong n) {
    auto v = catalan_construct_kbit(n, 1'024);
    std::array<ulong, 16> out{};
    for (size_t i = 0; i < out.size(); ++i)
        out[i] = (i < v.size() ? v[i] : 0ULL);
    return out;
}

std::array<ulong, 32> catalan_construct_2048(ulong n) {
    auto v = catalan_construct_kbit(n, 2'048);
    std::array<ulong, 32> out{};
    for (size_t i = 0; i < out.size(); ++i)
        out[i] = (i < v.size() ? v[i] : 0ULL);
    return out;
}

std::string limbs_to_decimal(const std::vector<ulong> &limbs) {
    std::vector<unsigned __int128> a(limbs.size());
    for (size_t i = 0; i < limbs.size(); ++i)
        a[i] = (unsigned __int128)limbs[i];

    if (std::all_of(a.begin(), a.end(), [](unsigned __int128 v) {
            return v == 0;
        }))
        return "0";

    std::string out;
    out.reserve(80);
    while (true) {
        unsigned __int128 carry = 0;
        for (int i = (int)a.size() - 1; i >= 0; --i) {
            unsigned __int128 cur = (carry << 64) | a[i];
            a[i]                  = cur / 10;
            carry                 = cur % 10;
        }
        out.push_back(char('0' + (unsigned)carry));
        bool all_zero = true;
        for (auto &v : a)
            if (v != 0) {
                all_zero = false;
                break;
            }
        if (all_zero) break;
    }
    std::reverse(out.begin(), out.end());
    return out;
}
