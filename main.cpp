#include "catalan.hpp"
#include "catalan_construction.hpp"
#include "sieve.hpp"

#include <iostream>
#include <map>

typedef uint64_t ulong;

template <typename T> void pretty_print_vector(std::vector<T> v) {
    auto len  = v.size();
    auto last = v[len - 1];
    if (v.size() > 100) {
        std::cout << "[";
        for (int i = 0; i < 50; i++) std::cout << v[i] << ", ";
        std::cout << "... ";
        for (int i = 0; i < 49; i++) std::cout << v[len - 50 + i] << ", ";
        std::cout << last << "]" << std::endl;
        return;
    }

    v.pop_back();
    std::cout << "[";
    for (auto p : v) std::cout << p << ", ";
    std::cout << last << "]";
    v.push_back(last);
}

void print_factorization(const std::map<ulong, ulong> &factorization, ulong n) {
    std::cout << "Prime factorization of Catalan number C(" << n << "):\n";
    std::cout << "C(" << n << ") = ";

    bool first = true;
    for (const auto &[prime, exponent] : factorization) {
        if (!first) std::cout << " * ";
        std::cout << prime;
        if (exponent > 1) { std::cout << "^" << exponent; }
        first = false;
    }
    std::cout << "\n\n";

    ulong total_prime_factors = 0;
    for (const auto &[prime, exponent] : factorization) {
        total_prime_factors += exponent;
    }

    std::cout << "Number of distinct primes: " << factorization.size() << "\n";
    std::cout << "Total prime factors (with multiplicity): "
              << total_prime_factors << "\n";
}

ulong simple_catalan_multiply(const std::map<ulong, ulong> &factorization) {
    ulong result = 1;

    for (const auto &[prime, exponent] : factorization) {
        ulong power = 1;
        for (ulong i = 0; i < exponent; i++) { power *= prime; }
        result *= power;
    }

    return result;
}

ulong catalan_direct(ulong n) {
    if (n == 0) return 1;

    ulong numerator   = 1;
    ulong denominator = 1;

    for (ulong k = 2; k <= n; k++) {
        numerator *= (n + k);
        denominator *= k;
    }

    return numerator / denominator;
}

int main() {
    // sieve test
    std::cout << "Hello First 20 million primes: ";
    auto primes = segmented_wheel_sieve(10'000'000);
    pretty_print_vector(primes);
    std::cout << "!" << std::endl;

    // factorization test
    for (ulong n = 0; n < 15; n++) {
        auto factorization = catalan_prime_factorization(n);
        print_factorization(factorization, n);

        if (n <= 15) {
            unsigned long direct = catalan_direct(n);
            unsigned long from_factorization =
                simple_catalan_multiply(factorization);

            std::cout << "Verification:\n";
            std::cout << "Direct computation: C(" << n << ") = " << direct
                      << "\n";
            std::cout << "From factorization: C(" << n
                      << ") = " << from_factorization << "\n";

            if (direct == from_factorization) {
                std::cout << "✓ Factorization verified!\n";
            } else {
                std::cout << "✗ Factorization error!\n";
            }
        }
    }
    std::cout << std::endl;

    // construction test
    ulong n  = 100;
    auto out = catalan_construct_2048(n);
    std::vector<ulong> v(std::begin(out), std::end(out));
    std::cout << "C(" << n << ") = " << limbs_to_decimal(v) << std::endl;

    return 0;
}
