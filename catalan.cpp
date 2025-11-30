#include "sieve.hpp"

#include <cmath>
#include <cstring>
#include <execution>
#include <immintrin.h>
#include <map>
#include <vector>

typedef uint64_t ulong;

ulong legendre_factorial(ulong n, ulong p) {
    ulong exponent = 0, power = p;

    while (power <= n) {
        exponent += n / power;
        // Check for overflow before multiplying
        if (power > n / p) break;
        power *= p;
    }

    return exponent;
}

void compute_catalan_exponents_parallel(
    const std::vector<ulong> &primes, ulong n, std::vector<ulong> &exponents
) {
    size_t num_primes = primes.size();
    exponents.assign(num_primes, 0);

    std::vector<size_t> idx(num_primes);
    std::iota(idx.begin(), idx.end(), 0);

    std::for_each(
        std::execution::par_unseq, idx.begin(), idx.end(), [&](size_t i) {
            ulong p = primes[i], exp_2n = legendre_factorial(2 * n, p),
                  exp_n   = legendre_factorial(n, p),
                  exp_n_1 = legendre_factorial(n + 1, p);
            // v_p(Catalan_n) = v_p((2n)!) - v_p((n+1)!) - v_p(n!)
            exponents[i] = exp_2n - exp_n_1 - exp_n;
        }
    );
}

std::map<ulong, ulong> catalan_prime_factorization(ulong n) {
    std::vector<ulong> primes = segmented_wheel_sieve(n);

    std::vector<ulong> exponents;
    compute_catalan_exponents_parallel(primes, n, exponents);

    std::map<ulong, ulong> factorization;
    for (ulong i = 0; i < primes.size(); i++)
        if (exponents[i] > 0) factorization[primes[i]] = exponents[i];

    return factorization;
}
