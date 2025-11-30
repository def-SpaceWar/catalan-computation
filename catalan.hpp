#ifndef catalan_hpp_INCLUDED
#define catalan_hpp_INCLUDED

#include <cstdint>
#include <map>

typedef uint64_t ulong;
std::map<ulong, ulong> catalan_prime_factorization(ulong n);

#endif // catalan_hpp_INCLUDED
