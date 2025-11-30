#ifndef sieve_hpp_INCLUDED
#define sieve_hpp_INCLUDED

#include <cstdint>
#include <vector>

typedef uint64_t ulong;
std::vector<ulong> segmented_wheel_sieve(ulong n);

#endif // sieve_hpp_INCLUDED
