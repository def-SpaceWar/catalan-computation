#ifndef catalan_construction_hpp_INCLUDED
#define catalan_construction_hpp_INCLUDED

#include <array>
#include <cstdint>
#include <string>
#include <vector>

typedef uint64_t ulong;
std::vector<ulong> catalan_construct_kbit(ulong n, size_t k_bits);
std::array<ulong, 4> catalan_construct_256(ulong n);
std::array<ulong, 8> catalan_construct_512(ulong n);
std::array<ulong, 16> catalan_construct_1024(ulong n);
std::array<ulong, 32> catalan_construct_2048(ulong n);
std::string limbs_to_decimal(const std::vector<ulong> &limbs);

#endif // catalan_construction_hpp_INCLUDED
