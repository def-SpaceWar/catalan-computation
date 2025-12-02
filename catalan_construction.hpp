#ifndef catalan_construction_hpp_INCLUDED
#define catalan_construction_hpp_INCLUDED

#include <array>
#include <cstdint>
#include <string>
#include <vector>

typedef uint64_t ulong;
std::vector<ulong> catalan_construct_kbit(ulong n, size_t k_bits);
std::string limbs_to_decimal(const std::vector<ulong> &limbs);

#endif // catalan_construction_hpp_INCLUDED
