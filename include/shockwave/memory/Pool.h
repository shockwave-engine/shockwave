#pragma once

#include <cstddef>
#include <memory>

#include "shockwave/memory/allocator.h"
#include "shockwave/utility/assert.h"

namespace sw
{

///
/// @tparam T The type of the objects to be stored
/// @tparam Alloc The allocator to use for allocating objects
/// @tparam generation_bits The number of bits for the generation in the index
/// @tparam initial_size The initial number of objects in the pool
template<class T,
         class Alloc                 = std::allocator<T>,
         std::size_t generation_bits = 8,
         std::size_t initial_size    = 128>
class Pool
{
public:

    class Index
    {
        static constexpr std::size_t total_bits = 8 * sizeof(std::size_t);
        static constexpr std::size_t index_bits = total_bits - generation_bits;

        static_assert(generation_bits < total_bits,
                      "Must keep at least one bit for the index");

        static constexpr std::size_t bitmask(std::size_t n_bits) {
            return (1 << n_bits) - 1;
        }

        template<std::size_t gen_offset = index_bits>
        static constexpr std::size_t make_index(std::size_t index,
                                                std::size_t generation) {
            return (generation << gen_offset) | index;
        }

        static_assert(bitmask(2) == 0b011, "");
        static_assert(make_index<8>(0xF0, 0x0F) == 0x0FF0, "");

    public:

        Index(std::size_t index, std::size_t generation)
            : value(make_index(index, generation)) {
            ensure(index == (index & bitmask(index_bits)), "index overflow");
            ensure(generation == (generation & bitmask(generation_bits)),
                   "generation overflow");
        }

        std::size_t index() const { return value & bitmask(index_bits); }

        std::size_t generation() const {
            return (value >> index_bits) & bitmask(generation_bits);
        }

        bool operator==(const Index&) const  = default;
        auto operator<=>(const Index&) const = default;

    private:

        std::size_t value;
    };


private:
};

}  // namespace sw
