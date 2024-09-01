#pragma once

#include <cstddef>
#include <cstdlib>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>
#include <queue>
#include <vector>
#include <ranges>

#include "shockwave/memory/allocator.h"
#include "shockwave/utility/assert.h"

namespace sw
{

/// A pool of objects, constructors and destructors are only called
/// once, regardless of how many times they are taken and given back.
///
/// Unless generation_bits == 0, the pool will keep track of the generations of
/// objects. This avoids pointing to a new object after the old one was given
/// back.
///
/// The pool will grow automatically, it cannot be shrinked unless reassigned.
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
    public:

        static constexpr std::size_t total_bits = 8 * sizeof(std::size_t);
        static constexpr std::size_t index_bits = total_bits - generation_bits;
        static constexpr std::size_t gen_bits   = generation_bits;

        static_assert(generation_bits < total_bits,
                      "Must keep at least one bit to encode existence");

        Index(std::size_t index, std::size_t generation)
            : value(make_index(index, generation)) {
            ensure(index == (index & bitmask<index_bits>()), "index overflow");
        }

        std::size_t index() const { return value & bitmask<index_bits>(); }

        std::size_t generation() const {
            if constexpr ( gen_bits == 0 ) {
                return 0;  // avoid shift by type width warning
            } else {
                return (value >> index_bits) & bitmask<generation_bits>();
            }
        }

        bool operator==(const Index&) const  = default;
        auto operator<=>(const Index&) const = default;

        friend std::ostream& operator<<(std::ostream& os, Index index) {
            return os << "{ index = " << index.index()
                      << ", generation = " << index.generation() << " }";
        }

    private:

        template<std::size_t n_bits>
        static constexpr std::size_t bitmask() {
            if constexpr ( n_bits == 0 ) {
                return 0;
            } else {
                return static_cast<std::size_t>(1) << (n_bits - 1)
                     | bitmask<n_bits - 1>();
            }
        }

        template<std::size_t gen_offset = index_bits>
        static constexpr std::size_t make_index(std::size_t index,
                                                std::size_t generation) {
            if constexpr ( gen_offset == total_bits ) {
                return index;  // avoid shift by width of type warning
            } else {
                return (generation << gen_offset) | index;
            }
        }

        static_assert(bitmask<2>() == 0b011, "");
        static_assert(make_index<8>(0xF0, 0x0F) == 0x0FF0, "");

        std::size_t value;
    };


    Pool(const Alloc& alloc = Alloc{})
        : pool{alloc},
          free_list{alloc},
          generations{alloc} {
        reallocate();
    }

    Index take() {
        if ( free_list.empty() ) {
            reallocate();
        }

        Index index = free_list.top();

        free_list.pop();
        generations[index.index()] = Generation{index.generation()};

        return index;
    }

    void give_back(Index index) {
        generations[index.index()] = Generation();
        free_list.emplace(index.index(), index.generation() + 1);
    }

    std::optional<const T*> get(Index index) const {
        if ( generations[index.index()] == index.generation() ) {
            return &pool[index.index()];
        } else {
            return std::nullopt;
        }
    }

    std::optional<T*> get(Index index) {
        if ( generations[index.index()] == index.generation() ) {
            return &pool[index.index()];
        } else {
            return std::nullopt;
        }
    }

private:

    void reallocate() {
        ensure(free_list.empty(),
               "Should not reallocate until all available indices are used");

        std::size_t before = pool.size();
        std::size_t new_size
            = pool.empty() ? initial_size : 2 * pool.capacity();

        pool.resize(new_size);
        generations.resize(new_size);

        for ( std::size_t i : std::views::iota(before, pool.size()) ) {
            free_list.push(Index{i, 0});
        }
    }

    struct Generation
    {
        static constexpr std::size_t not_set = 1 << generation_bits;

        bool operator==(std::size_t generation) const {
            return value == generation;
        }

        std::size_t value = not_set;
    };


    std::vector<T, Alloc> pool;

    std::priority_queue<Index,
                        std::vector<Index, sw::rebind_to<Alloc, Index>>,
                        std::greater<>>
        free_list;

    std::vector<Generation, rebind_to<Alloc, Generation>> generations;
};

}  // namespace sw
