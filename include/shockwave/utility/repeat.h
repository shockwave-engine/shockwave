#pragma once

#include <concepts>
#include <utility>
#include "shockwave/utility/types.h"

namespace sw::details
{

template<std::invocable<u32> Fn, u32 n, u32 i = 0>
struct repeat
{
    static constexpr void iterate(Fn&& fn) {
        fn(i);
        repeat<Fn, n, i + 1>::iterate(std::forward<Fn>(fn));
    }
};

template<std::invocable<u32> Fn, u32 n>
struct repeat<Fn, n, n>
{
    static constexpr void iterate(Fn&&) {}
};

}  // namespace sw::details

namespace sw
{

template<u32 n, std::invocable<u32> Fn>
void repeat(Fn&& fn) {
    details::repeat<Fn, n>::iterate(std::forward<Fn>(fn));
}

template<u32 n, u32 i, std::invocable<u32> Fn>
void repeat(Fn&& fn) {
    details::repeat<Fn, n, i>::iterate(std::forward<Fn>(fn));
}

}  // namespace sw
