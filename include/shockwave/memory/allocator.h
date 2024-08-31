#pragma once

#include <memory>

namespace sw
{

template<class Alloc>
struct rebind
{
    template<class T>
    using to = typename std::allocator_traits<Alloc>::template rebind_alloc<T>;
};

}  // namespace sw

static_assert(std::is_same_v<std::allocator<int>,
                             sw::rebind<std::allocator<float>>::to<int>>);
