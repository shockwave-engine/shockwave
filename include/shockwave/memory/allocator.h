#pragma once

#include <memory>

namespace sw
{

template<class Alloc, class T>
using rebind_to =
    typename std::allocator_traits<Alloc>::template rebind_alloc<T>;

}  // namespace sw

static_assert(std::is_same_v<std::allocator<int>,
                             sw::rebind_to<std::allocator<float>, int>>);
