#pragma once

#include <cassert>

// May be changed to an inline function throwing an exception with stacktrace
#define ensure(cond, msg) assert(cond&& msg)
