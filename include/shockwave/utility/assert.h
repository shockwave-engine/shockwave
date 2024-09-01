#pragma once

#include <exception>
#include <source_location>
#include <sstream>

#define ensure(cond, msg) \
    if ( !cond )          \
    throw ::sw::Exception(msg)

namespace sw
{

class Exception : public std::exception
{
public:

    Exception(const std::string&   msg,
              std::source_location loc = std::source_location::current())
        : std::exception{} {
        std::stringstream ss{};
        ss << "From " << loc.function_name() << " in " << loc.file_name() << ":"
           << loc.line() << '\n'
           << msg;
        msg_ = ss.str();
    }

    char const* what() const noexcept override { return msg_.c_str(); }

private:

    std::string msg_;
};

}  // namespace sw
