/** @file string_format.h */

#ifndef STRING_FORMAT_H_
#define STRING_FORMAT_H_

#include <string>
#include <memory>

template <typename... Args>
std::string string_format(const std::string &format, Args... args) {
    // wrapper for std::snprintf which sets up buffer of required size

    // determine required buffer size :
    int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1;
    if (size_s <= 0) {
        throw std::runtime_error("Error during formatting.");
    }

    // allocate buffer :
    auto size = static_cast<size_t>(size_s);
    std::unique_ptr<char[]> buf(new char[size]);

    // write formatted string to buffer :
    std::snprintf(buf.get(), size, format.c_str(), args...);

    return std::string(buf.get(), buf.get() + size - 1);
}

#endif