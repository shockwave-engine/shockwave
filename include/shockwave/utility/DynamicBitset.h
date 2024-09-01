#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace sw
{

template<class Alloc = std::allocator<uint8_t>>
class DynamicBitset
{
public:

    DynamicBitset(const Alloc& alloc = Alloc{}) : bytes{alloc} {}

    std::size_t size() const { return bytes.size() * 8; }
    void        resize(std::size_t n_bits) { bytes.resize(n_bytes(n_bits)); }

    bool get(std::size_t bit) const {
        Index index = get_index(bit);
        return bytes[index.byte] & (1 << index.bit);
    }

    void set(std::size_t bit) {
        Index index        = get_index(bit);
        bytes[index.byte] |= (1 << index.bit);
    }

    void unset(std::size_t bit) {
        Index index        = get_index(bit);
        bytes[index.byte] &= ~(1 << index.bit);
    }

private:

    struct Index
    {
        std::size_t byte;
        std::size_t bit;
    };

    static constexpr Index get_index(std::size_t bit) {
        return Index{
            .byte = bit / 8,
            .bit  = bit % 8,
        };
    }

    static constexpr std::size_t n_bytes(std::size_t n_bits) {
        Index index = get_index(n_bits);
        return index.bit == 0 ? index.byte : index.byte + 1;
    }

    std::vector<uint8_t, Alloc> bytes;
};

}  // namespace sw
