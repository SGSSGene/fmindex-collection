// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cstddef>
#include <cstdint>
#include <mmser/mmser.h>

namespace fmc {

struct VectorBool {
    struct Proxy {
        uint8_t& value;
        size_t   offset;
        operator bool() const {
            return (value >> offset) & 1;
        }
        void operator=(bool b) {
            //       value     clears value at offset |   sets value at offset
            value = (value & ~(uint8_t{1} << offset)) | (uint8_t{b} << offset);
        }
    };

    mmser::vector<uint8_t> values;
    size_t used_capacity{};


    void reserve(size_t i) {
        values.reserve((i+7)/8);
    }
    void resize(size_t i) {
        values.resize((i+7)/8, 0);
        used_capacity = i;
    }

    auto unused_capacity() const {
        return values.size()*8 - used_capacity;
    }
    void resize(size_t newSize, bool value) {
        values.resize((newSize+7)/8, 0);
        for (size_t i{used_capacity}; i < newSize; ++i) {
            at(i) = value;
        }
        used_capacity = newSize;
    }

    void push_back(bool b) {
        if (unused_capacity() == 0) {
            values.emplace_back(0);
        }
        at(size()) = b;
        used_capacity += 1;
    }

    auto size() const -> size_t {
        return used_capacity;
    }

    auto at(size_t idx) -> Proxy {
        size_t block = idx/8;
        size_t offset = idx%8;
        return Proxy{values[block], offset};
    }
    auto at(size_t idx) const -> bool {
        size_t block = idx/8;
        size_t offset = idx%8;
        return (values[block] >> offset) & 1;
    }
    auto operator[](this auto& self, size_t idx) {
        return self.at(idx);
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.used_capacity, self.values);
    }
};

}
