// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace fmindex_collection {

struct BitStack {
    uint64_t size{0};
    uint64_t ones{0};
    std::vector<uint8_t> data;
    void push(bool bit) {
        if (size % 8 == 0) {
            data.push_back(0);
        }
        data.back() = data.back() | (bit << (size % 8));
        size += 1;
        if (bit) {
            ones += 1;
        }
    }

    bool value(size_t idx) const {
        auto block  = idx / 8;
        auto bitNbr = idx % 8;
        return (data[block] & (1ull << bitNbr));
    }

    void append(std::vector<uint8_t>& buffer) const {
        buffer.reserve(buffer.size() + 16 + data.size());
        {
            buffer.resize(buffer.size() + 8);
            memcpy(buffer.data() + buffer.size() - 8, &size, sizeof(size));
        }
        {
            buffer.resize(buffer.size() + 8);
            memcpy(buffer.data() + buffer.size() - 8, &ones, sizeof(ones));
        }

        for (auto e : data) {
            buffer.push_back(e);
        }
        if (data.size() % 8 > 0) {
            for (size_t i{data.size() % 8}; i < 8; ++i) {
                buffer.push_back(0);
            }
        }
    }


    size_t read(uint8_t const* buffer, [[maybe_unused]]size_t len) {
        assert(len >= 16);

        memcpy(&size, buffer, 8); buffer += 8;
        memcpy(&ones, buffer, 8); buffer += 8;

        data.resize(size / 8 + (((size % 8) > 0)?1:0));
        memcpy(data.data(), buffer, data.size());
        if (data.size() % 8 > 0) {
            return 16 + data.size() + (8-(data.size() % 8));
        }
        return 16 + data.size();
    }
};

}
