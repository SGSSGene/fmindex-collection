// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cstddef>

template <typename T, typename T2=std::nullptr_t>
struct Restore {
    T* value;
    T oldValue;
    Restore(T& _value)
        : value{&_value}
        , oldValue{*value}
    {}
    Restore(T& _value, T2 newValue)
        : value{&_value}
        , oldValue{*value}
    {
        *value = static_cast<T>(newValue);
    }
    Restore(Restore const&) = delete;
    Restore(Restore&&) = delete;
    ~Restore() {
        *value = oldValue;
    }
};
