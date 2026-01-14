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
template <typename T, typename T2=std::nullptr_t>
struct RestoreAdd {
    T* value;
    T oldValue;
    RestoreAdd(T& _value, T2 newValue)
        : value{&_value}
        , oldValue{*value}
    {
        *value += static_cast<T>(newValue);
    }
    RestoreAdd(RestoreAdd const&) = delete;
    RestoreAdd(RestoreAdd&&) = delete;
    ~RestoreAdd() {
        *value = oldValue;
    }
};

template <typename T, typename T2=std::nullptr_t>
struct RestoreSub {
    T* value;
    T oldValue;
    RestoreSub(T& _value, T2 newValue)
        : value{&_value}
        , oldValue{*value}
    {
        *value -= static_cast<T>(newValue);
    }
    RestoreSub(RestoreSub const&) = delete;
    RestoreSub(RestoreSub&&) = delete;
    ~RestoreSub() {
        *value = oldValue;
    }
};
