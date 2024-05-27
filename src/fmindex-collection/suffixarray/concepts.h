// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/concepts.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <tuple>
#include <utility>

namespace fmindex_collection {

template <typename T>
concept SuffixArray_c = requires(T t, T t2, size_t idx) {
    /** Default constructible **/
    { T{} } -> std::same_as<T>;

    /** Moveable */
    { T{std::move(t)} } -> std::same_as<T>;

    /** Move assignable */
    { t = T{} } -> std::same_as<T&>;

    /* Return a value at a certain position
     *
     * \param idx - row index
     * \return - optional value of tuple, first value is the sequence number, the second the position inside the sequence
     */
    { t.value(idx) } -> std::same_as<std::optional<std::tuple<uint64_t, uint64_t>>>;

    /* Push another sampling suffix array entry
     *
     * \param seqNr, pos
     */
    { t.push_back(std::optional<std::tuple<size_t, size_t>>{}) };
};

}
