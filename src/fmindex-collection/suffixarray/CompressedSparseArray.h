// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/Bitvector2L.h"
#include "../bitvector/SparseRBBitvector.h"
#include "../string/FlattenedBitvectors2L.h"
#include "concepts.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <tuple>

#if __has_include(<cereal/types/tuple.hpp>)
#include <cereal/types/tuple.hpp>
#endif
#if __has_include(<cereal/types/vector.hpp>)
#include <cereal/types/vector.hpp>
#endif


namespace fmc::suffixarray {

/** A compressed sparse array
 *
 * Space efficient array over sparse data
 */
struct CompressedSparseArray {
    using String = string::FlattenedBitvectors_512_64k<4>; // 0, 2, 4, or 8 bytes
    std::vector<uint16_t> documentContent2B; // 2 byte
    std::vector<uint32_t> documentContent4B; // 4 byte
    std::vector<uint64_t> documentContent8B; // 8 byte
    String                indicator;

    CompressedSparseArray() = default;
    CompressedSparseArray(CompressedSparseArray const&) = delete;
    CompressedSparseArray(CompressedSparseArray&&) noexcept = default;

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, std::optional<uint64_t>>
    CompressedSparseArray(range_t const& _range)
        : indicator{_range | std::views::transform([&](std::optional<uint64_t> const& t) -> uint8_t {
            if (!t) return 0;
            if (*t < std::numeric_limits<uint16_t>::max()) {
                documentContent2B.push_back(*t);
                return 1;
            }
            if (*t < std::numeric_limits<uint32_t>::max()) {
                documentContent4B.push_back(*t);
                return 2;
            }
            documentContent8B.push_back(*t);
            return 3;
        })}
        {
    }

    auto operator=(CompressedSparseArray const&) -> CompressedSparseArray& = delete;
    auto operator=(CompressedSparseArray&&) noexcept -> CompressedSparseArray& = default;

    auto value(size_t idx) const -> std::optional<uint64_t> {
        assert(idx < indicator.size());
        auto symb = indicator.symbol(idx);
        if (symb == 0) return std::nullopt;
        auto r = indicator.rank(idx, symb);
        if (symb == 1) return documentContent2B[r];
        else if (symb == 2) return documentContent4B[r];
        return documentContent8B[r];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(documentContent2B, documentContent4B, documentContent8B, indicator);
    }
};

/** A compressed sparse array v2
 *
 * Space efficient array over sparse data
 */
struct CompressedSparseArrayV2 {
    using Bitvector = bitvector::Bitvector2L<512, 65536>;
    using String = string::FlattenedBitvectors_512_64k<4>; // 1, 2, 4, or 8 bytes
    std::vector<uint16_t> documentContent1B; // 1 byte
    std::vector<uint16_t> documentContent2B; // 2 byte
    std::vector<uint32_t> documentContent4B; // 4 byte
    std::vector<uint64_t> documentContent8B; // 8 byte
    Bitvector             indicator;
    String                indicatorType;

    CompressedSparseArrayV2() = default;
    CompressedSparseArrayV2(CompressedSparseArrayV2 const&) = delete;
    CompressedSparseArrayV2(CompressedSparseArrayV2&&) noexcept = default;

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, std::optional<uint64_t>>
    CompressedSparseArrayV2(range_t const& _range)
        : indicator{_range | std::views::transform([&](std::optional<uint64_t> const& t) -> bool {
            return static_cast<bool>(t);
        })}
        , indicatorType{_range
            | std::views::filter([&](std::optional<uint64_t> const& t) -> bool {
                return t.has_value();
            })
            | std::views::transform([&](std::optional<uint64_t> const& t) -> uint8_t {
                if (*t < std::numeric_limits<uint8_t>::max()) {
                    documentContent1B.push_back(*t);
                    return 0;
                }
                if (*t < std::numeric_limits<uint16_t>::max()) {
                    documentContent2B.push_back(*t);
                    return 1;
                }
                if (*t < std::numeric_limits<uint32_t>::max()) {
                    documentContent4B.push_back(*t);
                    return 2;
                }
                documentContent8B.push_back(*t);
                return 3;
        })}
        {
    }

    auto operator=(CompressedSparseArrayV2 const&) -> CompressedSparseArrayV2& = delete;
    auto operator=(CompressedSparseArrayV2&&) noexcept -> CompressedSparseArrayV2& = default;

    auto value(size_t idx) const -> std::optional<uint64_t> {
        assert(idx < indicator.size());
        if (!indicator.symbol(idx)) return std::nullopt;

        auto pos = indicator.rank(idx);

        auto symb = indicatorType.symbol(pos);
        auto r    = indicatorType.rank(pos, symb);

        if (symb == 0) return documentContent1B[r];
        if (symb == 1) return documentContent2B[r];
        else if (symb == 2) return documentContent4B[r];
        return documentContent8B[r];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(documentContent1B, documentContent2B, documentContent4B, documentContent8B, indicator, indicatorType);
    }
};

/** A compressed sparse array v3
 *
 * Space efficient array over sparse data
 */
struct CompressedSparseArrayV3 {
    using Bitvector = bitvector::SparseRBBitvector<2, bitvector::Bitvector2L<512, 65536>, bitvector::Bitvector2L<512, 65536>>;
    using String = string::FlattenedBitvectors_512_64k<4>; // 1, 2, 4, or 8 bytes
    std::vector<uint16_t> documentContent1B; // 1 byte
    std::vector<uint16_t> documentContent2B; // 2 byte
    std::vector<uint32_t> documentContent4B; // 4 byte
    std::vector<uint64_t> documentContent8B; // 8 byte
    Bitvector             indicator;
    String                indicatorType;

    CompressedSparseArrayV3() = default;
    CompressedSparseArrayV3(CompressedSparseArrayV3 const&) = delete;
    CompressedSparseArrayV3(CompressedSparseArrayV3&&) noexcept = default;

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, std::optional<uint64_t>>
    CompressedSparseArrayV3(range_t const& _range)
        : indicator{_range | std::views::transform([&](std::optional<uint64_t> const& t) -> bool {
            return static_cast<bool>(t);
        })}
        , indicatorType{_range
            | std::views::filter([&](std::optional<uint64_t> const& t) -> bool {
                return t.has_value();
            })
            | std::views::transform([&](std::optional<uint64_t> const& t) -> uint8_t {
                if (*t < std::numeric_limits<uint8_t>::max()) {
                    documentContent1B.push_back(*t);
                    return 0;
                }
                if (*t < std::numeric_limits<uint16_t>::max()) {
                    documentContent2B.push_back(*t);
                    return 1;
                }
                if (*t < std::numeric_limits<uint32_t>::max()) {
                    documentContent4B.push_back(*t);
                    return 2;
                }
                documentContent8B.push_back(*t);
                return 3;
        })}
        {
    }

    auto operator=(CompressedSparseArrayV3 const&) -> CompressedSparseArrayV3& = delete;
    auto operator=(CompressedSparseArrayV3&&) noexcept -> CompressedSparseArrayV3& = default;

    auto value(size_t idx) const -> std::optional<uint64_t> {
        assert(idx < indicator.size());
        if (!indicator.symbol(idx)) return std::nullopt;

        auto pos = indicator.rank(idx);

        auto symb = indicatorType.symbol(pos);
        auto r    = indicatorType.rank(pos, symb);

        if (symb == 0) return documentContent1B[r];
        if (symb == 1) return documentContent2B[r];
        else if (symb == 2) return documentContent4B[r];
        return documentContent8B[r];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(documentContent1B, documentContent2B, documentContent4B, documentContent8B, indicator, indicatorType);
    }
};

/** A compressed sparse array v4
 *
 * Space efficient array over sparse data
 */
struct CompressedSparseArrayV4 {
    using Bitvector = bitvector::Bitvector2L<512, 65536>;
    using String = string::FlattenedBitvectors_512_64k<2>; // 4, or 8 bytes
    std::vector<uint32_t> documentContent4B; // 4 byte
    std::vector<uint64_t> documentContent8B; // 8 byte
    Bitvector             indicator;
    String                indicatorType;

    CompressedSparseArrayV4() = default;
    CompressedSparseArrayV4(CompressedSparseArrayV4 const&) = delete;
    CompressedSparseArrayV4(CompressedSparseArrayV4&&) noexcept = default;

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, std::optional<uint64_t>>
    CompressedSparseArrayV4(range_t const& _range)
        : indicator{_range | std::views::transform([&](std::optional<uint64_t> const& t) -> bool {
            return static_cast<bool>(t);
        })}
        , indicatorType{_range
            | std::views::filter([&](std::optional<uint64_t> const& t) -> bool {
                return t.has_value();
            })
            | std::views::transform([&](std::optional<uint64_t> const& t) -> bool {
                if (*t < std::numeric_limits<uint32_t>::max()) {
                    documentContent4B.push_back(*t);
                    return false;
                }
                documentContent8B.push_back(*t);
                return true;
        })}
        {
    }

    auto operator=(CompressedSparseArrayV4 const&) -> CompressedSparseArrayV4& = delete;
    auto operator=(CompressedSparseArrayV4&&) noexcept -> CompressedSparseArrayV4& = default;

    auto value(size_t idx) const -> std::optional<uint64_t> {
        assert(idx < indicator.size());
        if (!indicator.symbol(idx)) return std::nullopt;

        auto pos = indicator.rank(idx);

        auto symb = indicatorType.symbol(pos);
        auto r    = indicatorType.rank(pos, symb);

        if (!symb) return documentContent4B[r];
        return documentContent8B[r];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(documentContent4B, documentContent8B, indicator, indicatorType);
    }
};

#if 1
/** A compressed sparse array V5
 *
 * Space efficient array over sparse data
 */
struct CompressedSparseArrayV5 {
    using Bitvector = bitvector::SparseRBBitvector<2, bitvector::Bitvector2L<512, 65536>, bitvector::Bitvector2L<512, 65536>>;
    using String = string::FlattenedBitvectors_512_64k<2>; // 4, or 8 bytes
    std::vector<uint32_t> documentContent4B; // 4 byte
    std::vector<uint64_t> documentContent8B; // 8 byte
    Bitvector             indicator;
    String                indicatorType;

    CompressedSparseArrayV5() = default;
    CompressedSparseArrayV5(CompressedSparseArrayV5 const&) = delete;
    CompressedSparseArrayV5(CompressedSparseArrayV5&&) noexcept = default;

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, std::optional<uint64_t>>
    CompressedSparseArrayV5(range_t const& _range)
        : indicator{_range | std::views::transform([&](std::optional<uint64_t> const& t) -> bool {
            return static_cast<bool>(t);
        })}
        , indicatorType{_range
            | std::views::filter([&](std::optional<uint64_t> const& t) -> bool {
                return t.has_value();
            })
            | std::views::transform([&](std::optional<uint64_t> const& t) -> bool {
                if (*t < std::numeric_limits<uint32_t>::max()) {
                    documentContent4B.push_back(*t);
                    return false;
                }
                documentContent8B.push_back(*t);
                return true;
        })}
        {
    }

    auto operator=(CompressedSparseArrayV5 const&) -> CompressedSparseArrayV5& = delete;
    auto operator=(CompressedSparseArrayV5&&) noexcept -> CompressedSparseArrayV5& = default;

    auto value(size_t idx) const -> std::optional<uint64_t> {
        assert(idx < indicator.size());
        if (!indicator.symbol(idx)) return std::nullopt;

        auto pos = indicator.rank(idx);

        auto symb = indicatorType.symbol(pos);
        auto r    = indicatorType.rank(pos, symb);

        if (!symb) return documentContent4B[r];
        return documentContent8B[r];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(documentContent4B, documentContent8B, indicator, indicatorType);
    }
};
#endif

}
