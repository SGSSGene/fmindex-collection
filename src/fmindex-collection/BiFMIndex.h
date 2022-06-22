#pragma once

#include "CSA.h"
#include "occtable/concepts.h"
#include "utils.h"

#include <algorithm>

namespace fmindex_collection {

template <OccTable Table, typename TCSA = CSA>
struct BiFMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    using TTable = Table;

    Table  occ;
    Table  occRev;
    TCSA   csa;

//private:
    BiFMIndex(std::vector<uint8_t> const& bwt, std::vector<uint8_t> const& bwtRev, TCSA _csa)
        : occ{bwt}
        , occRev{bwtRev}
        , csa{std::move(_csa)}
    {
        assert(bwt.size() == bwtRev.size());
        assert(occ.size() == occRev.size());
        if (bwt.size() != bwtRev.size()) {
            throw std::runtime_error("bwt don't have the same size: " + std::to_string(bwt.size()) + " " + std::to_string(bwtRev.size()));
        }
        if (occ.size() != occRev.size()) {
            throw std::runtime_error("occ don't have the same size: " + std::to_string(occ.size()) + " " + std::to_string(occRev.size()));
        }
        for (size_t sym{0}; sym < Sigma; ++sym) {
            if (occ.rank(occ.size(), sym) != occRev.rank(occ.size(), sym)) {
                throw std::runtime_error("wrong rank for the last entry");
            }
        }
        if constexpr (requires(Table t) {{ t.hasValue(size_t{}) }; }) {
            for (size_t i{0}; i < occ.size(); ++i) {
                if (csa.value(i).has_value()) {
                    occ.setValue(i);
                }
            }
        }
    }

public:
    /**!\brief Creates a BiFMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     */
    BiFMIndex(std::vector<std::vector<uint8_t>> _input, size_t samplingRate)
        : occ{cereal_tag{}}
        , occRev{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        auto [totalSize, inputText, inputSizes] = createSequences(_input);
        decltype(_input){}.swap(_input); // input memory can be deleted

        // create BurrowsWheelerTransform and CompressedSuffixArray
        auto [bwt, csa] = [&] () {
            auto sa  = createSA(inputText);
            auto bwt = createBWT(inputText, sa);
            auto csa = TCSA(std::move(sa), samplingRate, inputSizes);

            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        // create BurrowsWheelerTransform on reversed text
        auto bwtRev = [&]() {
            std::reverse(begin(inputText), end(inputText));
            auto saRev  = createSA(inputText);
            auto bwtRev = createBWT(inputText, saRev);
            return bwtRev;
        }();

        decltype(inputText){}.swap(inputText); // inputText memory can be deleted

        *this = BiFMIndex{bwt, bwtRev, std::move(csa)};
    }


    /*!\brief Specific c'tor for serialization use
     */
    BiFMIndex(cereal_tag)
        : occ{cereal_tag{}}
        , occRev{cereal_tag{}}
        , csa{cereal_tag{}}
    {}

    size_t memoryUsage() const requires OccTableMemoryUsage<Table> {
        return occ.memoryUsage() + occRev.memoryUsage() + csa.memoryUsage();
    }

    size_t size() const {
        return occ.size();
    }

    auto locate(size_t idx) const -> std::tuple<size_t, size_t> {
        if constexpr (requires(Table t) {{ t.hasValue(size_t{}) }; }) {
            bool v = occ.hasValue(idx);
            uint64_t steps{};
            while(!v) {
                idx = occ.rank_symbol(idx);
                steps += 1;
                v = occ.hasValue(idx);
            }
            auto [chr, pos] = csa.value(idx);
            return {chr, pos+steps};

        } else {
            auto opt = csa.value(idx);
            uint64_t steps{};
            while(!opt) {
                idx = [&]() {
                    if constexpr (requires(Table t) { { t.rank_symbol(size_t{}) }; }) {
                        return occ.rank_symbol(idx);
                    } else {
                        return occ.rank(idx, occ.symbol(idx));
                    }
                }();
                steps += 1;
                opt = csa.value(idx);
            }
            auto [chr, pos] = *opt;
            return {chr, pos+steps};
        }
    }

    auto locate(size_t idx, size_t maxSteps) const -> std::optional<std::tuple<size_t, size_t>> {
        auto opt = csa.value(idx);
        uint64_t steps{};
        for (;!opt and maxSteps > 0; --maxSteps) {
            idx = occ.rank(idx, occ.symbol(idx));
            steps += 1;
            opt = csa.value(idx);
        }
        if (opt) {
            std::get<1>(*opt) += steps;
        }
        return opt;
    }


    auto single_locate_step(size_t idx) const -> std::optional<std::tuple<size_t, size_t>> {
        return csa.value(idx);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, occRev, csa);
    }
};

}
