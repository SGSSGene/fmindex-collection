#pragma once

#include "CSA.h"
#include "occtable/concepts.h"
#include "utils.h"

namespace fmindex_collection {

/*
 * Same as the FMIndex, but build over the reverse text.
 *
 * TODO: has some correctional features...
 */
template <OccTable Table, typename TCSA = CSA>
struct ReverseFMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    using TTable = Table;

    Table  occ;
    TCSA   csa;

    ReverseFMIndex(std::vector<uint8_t> const& bwt, TCSA _csa)
        : occ{bwt}
        , csa{std::move(_csa)}
    {}

    ReverseFMIndex(Sequences auto const& _input, size_t samplingRate)
        : occ{cereal_tag{}}
        , csa{cereal_tag{}}
    {

        auto [totalSize, inputText, inputSizes] = createSequences(_input, true);

        auto [bwt, csa] = [&, &inputText=inputText, &inputSizes=inputSizes] () {
            auto sa  = createSA(inputText);
            auto bwt = createBWT(inputText, sa);
            auto csa = TCSA{std::move(sa), samplingRate, inputSizes, true};

            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        decltype(inputText){}.swap(inputText); // inputText memory can be deleted

        *this = ReverseFMIndex{bwt, std::move(csa)};
    }


    ReverseFMIndex(cereal_tag)
        : occ{cereal_tag{}}
        , csa{cereal_tag{}}
    {}

    size_t memoryUsage() const requires OccTableMemoryUsage<Table> {
        return occ.memoryUsage() + csa.memoryUsage();
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
            return {chr, pos-steps};

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
            return {chr, pos-steps};
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
            std::get<1>(*opt) -= steps;
        }
        return opt;
    }


    auto single_locate_step(size_t idx) const -> std::optional<std::tuple<size_t, size_t>> {
        return csa.value(idx);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, csa);
    }
};

}
