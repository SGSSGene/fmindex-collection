#pragma once

#include "CSA.h"
#include "occtable/concepts.h"
#include "utils.h"

namespace fmindex_collection {

template <OccTable Table>
struct FMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    using TTable = Table;

    Table  occ;
    CSA    csa;


    FMIndex(std::vector<uint8_t> const& bwt, CSA _csa)
        : occ{bwt}
        , csa{std::move(_csa)}
    {}

    FMIndex(std::vector<uint8_t> _input, size_t samplingRate)
        : occ{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        auto input = std::vector<std::vector<uint8_t>>{std::move(_input)};
        *this = FMIndex{std::move(input), samplingRate};
    }

    FMIndex(std::vector<std::vector<uint8_t>> _input, size_t samplingRate)
        : occ{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        size_t totalSize = std::accumulate(begin(_input), end(_input), size_t{0}, [](auto s, auto const& l) { return s + l.size() + 1; });

        auto input = std::vector<uint8_t>{};
        input.reserve(totalSize);

        auto inputSizes = std::vector<size_t>{};
        inputSizes.reserve(_input.size());

        for (auto const& l : _input) {
            input.insert(end(input), begin(l), end(l));
            input.emplace_back(0);
            inputSizes.emplace_back(l.size()+1);
        }
        decltype(_input){}.swap(_input); // input memory can be deleted

        auto [bwt, csa] = [&input, &samplingRate, &inputSizes, this] () {
            auto sa  = createSA(input);
            auto bwt = createBWT(input, sa);
            auto csa = CSA{std::move(sa), samplingRate, inputSizes};

            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        decltype(input){}.swap(input); // input memory can be deleted

        *this = FMIndex{bwt, std::move(csa)};
    }


    FMIndex(cereal_tag)
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
        auto opt = csa.value(idx);
        uint64_t steps{};
        while(!opt) {
            idx = occ.rank(idx, occ.symbol(idx));
            steps += 1;
            opt = csa.value(idx);
        }
        auto [chr, pos] = *opt;
        return std::make_tuple(chr, pos + steps);
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
        ar(occ, csa);
    }
};

}
