#pragma once

#include "CSA.h"
#include "occtable/concepts.h"

#include <cassert>
#include <divsufsort64.h>

namespace fmindex_collection {

template <OccTable Table>
struct BiFMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    Table occ;
    Table occRev;
    CSA   csa;

    BiFMIndex(std::vector<uint8_t> const& bwt, std::vector<uint8_t> const& bwtRev, CSA _csa)
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
    }
//    BiFMIndex(BiFMIndex&&) noexcept = default;

    static auto createSA(std::vector<uint8_t> const& input) -> std::vector<int64_t> {
        auto sa = std::vector<int64_t>{};
        sa.resize(input.size());
        auto error = divsufsort64(static_cast<uint8_t const*>(input.data()), sa.data(), input.size());
        if (error != 0) {
            throw std::runtime_error("some error while creating the suffix array");
        }
        return sa;
    }

    static auto createBWT(std::vector<uint8_t> const& input, std::vector<int64_t> const& sa) -> std::vector<uint8_t> {
        assert(input.size() == sa.size());
        std::vector<uint8_t> bwt;
        bwt.resize(input.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            bwt[i] = input[(sa[i] + input.size()- 1) % input.size()];
        }
        return bwt;
    }

    static auto createCSA(std::vector<int64_t> sa, size_t samplingRate) -> CSA {
        auto bitStack = fmindex_collection::BitStack{};
        auto ssa      = std::vector<uint64_t>{};
        if (samplingRate > 0) {
            ssa.reserve(sa.size() / samplingRate + 1);
        }
        for (size_t i{0}; i < sa.size(); ++i) {
            bool sample = samplingRate == 0 || i % samplingRate == 0;
            bitStack.push(sample);
            if (sample) {
                ssa.push_back(sa[i]);
            }
        }
        return CSA{std::move(ssa), bitStack};
    }


    BiFMIndex(std::vector<uint8_t> input, size_t samplingRate)
        : occ{cereal_tag{}}
        , occRev{cereal_tag{}}
        , csa{cereal_tag{}}
    {

        auto [bwt, csa] = [&] () {
            auto sa  = createSA(input);
            auto bwt = createBWT(input, sa);
            auto csa = createCSA(std::move(sa), samplingRate);

            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        auto bwtRev = [&]() {
            std::reverse(begin(input), end(input));
            auto saRev  = createSA(input);
            auto bwtRev = createBWT(input, saRev);
            return bwtRev;
        }();

        decltype(input){}.swap(input); // input memory can be deleted

        *this = BiFMIndex{bwt, bwtRev, std::move(csa)};
    }

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

    size_t locate(size_t idx) const {
        auto opt = csa.value(idx);
        uint64_t steps{};
        while(!opt) {
            idx = occ.rank(idx, occ.symbol(idx));
            steps += 1;
            opt = csa.value(idx);
        }
        if (opt.value() + steps >= size()) {
            return opt.value() + steps - size();
        }
        return opt.value() + steps;
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, occRev, csa);
    }
};

}
