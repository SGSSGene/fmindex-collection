#pragma once

#include "CSA.h"
#include "occtable/concepts.h"

#include <cassert>
#include <sdsl/divsufsort.hpp>
#include <numeric>

namespace fmindex_collection {

template <OccTable Table>
struct BiFMIndex {
    static size_t constexpr Sigma = Table::Sigma;
    static size_t constexpr BitsForPosition = 48;

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
        auto error = sdsl::divsufsort<int64_t>(static_cast<uint8_t const*>(input.data()), sa.data(), input.size());
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

    static auto createCSA(std::vector<int64_t> sa, size_t samplingRate, std::vector<size_t> const& _inputSizes) -> CSA {
        auto bitStack = fmindex_collection::BitStack{};
        auto ssa      = std::vector<uint64_t>{};
        if (samplingRate > 0) {
            ssa.reserve(sa.size() / samplingRate + 1);
        }
        for (size_t i{0}; i < sa.size(); ++i) {
            bool sample = (samplingRate == 0) || (sa[i] % samplingRate == 0);
            if (!sample) {
                size_t acc{};
                for (auto l : _inputSizes) {
                    if (static_cast<size_t>(sa[i]) == acc) {
                        sample = true;
                        break;
                    }
                    if (static_cast<size_t>(sa[i]) < acc) {
                        break;
                    }
                    acc += l;
                }
            }
            bitStack.push(sample);
            if (sample) {
                ssa.push_back(sa[i]);
            }
        }
        for (auto& s : ssa) {
            // find which reference id this is
            for (size_t i{0}; i < _inputSizes.size(); ++i) {
                if (s < _inputSizes[i]) {
                    s = s | (i << BitsForPosition);
                    break;
                }
                s -= _inputSizes[i];
            }
        }
        return CSA{std::move(ssa), bitStack};
    }

    BiFMIndex(std::vector<uint8_t> _input, size_t samplingRate)
        : occ{cereal_tag{}}
        , occRev{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        auto input = std::vector<std::vector<uint8_t>>{std::move(_input)};
        *this = BiFMIndex{std::move(input), samplingRate};
    }

    BiFMIndex(std::vector<std::vector<uint8_t>> _input, size_t samplingRate)
        : occ{cereal_tag{}}
        , occRev{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        assert(_input.size() < (1ul << (64-BitsForPosition)));
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


        auto [bwt, csa] = [&input, &samplingRate, &inputSizes] () {
            auto sa  = createSA(input);
            auto bwt = createBWT(input, sa);
            auto csa = createCSA(std::move(sa), samplingRate, inputSizes);

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

    auto locate(size_t idx) const -> std::tuple<size_t, size_t> {
        auto opt = csa.value(idx);
        uint64_t steps{};
        while(!opt) {
            idx = occ.rank(idx, occ.symbol(idx));
            steps += 1;
            opt = csa.value(idx);
        }
        constexpr size_t posMask = (1ul<<BitsForPosition)-1;
        auto chr = opt.value() >> BitsForPosition;
        auto pos = (opt.value() & posMask) + steps;

        return {chr, pos};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, occRev, csa);
    }
};

}
