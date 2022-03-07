#pragma once

#include "CSA.h"
#include "occtable/concepts.h"

#include <cassert>
#include <cmath>
#include <numeric>
#include <sdsl/divsufsort.hpp>

namespace fmindex_collection {

template <OccTable Table>
struct BiFMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    using TTable = Table;

    Table  occ;
    Table  occRev;
    CSA    csa;
    size_t bitsForPosition;
    size_t bitPositionMask;


    BiFMIndex(std::vector<uint8_t> const& bwt, std::vector<uint8_t> const& bwtRev, CSA _csa, size_t _bitsForPosition)
        : occ{bwt}
        , occRev{bwtRev}
        , csa{std::move(_csa)}
        , bitsForPosition{_bitsForPosition}
    {
        assert(_bitsForPosition < 64);
        bitPositionMask = (1ul<<bitsForPosition)-1;

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

    auto createCSA(std::vector<int64_t> sa, size_t samplingRate, std::vector<size_t> const& _inputSizes, std::function<size_t(size_t)> _label) -> CSA {
        auto bitStack = fmindex_collection::BitStack{};
        auto ssa      = std::vector<uint64_t>{};
        if (samplingRate > 0) {
            ssa.reserve(sa.size() / samplingRate + 1 + _inputSizes.size());
        }

        auto newLabels = std::vector<uint64_t>{};
        newLabels.resize(sa.size(), std::numeric_limits<uint64_t>::max());

        auto accIter = _inputSizes.begin();
        size_t subjId{};
        size_t subjPos{};

        size_t lastSamplingPos{};
        for (size_t i{0}; i < newLabels.size(); ++i, ++subjPos) {
            while (subjPos >= *accIter) {
                subjPos -= *accIter;
                ++subjId;
                ++accIter;
            }
            bool sample = (samplingRate == 0) || ((i-lastSamplingPos) % samplingRate == 0) || (subjPos == 0);
            if (sample) {
                lastSamplingPos = i;
                newLabels[i] = (subjPos) | (_label(subjId) << bitsForPosition);
            }
        }

        for (size_t i{0}; i < sa.size(); ++i) {
            bool sample = newLabels[sa[i]] != std::numeric_limits<uint64_t>::max();
            bitStack.push(sample);
            if (sample) {
                ssa.push_back(newLabels[sa[i]]);
            }
        }
        return CSA{std::move(ssa), bitStack, samplingRate};
    }

    BiFMIndex(std::vector<uint8_t> _input, size_t samplingRate)
        : occ{cereal_tag{}}
        , occRev{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        auto input = std::vector<std::vector<uint8_t>>{std::move(_input)};
        *this = BiFMIndex{std::move(input), samplingRate};
    }

    BiFMIndex(std::vector<std::vector<uint8_t>> _input, size_t samplingRate, std::vector<size_t> _labels = {})
        : occ{cereal_tag{}}
        , occRev{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        size_t bitsForSeqId = std::max(1ul, size_t(std::ceil(std::log2(_input.size()))));
        assert(bitsForSeqId < 64);

        bitsForPosition = 64 - bitsForSeqId;
        bitPositionMask = (1ul<<bitsForPosition)-1;

        size_t totalSize = std::accumulate(begin(_input), end(_input), size_t{0}, [](auto s, auto const& l) { return s + l.size() + 1; });

        auto input = std::vector<uint8_t>{};
        input.reserve(totalSize);

        auto inputSizes = std::vector<size_t>{};
        inputSizes.reserve(_input.size());

        for (auto const& l : _input) {
            if (l.size() >= std::pow(2, bitsForPosition)) {
                throw std::runtime_error("sequence are to long and to many. Of 64bit available " + std::to_string(bitsForPosition) + "bits are required for positions and " + std::to_string(bitsForSeqId) + "bits are required for the sequence id");
            }
            input.insert(end(input), begin(l), end(l));
            input.emplace_back(0);
            inputSizes.emplace_back(l.size()+1);
        }
        decltype(_input){}.swap(_input); // input memory can be deleted


        auto [bwt, csa] = [&input, &samplingRate, &inputSizes, &_labels, this] () {
            auto sa  = createSA(input);
            auto bwt = createBWT(input, sa);
            auto csa = [&]() {
                if (_labels.empty()) {
                    return createCSA(std::move(sa), samplingRate, inputSizes, [](size_t i) { return i; });
                } else {
                    return createCSA(std::move(sa), samplingRate, inputSizes, [&](size_t i) {
                        return _labels[i];
                    });
                }
            }();

            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        auto bwtRev = [&]() {
            std::reverse(begin(input), end(input));
            auto saRev  = createSA(input);
            auto bwtRev = createBWT(input, saRev);
            return bwtRev;
        }();

        decltype(input){}.swap(input); // input memory can be deleted

        *this = BiFMIndex{bwt, bwtRev, std::move(csa), bitsForPosition};
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
            idx = [&]() {
                if constexpr (requires(Table t) { { t.rank_symbol(size_t{}) }; }) {
                    return occ.rank_symbol(idx);
                } else {
                    return occ.rank(idx, occ.symbol(idx));
                }
            }();

//            idx = occ.rank(idx, occ.symbol(idx));
            steps += 1;
            opt = csa.value(idx);
        }
        auto chr = opt.value() >> bitsForPosition;
        auto pos = (opt.value() & bitPositionMask) + steps;

        return {chr, pos};
    }

    auto locate(size_t idx, size_t maxSteps) const -> std::optional<std::tuple<size_t, size_t>> {
        auto opt = csa.value(idx);
        uint64_t steps{};
        for (;!opt and maxSteps > 0; --maxSteps) {
            idx = occ.rank(idx, occ.symbol(idx));
            steps += 1;
            opt = csa.value(idx);
        }
        if (!opt) {
            return std::nullopt;
        }
        auto chr = opt.value() >> bitsForPosition;
        auto pos = (opt.value() & bitPositionMask) + steps;

        return std::make_tuple(chr, pos);
    }


    auto single_locate_step(size_t idx) const -> std::optional<std::tuple<size_t, size_t>> {
        auto opt = csa.value(idx);
        if (!opt) return std::nullopt;
        auto chr = opt.value() >> bitsForPosition;
        auto pos = (opt.value() & bitPositionMask);

        return std::make_tuple(chr, pos);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, occRev, csa);
        ar(bitsForPosition, bitPositionMask);
    }
};

}
