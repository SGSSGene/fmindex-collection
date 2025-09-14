// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <fmindex-collection/string/all.h>
#include <fmindex-collection/bitvector/all.h>
#include <fmt/format.h>
#include <filesystem>
#include <thread>
#include <fstream>

template <size_t Sigma, size_t Sparse=1>
static auto generateString(size_t l) {
    auto buffer = std::vector<uint8_t>{};
    buffer.reserve(l);
    for (size_t i{0}; i < l; ++i) {
        if (rand() % Sparse == 0) {
            buffer.push_back(rand() % Sigma);
        } else {
            buffer.push_back(0);
        }
    }
    return buffer;
}


static size_t my_pow10(size_t e) {
    if (e == 0) return 1;
    if (e == 1) return 10;
    auto v = my_pow10(e/2);
    v = v*v;
    if (e%2 == 1) {
        v *= 10;
    }
    return v;
}

template <typename Vector>
void analyse_bitvector(std::string label, std::vector<uint8_t> const& text) {
    auto vector = Vector{text};
    // Save index to disk
    {
        auto fs     = std::ofstream("tmp.index", std::ios::binary);
        auto archive = cereal::BinaryOutputArchive{fs};
        archive(vector);
    }

    //auto size = *ram;
    auto size = std::filesystem::file_size("tmp.index");

    auto bits_per_char = size * 8. / text.size();
    fmt::print("{:<24} {: 8.2f} bits per bit\n", label, bits_per_char);

    (void)vector;
}

static void analyse_bitvectors() {
    using namespace fmc::bitvector;


    double onePercentage;
    auto f = [&](auto const& text) {
        fmt::print("length: {}, percantage of ones: {}%\n", text.size(), onePercentage);
        analyse_bitvector<Bitvector>("Bitvector", text);
        analyse_bitvector<CompactBitvector>("CompactBitvector", text);
        analyse_bitvector<CompactBitvector4Blocks>("CompactBitvector4Blocks", text);
        analyse_bitvector<SparseRBBitvector<1>>("SparseRBBitvector  2", text);
        analyse_bitvector<SparseRBBitvector<2>>("SparseRBBitvector  4", text);
        analyse_bitvector<SparseRBBitvector<3>>("SparseRBBitvector  8", text);
        analyse_bitvector<SparseRBBitvector<4>>("SparseRBBitvector 16", text);
        analyse_bitvector<SparseRBBitvector<5>>("SparseRBBitvector 32", text);
        analyse_bitvector<SparseRBBitvector<6>>("SparseRBBitvector 64", text);
        analyse_bitvector<SparseRBBitvector<2, Bitvector, SparseRBBitvector<1>>>("SparseRBBitvector  4/2", text);
        analyse_bitvector<SparseRBBitvector<3, Bitvector, SparseRBBitvector<1>>>("SparseRBBitvector  8/2", text);
        analyse_bitvector<SparseRBBitvector<3, Bitvector, SparseRBBitvector<2>>>("SparseRBBitvector  8/4", text);

    };

    for (auto e : {7}) {
        {
            onePercentage=50;
            auto text = generateString<2>(my_pow10(e));
            f(text);
        }
        {
            onePercentage=25;
            auto sparseText25 = generateString<2, 2>(my_pow10(e));
            f(sparseText25);
        }
        {
            onePercentage=10;
            auto sparseText10 = generateString<2, 5>(my_pow10(e));
            f(sparseText10);
        }
        {
            onePercentage=5;
            auto sparseText5 = generateString<2, 10>(my_pow10(e));
            f(sparseText5);
        }
        {
            onePercentage=0.5;
            auto sparseText0_5 = generateString<2, 100>(my_pow10(e));
            f(sparseText0_5);
        }
    }
}

template <typename RV>
void analyse_string(std::string label, std::vector<uint8_t> const& text) {
    auto string = RV{text};
    // Save index to disk
    {
        auto fs     = std::ofstream("tmp.index", std::ios::binary);
        auto archive = cereal::BinaryOutputArchive{fs};
        archive(string);
    }

    //auto size = *ram;
    auto size = std::filesystem::file_size("tmp.index");

    auto bits_per_char = size * 8. / text.size();
    fmt::print("{:<24} {: 8.2f} bits per character\n", label, bits_per_char);

    (void)string;
}

template <int64_t BL, typename Bitvector = fmc::bitvector::Bitvector, typename Bitvector2 = fmc::bitvector::Bitvector>
using SparseBitvector = fmc::bitvector::SparseRBBitvector<BL, Bitvector, Bitvector2>;

static void analyse_strings() {
    using namespace fmc::string;
    for (auto e : {7}) {
        constexpr static size_t Sigma = 21;
        auto text = generateString<Sigma>(my_pow10(e));
//        f<Naive<Sigma>>(text);
        analyse_string<MultiBitvector<Sigma>>("MultiBitvector", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<1>>>("SparseMultiBitvector 2", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<2>>>("SparseMultiBitvector 4", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<3>>>("SparseMultiBitvector 8", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<4>>>("SparseMultiBitvector 16", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<5>>>("SparseMultiBitvector 32", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<6>>>("SparseMultiBitvector 64", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<2, fmc::bitvector::Bitvector, SparseBitvector<1>>>>("SparseMultiBitvector 4/2", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<3, fmc::bitvector::Bitvector, SparseBitvector<1>>>>("SparseMultiBitvector 8/2", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<3, fmc::bitvector::Bitvector, SparseBitvector<2>>>>("SparseMultiBitvector 8/4", text);
        analyse_string<MultiBitvector<Sigma, SparseBitvector<4, fmc::bitvector::Bitvector, SparseBitvector<2>>>>("SparseMultiBitvector 16/4", text);

        analyse_string<InterleavedBitvector8<Sigma>>("InterleavedBitvector8", text);
        analyse_string<InterleavedBitvector16<Sigma>>("InterleavedBitvector16", text);
        //analyse_string<InterleavedBitvector32<Sigma>>("InterleavedBitvector32", text);
        analyse_string<InterleavedEPR8<Sigma>>("InterleavedEPR8", text);
        analyse_string<InterleavedEPR16<Sigma>>("InterleavedEPR16", text);
        //analyse_string<InterleavedEPR32<Sigma>>("InterleavedEPR32", text);
        analyse_string<InterleavedEPRV2_8<Sigma>>("InterleavedEPRV2_8", text);
        analyse_string<InterleavedEPRV2_16<Sigma>>("InterleavedEPRV2_16", text);
        //analyse_string<InterleavedEPRV2_32<Sigma>>("InterleavedEPRV2_32", text);
        analyse_string<EPRV3_8<Sigma>> ("EPRV3_8", text);
        analyse_string<EPRV3_16<Sigma>>("EPRV3_16", text);
        //analyse_string<EPRV3_32<Sigma>>("EPRV3_32", text);
        analyse_string<EPRV4<Sigma>> ("EPRV4", text);
        analyse_string<EPRV5<Sigma>> ("EPRV5", text);
        analyse_string<DenseEPRV6<Sigma>> ("DenseEPRV6", text);
        analyse_string<InterleavedEPRV7<Sigma>>("InterleavedEPRV7", text);
        analyse_string<InterleavedWavelet<Sigma>>("InterleavedWavelet", text);
#if FMC_USE_SDSL
        analyse_string<Sdsl_wt_bldc<Sigma>>("sdsl-wavelet", text);
#endif
        analyse_string<Wavelet<Sigma>>("Wavelet", text);
    }
}

int main() {
    fmt::print("analyse bitvectors:\n");
    analyse_bitvectors();

    fmt::print("\nanalyse strings:\n");
    analyse_strings();
}
