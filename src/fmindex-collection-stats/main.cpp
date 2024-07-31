// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <fmindex-collection/rankvector/rankvector.h>
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
    using namespace fmindex_collection::bitvector;


    double onePercentage;
    auto f = [&](auto const& text) {
        fmt::print("length: {}, percantage of ones: {}%\n", text.size(), onePercentage);
        analyse_bitvector<Bitvector>("Bitvector", text);
        analyse_bitvector<CompactBitvector>("CompactBitvector", text);
        analyse_bitvector<CompactBitvector4Blocks>("CompactBitvector4Blocks", text);
        analyse_bitvector<SparseBLEBitvector<1>>("SparseBLEBitvector  2", text);
        analyse_bitvector<SparseBLEBitvector<2>>("SparseBLEBitvector  4", text);
        analyse_bitvector<SparseBLEBitvector<3>>("SparseBLEBitvector  8", text);
        analyse_bitvector<SparseBLEBitvector<4>>("SparseBLEBitvector 16", text);
        analyse_bitvector<SparseBLEBitvector<5>>("SparseBLEBitvector 32", text);
        analyse_bitvector<SparseBLEBitvector<6>>("SparseBLEBitvector 64", text);
        analyse_bitvector<SparseBLEBitvector<1, SparseBLEBitvector<-1>>>("SparseBLEBitvector  2/2", text);
        analyse_bitvector<SparseBLEBitvector<2, SparseBLEBitvector<-1>>>("SparseBLEBitvector  4/2", text);
        analyse_bitvector<SparseBLEBitvector<2, SparseBLEBitvector<-2>>>("SparseBLEBitvector  4/4", text);
        analyse_bitvector<SparseBLEBitvector<3, SparseBLEBitvector<-1>>>("SparseBLEBitvector  8/2", text);
        analyse_bitvector<SparseBLEBitvector<3, SparseBLEBitvector<-2>>>("SparseBLEBitvector  8/4", text);
        analyse_bitvector<SparseBLEBitvector<3, SparseBLEBitvector<-3>>>("SparseBLEBitvector  8/8", text);

        analyse_bitvector<SparseBLEBitvector<2, Bitvector, SparseBLEBitvector<1>>>("SparseBLEBitvector  4/-2", text);
        analyse_bitvector<SparseBLEBitvector<3, Bitvector, SparseBLEBitvector<1>>>("SparseBLEBitvector  8/-2", text);
        analyse_bitvector<SparseBLEBitvector<3, Bitvector, SparseBLEBitvector<2>>>("SparseBLEBitvector  8/-4", text);

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
void analyse_rankvector(std::string label, std::vector<uint8_t> const& text) {
    auto rankvector = RV{text};
    // Save index to disk
    {
        auto fs     = std::ofstream("tmp.index", std::ios::binary);
        auto archive = cereal::BinaryOutputArchive{fs};
        archive(rankvector);
    }

    //auto size = *ram;
    auto size = std::filesystem::file_size("tmp.index");

    auto bits_per_char = size * 8. / text.size();
    fmt::print("{:<24} {: 8.2f} bits per character\n", label, bits_per_char);

    (void)rankvector;
}

template <int64_t BL, typename Bitvector = fmindex_collection::bitvector::Bitvector, typename Bitvector2 = fmindex_collection::bitvector::Bitvector>
using SparseBitvector = fmindex_collection::bitvector::SparseBLEBitvector<BL, Bitvector, Bitvector2>;

static void analyse_rankvectors() {
    using namespace fmindex_collection::rankvector;
    for (auto e : {7}) {
        constexpr static size_t Sigma = 21;
        auto text = generateString<Sigma>(my_pow10(e));
//        f<Naive<Sigma>>(text);
        analyse_rankvector<MultiBitvector<Sigma>>("MultiBitvector", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<1>>>("SparseMultiBitvector 2", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<2>>>("SparseMultiBitvector 4", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<3>>>("SparseMultiBitvector 8", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<4>>>("SparseMultiBitvector 16", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<5>>>("SparseMultiBitvector 32", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<6>>>("SparseMultiBitvector 64", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<1, SparseBitvector<-1>>>>("SparseMultiBitvector 2/2", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<2, SparseBitvector<-1>>>>("SparseMultiBitvector 4/2", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<2, SparseBitvector<-2>>>>("SparseMultiBitvector 4/4", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<2, fmindex_collection::bitvector::Bitvector, SparseBitvector<1>>>>("SparseMultiBitvector 4/-2", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<3, fmindex_collection::bitvector::Bitvector, SparseBitvector<1>>>>("SparseMultiBitvector 8/-2", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<3, fmindex_collection::bitvector::Bitvector, SparseBitvector<2>>>>("SparseMultiBitvector 8/-4", text);
        analyse_rankvector<MultiBitvector<Sigma, SparseBitvector<4, fmindex_collection::bitvector::Bitvector, SparseBitvector<2>>>>("SparseMultiBitvector 16/-4", text);

        analyse_rankvector<InterleavedBitvector8<Sigma>>("InterleavedBitvector8", text);
        analyse_rankvector<InterleavedBitvector16<Sigma>>("InterleavedBitvector16", text);
        //analyse_rankvector<InterleavedBitvector32<Sigma>>("InterleavedBitvector32", text);
        analyse_rankvector<InterleavedEPR8<Sigma>>("InterleavedEPR8", text);
        analyse_rankvector<InterleavedEPR16<Sigma>>("InterleavedEPR16", text);
        //analyse_rankvector<InterleavedEPR32<Sigma>>("InterleavedEPR32", text);
        analyse_rankvector<InterleavedEPRV2_8<Sigma>>("InterleavedEPRV2_8", text);
        analyse_rankvector<InterleavedEPRV2_16<Sigma>>("InterleavedEPRV2_16", text);
        //analyse_rankvector<InterleavedEPRV2_32<Sigma>>("InterleavedEPRV2_32", text);
        analyse_rankvector<EPRV3_8<Sigma>> ("EPRV3_8", text);
        analyse_rankvector<EPRV3_16<Sigma>>("EPRV3_16", text);
        //analyse_rankvector<EPRV3_32<Sigma>>("EPRV3_32", text);
        analyse_rankvector<EPRV4<Sigma>> ("EPRV4", text);
        analyse_rankvector<EPRV5<Sigma>> ("EPRV5", text);
        analyse_rankvector<DenseEPRV6<Sigma>> ("DenseEPRV6", text);
        analyse_rankvector<InterleavedEPRV7<Sigma>>("InterleavedEPRV7", text);
        analyse_rankvector<InterleavedWavelet<Sigma>>("InterleavedWavelet", text);
#if FMC_USE_SDSL
        analyse_rankvector<Sdsl_wt_bldc<Sigma>>("sdsl-wavelet", text);
#endif
        analyse_rankvector<Wavelet<Sigma>>("Wavelet", text);
    }
}

int main() {
    fmt::print("analyse bitvectors:\n");
    analyse_bitvectors();

    fmt::print("\nanalyse rankvectors:\n");
    analyse_rankvectors();
}
