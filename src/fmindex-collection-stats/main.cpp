// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include <fmindex-collection/rankvector/rankvector.h>
#include <fmindex-collection/bitvector/all.h>
#include <fmt/format.h>
#include <filesystem>
#include <thread>
#include <fstream>

template <size_t Sigma>
static auto generateString(size_t l) {
    auto buffer = std::vector<uint8_t>{};
    buffer.reserve(l);
    for (size_t i{0}; i < l; ++i) {
        buffer.push_back(rand() % Sigma);
    }
    return buffer;
}

static size_t pow10(size_t e) {
    if (e == 0) return 1;
    if (e == 1) return 10;
    auto v = pow10(e/2);
    v = v*v;
    if (e%2 == 1) {
        v *= 10;
    }
    return v;
}

template <typename RV>
void f(std::string label, std::vector<uint8_t> const& text) {
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


int main() {
    using namespace fmindex_collection::rankvector;
    for (auto e : {9}) {
        constexpr static size_t Sigma = 5;
        auto text = generateString<Sigma>(pow10(e));
//        f<Naive<Sigma>>(text);
        f<MultiBitvector<Sigma>>("MultiBitvector", text);
        f<InterleavedBitvector8<Sigma>>("InterleavedBitvector8", text);
        f<InterleavedBitvector16<Sigma>>("InterleavedBitvector16", text);
        f<InterleavedBitvector32<Sigma>>("InterleavedBitvector32", text);
        f<InterleavedEPR8<Sigma>>("InterleavedEPR8", text);
        f<InterleavedEPR16<Sigma>>("InterleavedEPR16", text);
        f<InterleavedEPR32<Sigma>>("InterleavedEPR32", text);
        f<InterleavedEPRV2_8<Sigma>>("InterleavedEPRV2_8", text);
        f<InterleavedEPRV2_16<Sigma>>("InterleavedEPRV2_16", text);
        f<InterleavedEPRV2_32<Sigma>>("InterleavedEPRV2_32", text);
        f<EPRV3_8<Sigma>> ("EPRV3_8", text);
        f<EPRV3_16<Sigma>>("EPRV3_16", text);
        f<EPRV3_32<Sigma>>("EPRV3_32", text);
        f<EPRV4<Sigma>> ("EPRV4", text);
        f<EPRV5<Sigma>> ("EPRV5", text);
        f<DenseEPRV6<Sigma>> ("DenseEPRV6", text);
        f<InterleavedEPRV7<Sigma>>("InterleavedEPRV7", text);
        f<InterleavedWavelet<Sigma>>("InterleavedWavelet", text);
        f<Sdsl_wt_bldc<Sigma>>("sdsl-wavelet", text);
        f<Wavelet<Sigma>>("Wavelet", text);
    }
}
