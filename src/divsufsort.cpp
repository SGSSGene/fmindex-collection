#include "FMIndex_SimpleOcc.h"
#include "FMIndex_CompactOcc.h"
#include "FMIndex_CompactOcc_Aligned.h"
#include "FMIndex_CompactOcc_Prefix.h"
#include "FMIndex_CompactOcc2.h"
#include "FMIndex_CompactOcc2_Aligned.h"
#include "FMIndex_Bitvector.h"
#include "FMIndex_Bitvector_Prefix.h"

#include "random.h"
#include "StopWatch.h"

#include <divsufsort64.h>

#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


inline auto construct_bwt_from_sa(std::vector<int64_t> const& sa, std::string_view const& text) -> std::vector<uint8_t> {
    assert(sa.size() == text.size());
    std::vector<uint8_t> bwt;
    bwt.resize(text.size());
    for (size_t i{0}; i < sa.size(); ++i) {
        bwt[i] = text[(sa[i] + text.size()- 1) % text.size()];
    }
    return bwt;
}
inline auto readFile(std::filesystem::path const& file) -> std::vector<uint8_t> {
    auto ifs = std::ifstream{file, std::ios::binary};
    ifs.seekg(0, std::ios::end);
    std::size_t size = ifs.tellg();
    auto buffer = std::vector<uint8_t>(size);
    ifs.seekg(0, std::ios::beg);
    ifs.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
    return buffer;
}


template <FMIndex Index, typename T>
void printFMIndex(Index const& index, T const& bwt) {
    if (bwt.size() > 16) return;
    for (size_t i{0}; i < bwt.size(); ++i) {
        std::cout << i << " " << int(bwt[i]) << "\t";
        for (size_t j{0}; j < index.Sigma; ++j) {
            std::cout << " " << index.rank(j, i);
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";

    for (size_t i{0}; i < bwt.size(); ++i) {
        std::cout << i << " " << int(bwt[i]) << "\t";
        for (size_t j{0}; j < index.Sigma; ++j) {
            std::cout << " " << index.prefix_rank(j, i);
        }
        std::cout << "\n";
    }
}
template <FMIndex Index, typename T>
void benchmarkFMIndex(Index const& index, T const& bwt) {

    StopWatch watch;
    xorshf96_reset();
    uint64_t a{};
    for (size_t i{0}; i < 1'000'000; ++i) {
        auto symb = xorshf96() % Index::Sigma;
        auto row = xorshf96() % index.size();
        a += index.rank(symb, row);
    }
    auto time_accesstime = watch.reset();
    std::cout << "FMIndex - rank time: "<< time_accesstime << "s\n";
    std::cout << "a: " << a << "\n";

    a = 0;
    for (size_t i{0}; i < 1'000'000; ++i) {
        auto symb = xorshf96() % Index::Sigma;
        auto row = xorshf96() % index.size();
        a += index.prefix_rank(symb, row);
    }
    auto time_accesstime2 = watch.reset();
    std::cout << "FMIndex - rank2 time: "<< time_accesstime2 << "s\n";
    std::cout << "a: " << a << "\n";

    {
        uint64_t jumps{1};
        uint64_t pos = index.rank(bwt[0], 0);
        while (pos != 0 && jumps/2 < bwt.size()) {
            jumps += 1;
            pos = index.rank(bwt[pos], pos);
        }

        auto time_fullbackwardsjump = watch.reset();
        std::cout << "FMIndex - backwards: "<< time_fullbackwardsjump << "s\n";
        std::cout << "jumps: " << jumps << "/" << bwt.size() << "\n";
    }
    {
        uint64_t a{};
        uint64_t jumps{1};
        uint64_t pos = index.rank(bwt[0], 0);
        while (pos != 0 && jumps/2 < bwt.size()) {
            jumps += 1;
            pos = index.rank(bwt[pos], pos);
            a += index.prefix_rank(bwt[pos], pos);
        }

        auto time_fullbackwardsjump = watch.reset();
        std::cout << "FMIndex - backwards: "<< time_fullbackwardsjump << "s\n";
        std::cout << "jumps: " << jumps << "/" << bwt.size() << " " << a << "\n";
    }

}

template <FMIndex Index, typename T>
void constructIndex(std::string name, T const& bwt) {
    std::cout << "measuring " << name << "\n";
    size_t s = Index::expectedMemoryUsage(bwt.size());
    std::cout << "Expected memory consumption: " << s/(1024*1024) << "MB\n";

    if (s < 1024*1024*1024*8ul) {
        StopWatch watch;
        auto index = Index{bwt};

        auto time_fmconstruction = watch.reset();
        std::cout << "FMIndex - construction time: "<< time_fmconstruction << "s\n";
        printFMIndex(index, bwt);
        benchmarkFMIndex(index, bwt);
    }
    std::cout << "===\n\n";
}

int main() {
    StopWatch watch;

    constexpr size_t Sigma = 6;
/*    std::string text;
    text.resize(1ul<<30, '\0');
    for (size_t i{0}; i < text.size(); ++i) {
        text[i] = (xorshf96() % (Sigma-1))+1;
    }
    text.back() = '\0';

    auto time_generation = watch.reset();
    std::cout << "text size: " << text.size()/1024/1024 << " million chars, "<<  text.size()/1024/1024*std::log(Sigma)/std::log(2) << " million bits\n";
    std::cout << "text generation: "<< time_generation << "s\n";*/

/*    auto bwt = [&]() {
        std::vector<int64_t> sa;
        sa.resize(text.size());
        auto error = divsufsort64((uint8_t const*)text.data(), sa.data(), text.size());
        if (error != 0) {
            throw std::runtime_error("some error while creating the suffix array");
        }

        auto time_saconstruction = watch.reset();
        std::cout << "sa - construction time: "<< time_saconstruction << "s\n";

        return construct_bwt_from_sa(sa, text);
    }();*/
    auto bwt = readFile("/home/gene/hg38/text.dna5.bwt");

    auto time_bwtconstruction = watch.reset();
    std::cout << "bwt - construction time: "<< time_bwtconstruction << "s, length: " << bwt.size() << "\n";

    constructIndex<simpleocc::FMIndex<Sigma>>("simple", bwt);
    constructIndex<compactocc::FMIndex<Sigma>>("compact", bwt);
    constructIndex<compactocc_align::FMIndex<Sigma>>("compact_aligned", bwt);
    constructIndex<compactocc_prefix::FMIndex<Sigma>>("compact_prefix", bwt);
    constructIndex<compactocc2::FMIndex<Sigma>>("compact2", bwt);
    constructIndex<compactocc2_align::FMIndex<Sigma>>("compact2_aligned", bwt);
    constructIndex<bitvectorocc::FMIndex<Sigma>>("bitvector", bwt);
    constructIndex<bitvectorocc_prefix::FMIndex<Sigma>>("bitvector_prefix", bwt);

    return 0;
}


