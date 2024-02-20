#include <fmindex-collection/rankvector/rankvector.h>
#include <fmt/format.h>
#include <sys/resource.h>
#include <thread>
#include <fstream>

// Estimate current ram usage
static size_t ramUsage() {
    {
        auto fs = std::ofstream{"/proc/self/clear_refs"};
        fs << 5;
    }

    rusage r{};
    getrusage(RUSAGE_SELF, &r);
    return r.ru_maxrss * 1024;

/*    auto fs = std::ifstream{"/proc/self/smaps_rollup"};
    auto line = std::string{};
    for (size_t i{0}; i < 1; ++i) {
        std::getline(fs, line);
    }

    size_t value;
    fs >> line;
    fs >> value;
    fmt::print("{}: {}\n", line, value);
    return value * 1024;*/
}

struct ScopedRam {
    size_t start = ramUsage();
    ScopedRam() {
//        auto fs = std::ofstream{"/proc/self/clear_refs"};
//        fs << 1;
    }
    auto operator*() const -> int64_t {
        return (int64_t)ramUsage() - start;
    }
    void print() const {
        fmt::print("ram usage: {}\n", to_string());
    }
    auto to_string() const -> std::string {
        auto v = **this;
        auto suffix = std::string{};
        if (v > int64_t{1024}*1024*1024*100) {
            v = v/1024/1024/1024;
            suffix = "gb";
        } else if (v > 1024*1024*100) {
            v = v/1024/1024;
            suffix = "mb";
        } else if (v > 1024*100) {
            v = v/1024;
            suffix = "kb";
        }
        return fmt::format("{}{}", v, suffix);
    }
};

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
