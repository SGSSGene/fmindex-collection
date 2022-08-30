#pragma once

#include <set>
#include <stdexcept>
#include <string>
#include <vector>

struct Config {
    std::string generator = "h2-k2";
    bool generator_dyn = false;
    size_t maxQueries{};
    size_t readLength{};
    std::filesystem::path saveOutput;
    size_t minK{0}, maxK{6}, k_stepSize{1};
    bool reverse{true};
    bool help{false};
    std::set<std::string> extensions;

    std::vector<std::string> algorithms;

    std::string queryPath{};
    std::string indexPath{};

    enum class Mode {
        All,
        BestHits,
    };
    Mode mode {Mode::All};
    size_t maxHitsPerQuery{0};
};

auto loadConfig(int argc, char const* const* argv) {
    Config config;
    for (int i{1}; i < argc; ++i) {
        if (argv[i] == std::string{"--query"} and i+1 < argc) {
            ++i;
            config.queryPath = argv[i];
        } else if (argv[i] == std::string{"--index"} and i+1 < argc) {
            ++i;
            config.indexPath = argv[i];
        } else if (argv[i] == std::string{"--algo"} and i+1 < argc) {
            ++i;
            config.algorithms.emplace_back(argv[i]);
        } else if (argv[i] == std::string{"--ext"} and i+1 < argc) {
            ++i;
            config.extensions.emplace(argv[i]);
        } else if (argv[i] == std::string{"--gen"} and i+1 < argc) {
            ++i;
            config.generator = argv[i];
            if (config.generator.size() > 4 and config.generator.substr(config.generator.size()-4) == "_dyn") {
                config.generator = config.generator.substr(0, config.generator.size()-4);
                config.generator_dyn = true;
            }
        } else if (argv[i] == std::string{"--queries"} and i+1 < argc) {
            ++i;
            config.maxQueries = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--read_length"} and i+1 < argc) {
            ++i;
            config.readLength = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--save_output"} and i+1 < argc) {
            ++i;
            config.saveOutput = argv[i];
        } else if (argv[i] == std::string{"--min_k"} and i+1 < argc) {
            ++i;
            config.minK = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--max_k"} and i+1 < argc) {
            ++i;
            config.maxK = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--stepSize_k"} and i+1 < argc) {
            ++i;
            config.k_stepSize = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--no-reverse"}) {
            config.reverse = false;
        } else if (argv[i] == std::string{"--help"}) {
            config.help = true;
        } else if (argv[i] == std::string{"--mode"} and i+1 < argc) {
            ++i;
            auto s = std::string{argv[i]};
            if (s == "all") {
                config.mode = Config::Mode::All;
            } else if (s == "besthits") {
                config.mode = Config::Mode::BestHits;
            } else {
                throw std::runtime_error("invalid mode \"" + s + "\", must be any of \"all\", \"besthits\"");
            }
        } else if (argv[i] == std::string{"--maxhitperquery"} and i+1 < argc) {
            ++i;
            config.maxHitsPerQuery = std::stod(argv[i]);
        } else {
            throw std::runtime_error("unknown commandline " + std::string{argv[i]});
        }
    }
    return config;
}
