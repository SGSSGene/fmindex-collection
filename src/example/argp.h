#pragma once

struct Config {
    std::string generator = "h2";
    size_t maxQueries{};
    size_t readLength{};
    bool saveOutput{false};
    size_t minK{0}, maxK{6}, k_stepSize{1};
    bool reverse{true};
    bool help{false};

    std::vector<std::string> algorithms;

    std::string queryPath{};
    std::string indexPath{};
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
        } else if (argv[i] == std::string{"--gen"} and i+1 < argc) {
            ++i;
            config.generator = argv[i];
        } else if (argv[i] == std::string{"--queries"} and i+1 < argc) {
            ++i;
            config.maxQueries = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--read_length"} and i+1 < argc) {
            ++i;
            config.readLength = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--save_output"}) {
            config.saveOutput = true;
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
        } else {
            throw std::runtime_error("unknown commandline " + std::string{argv[i]});
        }
    }
    return config;
}
