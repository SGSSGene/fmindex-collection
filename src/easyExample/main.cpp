#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/search/all.h>
#include <search_schemes/generator/all.h>
#include <search_schemes/expand.h>

#include <cereal/archives/binary.hpp>
#include <fmt/format.h>
#include <unordered_set>

using namespace fmindex_collection;

constexpr size_t Sigma = 5;

template <size_t Sigma>
using Table = occtable::compact2::OccTable<Sigma>;


int main(int argc, char const* const* argv) {
    auto reference = std::vector<uint8_t>{};
    auto index = BiFMIndex<Table<Sigma>>{std::move(reference), 16};

    auto queries = std::vector<std::vector<uint8_t>>{};

    auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, 1), queries[0].size());

    search_pseudo::search<true>(index, queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {
        fmt::print("found something {}\n", cursor.count());
    });
    return 0;
}


