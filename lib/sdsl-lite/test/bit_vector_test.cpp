#include <string>

#include <sdsl/bit_vectors.hpp> // for rrr_vector

#include "common.hpp"

#include <gtest/gtest.h>

namespace
{

using namespace sdsl;
using namespace std;

string test_file;
string temp_file;
string temp_dir;

template <class T>
class bit_vector_test : public ::testing::Test
{};

template <class T>
class bit_vector_test_bv_only : public ::testing::Test
{};

using testing::Types;

#ifdef FULL_TEST_SUITE

typedef Types<bit_vector,
              bit_vector_il<>,
              bit_vector_il<64>,
              bit_vector_il<128>,
              bit_vector_il<256>,
              bit_vector_il<1024>,
              rrr_vector<64>,
              rrr_vector<256>,
              rrr_vector<129>,
              rrr_vector<192>,
              rrr_vector<255>,
              rrr_vector<15>,
              rrr_vector<31>,
              rrr_vector<63>,
              rrr_vector<83>,
              rrr_vector<127>,
              rrr_vector<128>,
              sd_vector<>,
              sd_vector<rrr_vector<63>>,
              hyb_vector<>>
                                                  Implementations;

#else

typedef Types<bit_vector, bit_vector_il<>, rrr_vector<>, sd_vector<>, hyb_vector<>> Implementations;

#endif

typedef Types<bit_vector> Implementations_BV_Only;

TYPED_TEST_SUITE(bit_vector_test, Implementations, );

//! Test operator==
TYPED_TEST(bit_vector_test, equality_operator)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam c_bv(bv);
    ASSERT_TRUE(store_to_file(c_bv, temp_file));
    TypeParam cc_bv;
    ASSERT_TRUE(load_from_file(cc_bv, temp_file));
    ASSERT_EQ(c_bv.size(), cc_bv.size());
    ASSERT_EQ(c_bv, cc_bv);
}

//! Test operator[]
TYPED_TEST(bit_vector_test, access)
{
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam c_bv(bv);
    ASSERT_EQ(bv.size(), c_bv.size());
    for (uint64_t j = 0; j < bv.size(); ++j) { ASSERT_EQ((bool)(bv[j]), (bool)(c_bv[j])); }
    TypeParam mo_bv = TypeParam(bv);
    ASSERT_EQ(bv.size(), mo_bv.size());
    for (uint64_t j = 0; j < bv.size(); ++j) { ASSERT_EQ((bool)(bv[j]), (bool)(c_bv[j])); }
}

TYPED_TEST(bit_vector_test, get_int)
{
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam c_bv(bv);
    ASSERT_EQ(bv.size(), c_bv.size());
    uint8_t len = 63;
    for (uint64_t j = 0; j + len < bv.size(); j += len) { ASSERT_EQ(bv.get_int(j, len), c_bv.get_int(j, len)); }
}

TYPED_TEST(bit_vector_test, get_int_all_block_sizes)
{
    bit_vector bv(10000, 0);
    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, 9);
    auto dice = bind(distribution, rng);
    for (size_t i = 1001; i < bv.size(); ++i)
    {
        if (0 == dice()) bv[i] = 1;
    }

    TypeParam c_bv(bv);
    for (uint8_t len = 1; len <= 64; ++len)
    {
        for (size_t i = 0; i + len <= bv.size(); ++i)
        {
            ASSERT_EQ(bv.get_int(i, len), c_bv.get_int(i, len)) << "i=" << i << " len=" << (int)len << endl;
        }
    }
}

TYPED_TEST(bit_vector_test, swap)
{
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam c_bv(bv);
    ASSERT_EQ(bv.size(), c_bv.size());
    TypeParam bv_empty;
    ASSERT_EQ((uint64_t)0, bv_empty.size());
    std::swap(bv_empty, c_bv);
    ASSERT_EQ((uint64_t)0, c_bv.size());
    ASSERT_EQ(bv.size(), bv_empty.size());
    for (uint64_t j = 0; j < bv.size(); ++j) { ASSERT_EQ((bool)(bv[j]), (bool)(bv_empty[j])); }
}

TYPED_TEST(bit_vector_test, delete_)
{
    sdsl::remove(temp_file);
}

TYPED_TEST_SUITE(bit_vector_test_bv_only, Implementations_BV_Only, );

// clang-format off
#define LFSR_START 0x00000001    // linear-feedback shift register with
// clang-format on
#define LFSR_FEEDBACK 0x0110F65C // .. period 33554431 = 31*601*1801
#define LFSR_NEXT(x) (((x) >> 1) ^ (((x)&1) * LFSR_FEEDBACK))
// nota bene: LFSR output has ~50% 1s, will bias compression types like RRR

TYPED_TEST(bit_vector_test_bv_only, and_with)
{
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam bv1(bv);

    TypeParam bv2(bv1.size(), 0);
    uint32_t lfsr = LFSR_START;
    for (size_t i = 0; i < bv1.size(); ++i)
    {
        lfsr = LFSR_NEXT(lfsr);
        bv2[i] = lfsr & 1;
    }

    bv2 &= bv1;

    lfsr = LFSR_START;
    for (size_t i = 0; i < bv1.size(); ++i)
    {
        lfsr = LFSR_NEXT(lfsr);
        ASSERT_EQ((bit_vector::value_type)bv2[i], bv1[i] & (lfsr & 1)) << "i=" << i << endl;
    }
}

TYPED_TEST(bit_vector_test_bv_only, or_with)
{
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam bv1(bv);

    TypeParam bv2(bv1.size(), 0);
    uint32_t lfsr = LFSR_START;
    for (size_t i = 0; i < bv1.size(); ++i)
    {
        lfsr = LFSR_NEXT(lfsr);
        bv2[i] = lfsr & 1;
    }

    bv2 |= bv1;

    lfsr = LFSR_START;
    for (size_t i = 0; i < bv1.size(); ++i)
    {
        lfsr = LFSR_NEXT(lfsr);
        ASSERT_EQ((bit_vector::value_type)bv2[i], bv1[i] | (lfsr & 1)) << "i=" << i << endl;
    }
}

TYPED_TEST(bit_vector_test_bv_only, xor_with)
{
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam bv1(bv);

    TypeParam bv2(bv1.size(), 0);
    uint32_t lfsr = LFSR_START;
    for (size_t i = 0; i < bv1.size(); ++i)
    {
        lfsr = LFSR_NEXT(lfsr);
        bv2[i] = lfsr & 1;
    }

    bv2 ^= bv1;

    lfsr = LFSR_START;
    for (size_t i = 0; i < bv1.size(); ++i)
    {
        lfsr = LFSR_NEXT(lfsr);
        ASSERT_EQ((bit_vector::value_type)bv2[i], bv1[i] ^ (lfsr & 1)) << "i=" << i << endl;
    }
}

#if SDSL_HAS_CEREAL
template <typename in_archive_t, typename out_archive_t, typename TypeParam>
void do_serialisation(TypeParam const & l)
{
    {
        std::ofstream os{ temp_file, std::ios::binary };
        out_archive_t oarchive{ os };
        oarchive(l);
    }

    TypeParam in_l{};
    {
        std::ifstream is{ temp_file, std::ios::binary };
        in_archive_t iarchive{ is };
        iarchive(in_l);
    }
    EXPECT_EQ(l, in_l);
}

TYPED_TEST(bit_vector_test, cereal)
{
    if (temp_dir != "@/")
    {
        bit_vector bv;
        ASSERT_TRUE(load_from_file(bv, test_file));
        TypeParam bv1(bv);
        ;

        do_serialisation<cereal::BinaryInputArchive, cereal::BinaryOutputArchive>(bv1);
        do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(bv1);
        do_serialisation<cereal::JSONInputArchive, cereal::JSONOutputArchive>(bv1);
        do_serialisation<cereal::XMLInputArchive, cereal::XMLOutputArchive>(bv1);
    }
}
#endif // SDSL_HAS_CEREAL

} // end namespace

int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    if (init_2_arg_test(argc, argv, "BV", test_file, temp_dir, temp_file) != 0) { return 1; }
    return RUN_ALL_TESTS();
}
