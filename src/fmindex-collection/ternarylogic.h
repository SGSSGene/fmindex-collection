#pragma once

#include <array>
#include <bitset>
#include <cstdint>
#include <functional>

template <size_t R, size_t N, typename T=std::bitset<N>>
auto ternarylogic_impl(T const& a, T const& b, T const& c) -> T {
    if constexpr (R == 0x00) {
        return T{};
    } else if constexpr (R == 0x01) {
        return ~(a | b | c);
    } else if constexpr (R == 0x02) {
        return c & ~(a | b);
    } else if constexpr (R == 0x03) {
        return ~(a | b);
    } else if constexpr (R == 0x04) {
        return b & ~(a | c);
    } else if constexpr (R == 0x05) {
        return ~(a | c);
    } else if constexpr (R == 0x06) {
        return ~a & (b ^~ c);
    } else if constexpr (R == 0x07) {
        return ~a & (~b | ~c);
    } else if constexpr (R == 0x08) {
        return b & c & ~a;
    } else if constexpr (R == 0x09) {
        return ~a & (b ^ c);
    } else if constexpr (R == 0x0a) {
        return c & ~a;
    } else if constexpr (R == 0x0b) {
        return ~a & (c | ~b);
    } else if constexpr (R == 0x0c) {
        return b & ~a;
    } else if constexpr (R == 0x0d) {
        return ~a & (b | ~c);
    } else if constexpr (R == 0x0e) {
        return ~a & (b | c);
    } else if constexpr (R == 0x0f) {
        return ~a;
    } else if constexpr (R == 0x10) {
        return a & ~b & ~c;
    } else if constexpr (R == 0x11) {
        return ~b & ~c;
    } else if constexpr (R == 0x12) {
        return ~b & (a | c) & (~a | ~c);
    } else if constexpr (R == 0x13) {
        return ~b & (~a | ~c);
    } else if constexpr (R == 0x14) {
        return ~c & (a | b) & (~a | ~b);
    } else if constexpr (R == 0x15) {
        return ~c & (~a | ~b);
    } else if constexpr (R == 0x16) {
        return a ^ b ^ c ^ (a & b & c);
    } else if constexpr (R == 0x17) {
        return (~a | ~b) & (~a | ~c) & (~b | ~c);
    } else if constexpr (R == 0x18) {
        return (b & c & ~a) | (a & ~b & ~c);
    } else if constexpr (R == 0x19) {
        return (~b & ~c) | (b & c & ~a);
    } else if constexpr (R == 0x1a) {
        return (c & ~a) | (a & ~b & ~c);
    } else if constexpr (R == 0x1b) {
        return (c | ~b) & (~a | ~c);
    } else if constexpr (R == 0x1c) {
        return (b & ~a) | (a & ~b & ~c);
    } else if constexpr (R == 0x1d) {
        return (b | ~c) & (~a | ~b);
    } else if constexpr (R == 0x1e) {
        return a ^ b ^ c ^ (b & c);
    } else if constexpr (R == 0x1f) {
        return ~a | (~b & ~c);
    } else if constexpr (R == 0x20) {
        return a & c & ~b;
    } else if constexpr (R == 0x21) {
        return ~b & (a | ~c) & (c | ~a);
    } else if constexpr (R == 0x22) {
        return c & ~b;
    } else if constexpr (R == 0x23) {
        return ~b & (c | ~a);
    } else if constexpr (R == 0x24) {
        return (a & c & ~b) | (b & ~a & ~c);
    } else if constexpr (R == 0x25) {
        return (~a & ~c) | (a & c & ~b);
    } else if constexpr (R == 0x26) {
        return (c & ~b) | (b & ~a & ~c);
    } else if constexpr (R == 0x27) {
        return (c | ~a) & (~b | ~c);
    } else if constexpr (R == 0x28) {
        return (a & c) ^ (b & c);
    } else if constexpr (R == 0x29) {
        return a ^ b ^ c ^ ~T{} ^ (a & b) ^ (a & b & c);
    } else if constexpr (R == 0x2a) {
        return c & (~a | ~b);
    } else if constexpr (R == 0x2b) {
        return (c | ~a) & (c | ~b) & (~a | ~b);
    } else if constexpr (R == 0x2c) {
        return (b & ~a) | (a & c & ~b);
    } else if constexpr (R == 0x2d) {
        return a ^ c ^ ~T{} ^ (b & c);
    } else if constexpr (R == 0x2e) {
        return (b | c) & (~a | ~b);
    } else if constexpr (R == 0x2f) {
        return ~a | (c & ~b);
    } else if constexpr (R == 0x30) {
        return a & ~b;
    } else if constexpr (R == 0x31) {
        return ~b & (a | ~c);
    } else if constexpr (R == 0x32) {
        return ~b & (a | c);
    } else if constexpr (R == 0x33) {
        return ~b;
    } else if constexpr (R == 0x34) {
        return (a & ~b) | (b & ~a & ~c);
    } else if constexpr (R == 0x35) {
        return (a | ~c) & (~a | ~b);
    } else if constexpr (R == 0x36) {
        return a ^ b ^ c ^ (a & c);
    } else if constexpr (R == 0x37) {
        return ~b | (~a & ~c);
    } else if constexpr (R == 0x38) {
        return (a & ~b) | (b & c & ~a);
    } else if constexpr (R == 0x39) {
        return b ^ c ^ ~T{} ^ (a & c);
    } else if constexpr (R == 0x3a) {
        return (a | c) & (~a | ~b);
    } else if constexpr (R == 0x3b) {
        return ~b | (c & ~a);
    } else if constexpr (R == 0x3c) {
        return a ^ b;
    } else if constexpr (R == 0x3d) {
        return (~a | ~b) & (a | b | ~c);
    } else if constexpr (R == 0x3e) {
        return (a | b | c) & (~a | ~b);
    } else if constexpr (R == 0x3f) {
        return ~a | ~b;
    } else if constexpr (R == 0x40) {
        return a & b & ~c;
    } else if constexpr (R == 0x41) {
        return ~c & (a | ~b) & (b | ~a);
    } else if constexpr (R == 0x42) {
        return (a & b & ~c) | (c & ~a & ~b);
    } else if constexpr (R == 0x43) {
        return (~a & ~b) | (a & b & ~c);
    } else if constexpr (R == 0x44) {
        return b & ~c;
    } else if constexpr (R == 0x45) {
        return ~c & (b | ~a);
    } else if constexpr (R == 0x46) {
        return (b & ~c) | (c & ~a & ~b);
    } else if constexpr (R == 0x47) {
        return (b | ~a) & (~b | ~c);
    } else if constexpr (R == 0x48) {
        return (a & b) ^ (b & c);
    } else if constexpr (R == 0x49) {
        return a ^ b ^ c ^ ~T{} ^ (a & c) ^ (a & b & c);
    } else if constexpr (R == 0x4a) {
        return (c & ~a) | (a & b & ~c);
    } else if constexpr (R == 0x4b) {
        return a ^ b ^ ~T{} ^ (b & c);
    } else if constexpr (R == 0x4c) {
        return b & (~a | ~c);
    } else if constexpr (R == 0x4d) {
        return (b | ~a) & (b | ~c) & (~a | ~c);
    } else if constexpr (R == 0x4e) {
        return (b | c) & (~a | ~c);
    } else if constexpr (R == 0x4f) {
        return ~a | (b & ~c);
    } else if constexpr (R == 0x50) {
        return a & ~c;
    } else if constexpr (R == 0x51) {
        return ~c & (a | ~b);
    } else if constexpr (R == 0x52) {
        return (a & ~c) | (c & ~a & ~b);
    } else if constexpr (R == 0x53) {
        return (a | ~b) & (~a | ~c);
    } else if constexpr (R == 0x54) {
        return ~c & (a | b);
    } else if constexpr (R == 0x55) {
        return ~c;
    } else if constexpr (R == 0x56) {
        return a ^ b ^ c ^ (a & b);
    } else if constexpr (R == 0x57) {
        return ~c | (~a & ~b);
    } else if constexpr (R == 0x58) {
        return (a & ~c) | (b & c & ~a);
    } else if constexpr (R == 0x59) {
        return b ^ c ^ ~T{} ^ (a & b);
    } else if constexpr (R == 0x5a) {
        return a ^ c;
    } else if constexpr (R == 0x5b) {
        return (~a | ~c) & (a | c | ~b);
    } else if constexpr (R == 0x5c) {
        return (a | b) & (~a | ~c);
    } else if constexpr (R == 0x5d) {
        return ~c | (b & ~a);
    } else if constexpr (R == 0x5e) {
        return (a | b | c) & (~a | ~c);
    } else if constexpr (R == 0x5f) {
        return ~a | ~c;
    } else if constexpr (R == 0x60) {
        return (a & b) ^ (a & c);
    } else if constexpr (R == 0x61) {
        return a ^ b ^ c ^ ~T{} ^ (b & c) ^ (a & b & c);
    } else if constexpr (R == 0x62) {
        return (c & ~b) | (a & b & ~c);
    } else if constexpr (R == 0x63) {
        return a ^ b ^ ~T{} ^ (a & c);
    } else if constexpr (R == 0x64) {
        return (b & ~c) | (a & c & ~b);
    } else if constexpr (R == 0x65) {
        return a ^ c ^ ~T{} ^ (a & b);
    } else if constexpr (R == 0x66) {
        return b ^ c;
    } else if constexpr (R == 0x67) {
        return (~b | ~c) & (b | c | ~a);
    } else if constexpr (R == 0x68) {
        return (a & b) ^ (a & c) ^ (b & c) ^ (a & b & c);
    } else if constexpr (R == 0x69) {
        return a ^ b ^ c ^ ~T{};
    } else if constexpr (R == 0x6a) {
        return c ^ (a & b);
    } else if constexpr (R == 0x6b) {
        return (a | c | ~b) & (b | c | ~a) & (~a | ~b | ~c);
    } else if constexpr (R == 0x6c) {
        return b ^ (a & c);
    } else if constexpr (R == 0x6d) {
        return (a | b | ~c) & (b | c | ~a) & (~a | ~b | ~c);
    } else if constexpr (R == 0x6e) {
        return (b | c) & (~a | ~b | ~c);
    } else if constexpr (R == 0x6f) {
        return ~a | (b & ~c) | (c & ~b);
    } else if constexpr (R == 0x70) {
        return a & (~b | ~c);
    } else if constexpr (R == 0x71) {
        return (a | ~b) & (a | ~c) & (~b | ~c);
    } else if constexpr (R == 0x72) {
        return (a | c) & (~b | ~c);
    } else if constexpr (R == 0x73) {
        return ~b | (a & ~c);
    } else if constexpr (R == 0x74) {
        return (a | b) & (~b | ~c);
    } else if constexpr (R == 0x75) {
        return ~c | (a & ~b);
    } else if constexpr (R == 0x76) {
        return (a | b | c) & (~b | ~c);
    } else if constexpr (R == 0x77) {
        return ~b | ~c;
    } else if constexpr (R == 0x78) {
        return a ^ (b & c);
    } else if constexpr (R == 0x79) {
        return (a | b | ~c) & (a | c | ~b) & (~a | ~b | ~c);
    } else if constexpr (R == 0x7a) {
        return (a | c) & (~a | ~b | ~c);
    } else if constexpr (R == 0x7b) {
        return ~b | (a & ~c) | (c & ~a);
    } else if constexpr (R == 0x7c) {
        return (a | b) & (~a | ~b | ~c);
    } else if constexpr (R == 0x7d) {
        return ~c | (a & ~b) | (b & ~a);
    } else if constexpr (R == 0x7e) {
        return (a | b | c) & (~a | ~b | ~c);
    } else if constexpr (R == 0x7f) {
        return ~a | ~b | ~c;
    } else if constexpr (R == 0x80) {
        return a & b & c;
    } else if constexpr (R == 0x81) {
        return (a & b & c) | (~a & ~b & ~c);
    } else if constexpr (R == 0x82) {
        return c ^ (a & c) ^ (b & c);
    } else if constexpr (R == 0x83) {
        return (a & b & c) | (~a & ~b);
    } else if constexpr (R == 0x84) {
        return b ^ (a & b) ^ (b & c);
    } else if constexpr (R == 0x85) {
        return (a & b & c) | (~a & ~c);
    } else if constexpr (R == 0x86) {
        return b ^ c ^ (a & b) ^ (a & c) ^ (a & b & c);
    } else if constexpr (R == 0x87) {
        return a ^ ~T{} ^ (b & c);
    } else if constexpr (R == 0x88) {
        return b & c;
    } else if constexpr (R == 0x89) {
        return (b & c) | (~a & ~b & ~c);
    } else if constexpr (R == 0x8a) {
        return c & (b | ~a);
    } else if constexpr (R == 0x8b) {
        return (b | ~a) & (c | ~b);
    } else if constexpr (R == 0x8c) {
        return b & (c | ~a);
    } else if constexpr (R == 0x8d) {
        return (b | ~c) & (c | ~a);
    } else if constexpr (R == 0x8e) {
        return (b | c) & (b | ~a) & (c | ~a);
    } else if constexpr (R == 0x8f) {
        return ~a | (b & c);
    } else if constexpr (R == 0x90) {
        return a ^ (a & b) ^ (a & c);
    } else if constexpr (R == 0x91) {
        return (a & b & c) | (~b & ~c);
    } else if constexpr (R == 0x92) {
        return a ^ c ^ (a & b) ^ (b & c) ^ (a & b & c);
    } else if constexpr (R == 0x93) {
        return b ^ ~T{} ^ (a & c);
    } else if constexpr (R == 0x94) {
        return a ^ b ^ (a & c) ^ (b & c) ^ (a & b & c);
    } else if constexpr (R == 0x95) {
        return c ^ ~T{} ^ (a & b);
    } else if constexpr (R == 0x96) {
        return a ^ b ^ c;
    } else if constexpr (R == 0x97) {
        return (a | ~b | ~c) & (b | ~a | ~c) & (c | ~a | ~b);
    } else if constexpr (R == 0x98) {
        return (b & c) | (a & ~b & ~c);
    } else if constexpr (R == 0x99) {
        return b ^ c ^ ~T{};
    } else if constexpr (R == 0x9a) {
        return a ^ c ^ (a & b);
    } else if constexpr (R == 0x9b) {
        return (c | ~b) & (b | ~a | ~c);
    } else if constexpr (R == 0x9c) {
        return a ^ b ^ (a & c);
    } else if constexpr (R == 0x9d) {
        return (b | ~c) & (c | ~a | ~b);
    } else if constexpr (R == 0x9e) {
        return a ^ b ^ c ^ (b & c) ^ (a & b & c);
    } else if constexpr (R == 0x9f) {
        return ~a | (b & c) | (~b & ~c);
    } else if constexpr (R == 0xa0) {
        return a & c;
    } else if constexpr (R == 0xa1) {
        return (a & c) | (~a & ~b & ~c);
    } else if constexpr (R == 0xa2) {
        return c & (a | ~b);
    } else if constexpr (R == 0xa3) {
        return (a | ~b) & (c | ~a);
    } else if constexpr (R == 0xa4) {
        return (a & c) | (b & ~a & ~c);
    } else if constexpr (R == 0xa5) {
        return a ^ c ^ ~T{};
    } else if constexpr (R == 0xa6) {
        return b ^ c ^ (a & b);
    } else if constexpr (R == 0xa7) {
        return (c | ~a) & (a | ~b | ~c);
    } else if constexpr (R == 0xa8) {
        return c & (a | b);
    } else if constexpr (R == 0xa9) {
        return a ^ b ^ c ^ ~T{} ^ (a & b);
    } else if constexpr (R == 0xaa) {
        return c;
    } else if constexpr (R == 0xab) {
        return c | (~a & ~b);
    } else if constexpr (R == 0xac) {
        return (a | b) & (c | ~a);
    } else if constexpr (R == 0xad) {
        return (c | ~a) & (a | b | ~c);
    } else if constexpr (R == 0xae) {
        return c | (b & ~a);
    } else if constexpr (R == 0xaf) {
        return c | ~a;
    } else if constexpr (R == 0xb0) {
        return a & (c | ~b);
    } else if constexpr (R == 0xb1) {
        return (a | ~c) & (c | ~b);
    } else if constexpr (R == 0xb2) {
        return (a | c) & (a | ~b) & (c | ~b);
    } else if constexpr (R == 0xb3) {
        return ~b | (a & c);
    } else if constexpr (R == 0xb4) {
        return a ^ b ^ (b & c);
    } else if constexpr (R == 0xb5) {
        return (a | ~c) & (c | ~a | ~b);
    } else if constexpr (R == 0xb6) {
        return a ^ b ^ c ^ (a & c) ^ (a & b & c);
    } else if constexpr (R == 0xb7) {
        return ~b | (a & c) | (~a & ~c);
    } else if constexpr (R == 0xb8) {
        return (a | b) & (c | ~b);
    } else if constexpr (R == 0xb9) {
        return (c | ~b) & (a | b | ~c);
    } else if constexpr (R == 0xba) {
        return c | (a & ~b);
    } else if constexpr (R == 0xbb) {
        return c | ~b;
    } else if constexpr (R == 0xbc) {
        return a ^ b ^ (a & b & c);
    } else if constexpr (R == 0xbd) {
        return (a | b | ~c) & (c | ~a | ~b);
    } else if constexpr (R == 0xbe) {
        return c | (a & ~b) | (b & ~a);
    } else if constexpr (R == 0xbf) {
        return c | ~a | ~b;
    } else if constexpr (R == 0xc0) {
        return a & b;
    } else if constexpr (R == 0xc1) {
        return (a & b) | (~a & ~b & ~c);
    } else if constexpr (R == 0xc2) {
        return (a & b) | (c & ~a & ~b);
    } else if constexpr (R == 0xc3) {
        return a ^ b ^ ~T{};
    } else if constexpr (R == 0xc4) {
        return b & (a | ~c);
    } else if constexpr (R == 0xc5) {
        return (a | ~c) & (b | ~a);
    } else if constexpr (R == 0xc6) {
        return b ^ c ^ (a & c);
    } else if constexpr (R == 0xc7) {
        return (b | ~a) & (a | ~b | ~c);
    } else if constexpr (R == 0xc8) {
        return b & (a | c);
    } else if constexpr (R == 0xc9) {
        return a ^ b ^ c ^ ~T{} ^ (a & c);
    } else if constexpr (R == 0xca) {
        return (a | c) & (b | ~a);
    } else if constexpr (R == 0xcb) {
        return (b | ~a) & (a | c | ~b);
    } else if constexpr (R == 0xcc) {
        return b;
    } else if constexpr (R == 0xcd) {
        return b | (~a & ~c);
    } else if constexpr (R == 0xce) {
        return b | (c & ~a);
    } else if constexpr (R == 0xcf) {
        return b | ~a;
    } else if constexpr (R == 0xd0) {
        return a & (b | ~c);
    } else if constexpr (R == 0xd1) {
        return (a | ~b) & (b | ~c);
    } else if constexpr (R == 0xd2) {
        return a ^ c ^ (b & c);
    } else if constexpr (R == 0xd3) {
        return (a | ~b) & (b | ~a | ~c);
    } else if constexpr (R == 0xd4) {
        return (a | b) & (a | ~c) & (b | ~c);
    } else if constexpr (R == 0xd5) {
        return ~c | (a & b);
    } else if constexpr (R == 0xd6) {
        return a ^ b ^ c ^ (a & b) ^ (a & b & c);
    } else if constexpr (R == 0xd7) {
        return ~c | (a & b) | (~a & ~b);
    } else if constexpr (R == 0xd8) {
        return (a | c) & (b | ~c);
    } else if constexpr (R == 0xd9) {
        return (b | ~c) & (a | c | ~b);
    } else if constexpr (R == 0xda) {
        return a ^ c ^ (a & b & c);
    } else if constexpr (R == 0xdb) {
        return (a | c | ~b) & (b | ~a | ~c);
    } else if constexpr (R == 0xdc) {
        return b | (a & ~c);
    } else if constexpr (R == 0xdd) {
        return b | ~c;
    } else if constexpr (R == 0xde) {
        return b | (a & ~c) | (c & ~a);
    } else if constexpr (R == 0xdf) {
        return b | ~a | ~c;
    } else if constexpr (R == 0xe0) {
        return a & (b | c);
    } else if constexpr (R == 0xe1) {
        return a ^ b ^ c ^ ~T{} ^ (b & c);
    } else if constexpr (R == 0xe2) {
        return (b | c) & (a | ~b);
    } else if constexpr (R == 0xe3) {
        return (a | ~b) & (b | c | ~a);
    } else if constexpr (R == 0xe4) {
        return (b | c) & (a | ~c);
    } else if constexpr (R == 0xe5) {
        return (a | ~c) & (b | c | ~a);
    } else if constexpr (R == 0xe6) {
        return b ^ c ^ (a & b & c);
    } else if constexpr (R == 0xe7) {
        return (b | c | ~a) & (a | ~b | ~c);
    } else if constexpr (R == 0xe8) {
        return (a | b) & (a | c) & (b | c);
    } else if constexpr (R == 0xe9) {
        return a ^ b ^ c ^ ~T{} ^ (a & b & c);
    } else if constexpr (R == 0xea) {
        return c | (a & b);
    } else if constexpr (R == 0xeb) {
        return c | (a & b) | (~a & ~b);
    } else if constexpr (R == 0xec) {
        return b | (a & c);
    } else if constexpr (R == 0xed) {
        return b | (a & c) | (~a & ~c);
    } else if constexpr (R == 0xee) {
        return b | c;
    } else if constexpr (R == 0xef) {
        return b | c | ~a;
    } else if constexpr (R == 0xf0) {
        return a;
    } else if constexpr (R == 0xf1) {
        return a | (~b & ~c);
    } else if constexpr (R == 0xf2) {
        return a | (c & ~b);
    } else if constexpr (R == 0xf3) {
        return a | ~b;
    } else if constexpr (R == 0xf4) {
        return a | (b & ~c);
    } else if constexpr (R == 0xf5) {
        return a | ~c;
    } else if constexpr (R == 0xf6) {
        return a | (b & ~c) | (c & ~b);
    } else if constexpr (R == 0xf7) {
        return a | ~b | ~c;
    } else if constexpr (R == 0xf8) {
        return a | (b & c);
    } else if constexpr (R == 0xf9) {
        return a | (b & c) | (~b & ~c);
    } else if constexpr (R == 0xfa) {
        return a | c;
    } else if constexpr (R == 0xfb) {
        return a | c | ~b;
    } else if constexpr (R == 0xfc) {
        return a | b;
    } else if constexpr (R == 0xfd) {
        return a | b | ~c;
    } else if constexpr (R == 0xfe) {
        return a | b | c;
    } else if constexpr (R == 0xff) {
        return ~T{};
    }
}

template <size_t N, typename T=std::bitset<N>>
auto ternarylogic_impl2(uint8_t R, T const& a, T const& b, T const& c) -> T {
    switch(R) {
    case 0x00: return T{};
    case 0x01: return ~a & ~b & ~c;
    case 0x02: return c & ~a & ~b;
    case 0x03: return ~a & ~b;
    case 0x04: return b & ~a & ~c;
    case 0x05: return ~a & ~c;
    case 0x06: return ~a & (b | c) & (~b | ~c);
    case 0x07: return ~a & (~b | ~c);
    case 0x08: return b & c & ~a;
    case 0x09: return ~a & (b | ~c) & (c | ~b);
    case 0x0a: return c & ~a;
    case 0x0b: return ~a & (c | ~b);
    case 0x0c: return b & ~a;
    case 0x0d: return ~a & (b | ~c);
    case 0x0e: return ~a & (b | c);
    case 0x0f: return ~a;
    case 0x10: return a & ~b & ~c;
    case 0x11: return ~b & ~c;
    case 0x12: return ~b & (a | c) & (~a | ~c);
    case 0x13: return ~b & (~a | ~c);
    case 0x14: return ~c & (a | b) & (~a | ~b);
    case 0x15: return ~c & (~a | ~b);
    case 0x16: return a ^ b ^ c ^ (a & b & c);
    case 0x17: return (~a | ~b) & (~a | ~c) & (~b | ~c);
    case 0x18: return (b & c & ~a) | (a & ~b & ~c);
    case 0x19: return (~b & ~c) | (b & c & ~a);
    case 0x1a: return (c & ~a) | (a & ~b & ~c);
    case 0x1b: return (c | ~b) & (~a | ~c);
    case 0x1c: return (b & ~a) | (a & ~b & ~c);
    case 0x1d: return (b | ~c) & (~a | ~b);
    case 0x1e: return a ^ b ^ c ^ (b & c);
    case 0x1f: return ~a | (~b & ~c);
    case 0x20: return a & c & ~b;
    case 0x21: return ~b & (a | ~c) & (c | ~a);
    case 0x22: return c & ~b;
    case 0x23: return ~b & (c | ~a);
    case 0x24: return (a & c & ~b) | (b & ~a & ~c);
    case 0x25: return (~a & ~c) | (a & c & ~b);
    case 0x26: return (c & ~b) | (b & ~a & ~c);
    case 0x27: return (c | ~a) & (~b | ~c);
    case 0x28: return (a & c) ^ (b & c);
    case 0x29: return a ^ b ^ c ^ ~T{} ^ (a & b) ^ (a & b & c);
    case 0x2a: return c & (~a | ~b);
    case 0x2b: return (c | ~a) & (c | ~b) & (~a | ~b);
    case 0x2c: return (b & ~a) | (a & c & ~b);
    case 0x2d: return a ^ c ^ ~T{} ^ (b & c);
    case 0x2e: return (b | c) & (~a | ~b);
    case 0x2f: return ~a | (c & ~b);
    case 0x30: return a & ~b;
    case 0x31: return ~b & (a | ~c);
    case 0x32: return ~b & (a | c);
    case 0x33: return ~b;
    case 0x34: return (a & ~b) | (b & ~a & ~c);
    case 0x35: return (a | ~c) & (~a | ~b);
    case 0x36: return a ^ b ^ c ^ (a & c);
    case 0x37: return ~b | (~a & ~c);
    case 0x38: return (a & ~b) | (b & c & ~a);
    case 0x39: return b ^ c ^ ~T{} ^ (a & c);
    case 0x3a: return (a | c) & (~a | ~b);
    case 0x3b: return ~b | (c & ~a);
    case 0x3c: return a ^ b;
    case 0x3d: return (~a | ~b) & (a | b | ~c);
    case 0x3e: return (a | b | c) & (~a | ~b);
    case 0x3f: return ~a | ~b;
    case 0x40: return a & b & ~c;
    case 0x41: return ~c & (a | ~b) & (b | ~a);
    case 0x42: return (a & b & ~c) | (c & ~a & ~b);
    case 0x43: return (~a & ~b) | (a & b & ~c);
    case 0x44: return b & ~c;
    case 0x45: return ~c & (b | ~a);
    case 0x46: return (b & ~c) | (c & ~a & ~b);
    case 0x47: return (b | ~a) & (~b | ~c);
    case 0x48: return (a & b) ^ (b & c);
    case 0x49: return a ^ b ^ c ^ ~T{} ^ (a & c) ^ (a & b & c);
    case 0x4a: return (c & ~a) | (a & b & ~c);
    case 0x4b: return a ^ b ^ ~T{} ^ (b & c);
    case 0x4c: return b & (~a | ~c);
    case 0x4d: return (b | ~a) & (b | ~c) & (~a | ~c);
    case 0x4e: return (b | c) & (~a | ~c);
    case 0x4f: return ~a | (b & ~c);
    case 0x50: return a & ~c;
    case 0x51: return ~c & (a | ~b);
    case 0x52: return (a & ~c) | (c & ~a & ~b);
    case 0x53: return (a | ~b) & (~a | ~c);
    case 0x54: return ~c & (a | b);
    case 0x55: return ~c;
    case 0x56: return a ^ b ^ c ^ (a & b);
    case 0x57: return ~c | (~a & ~b);
    case 0x58: return (a & ~c) | (b & c & ~a);
    case 0x59: return b ^ c ^ ~T{} ^ (a & b);
    case 0x5a: return a ^ c;
    case 0x5b: return (~a | ~c) & (a | c | ~b);
    case 0x5c: return (a | b) & (~a | ~c);
    case 0x5d: return ~c | (b & ~a);
    case 0x5e: return (a | b | c) & (~a | ~c);
    case 0x5f: return ~a | ~c;
    case 0x60: return (a & b) ^ (a & c);
    case 0x61: return a ^ b ^ c ^ ~T{} ^ (b & c) ^ (a & b & c);
    case 0x62: return (c & ~b) | (a & b & ~c);
    case 0x63: return a ^ b ^ ~T{} ^ (a & c);
    case 0x64: return (b & ~c) | (a & c & ~b);
    case 0x65: return a ^ c ^ ~T{} ^ (a & b);
    case 0x66: return b ^ c;
    case 0x67: return (~b | ~c) & (b | c | ~a);
    case 0x68: return (a & b) ^ (a & c) ^ (b & c) ^ (a & b & c);
    case 0x69: return a ^ b ^ c ^ ~T{};
    case 0x6a: return c ^ (a & b);
    case 0x6b: return (a | c | ~b) & (b | c | ~a) & (~a | ~b | ~c);
    case 0x6c: return b ^ (a & c);
    case 0x6d: return (a | b | ~c) & (b | c | ~a) & (~a | ~b | ~c);
    case 0x6e: return (b | c) & (~a | ~b | ~c);
    case 0x6f: return ~a | (b & ~c) | (c & ~b);
    case 0x70: return a & (~b | ~c);
    case 0x71: return (a | ~b) & (a | ~c) & (~b | ~c);
    case 0x72: return (a | c) & (~b | ~c);
    case 0x73: return ~b | (a & ~c);
    case 0x74: return (a | b) & (~b | ~c);
    case 0x75: return ~c | (a & ~b);
    case 0x76: return (a | b | c) & (~b | ~c);
    case 0x77: return ~b | ~c;
    case 0x78: return a ^ (b & c);
    case 0x79: return (a | b | ~c) & (a | c | ~b) & (~a | ~b | ~c);
    case 0x7a: return (a | c) & (~a | ~b | ~c);
    case 0x7b: return ~b | (a & ~c) | (c & ~a);
    case 0x7c: return (a | b) & (~a | ~b | ~c);
    case 0x7d: return ~c | (a & ~b) | (b & ~a);
    case 0x7e: return (a | b | c) & (~a | ~b | ~c);
    case 0x7f: return ~a | ~b | ~c;
    case 0x80: return a & b & c;
    case 0x81: return (a & b & c) | (~a & ~b & ~c);
    case 0x82: return c ^ (a & c) ^ (b & c);
    case 0x83: return (a & b & c) | (~a & ~b);
    case 0x84: return b ^ (a & b) ^ (b & c);
    case 0x85: return (a & b & c) | (~a & ~c);
    case 0x86: return b ^ c ^ (a & b) ^ (a & c) ^ (a & b & c);
    case 0x87: return a ^ ~T{} ^ (b & c);
    case 0x88: return b & c;
    case 0x89: return (b & c) | (~a & ~b & ~c);
    case 0x8a: return c & (b | ~a);
    case 0x8b: return (b | ~a) & (c | ~b);
    case 0x8c: return b & (c | ~a);
    case 0x8d: return (b | ~c) & (c | ~a);
    case 0x8e: return (b | c) & (b | ~a) & (c | ~a);
    case 0x8f: return ~a | (b & c);
    case 0x90: return a ^ (a & b) ^ (a & c);
    case 0x91: return (a & b & c) | (~b & ~c);
    case 0x92: return a ^ c ^ (a & b) ^ (b & c) ^ (a & b & c);
    case 0x93: return b ^ ~T{} ^ (a & c);
    case 0x94: return a ^ b ^ (a & c) ^ (b & c) ^ (a & b & c);
    case 0x95: return c ^ ~T{} ^ (a & b);
    case 0x96: return a ^ b ^ c;
    case 0x97: return (a | ~b | ~c) & (b | ~a | ~c) & (c | ~a | ~b);
    case 0x98: return (b & c) | (a & ~b & ~c);
    case 0x99: return b ^ c ^ ~T{};
    case 0x9a: return a ^ c ^ (a & b);
    case 0x9b: return (c | ~b) & (b | ~a | ~c);
    case 0x9c: return a ^ b ^ (a & c);
    case 0x9d: return (b | ~c) & (c | ~a | ~b);
    case 0x9e: return a ^ b ^ c ^ (b & c) ^ (a & b & c);
    case 0x9f: return ~a | (b & c) | (~b & ~c);
    case 0xa0: return a & c;
    case 0xa1: return (a & c) | (~a & ~b & ~c);
    case 0xa2: return c & (a | ~b);
    case 0xa3: return (a | ~b) & (c | ~a);
    case 0xa4: return (a & c) | (b & ~a & ~c);
    case 0xa5: return a ^ c ^ ~T{};
    case 0xa6: return b ^ c ^ (a & b);
    case 0xa7: return (c | ~a) & (a | ~b | ~c);
    case 0xa8: return c & (a | b);
    case 0xa9: return a ^ b ^ c ^ ~T{} ^ (a & b);
    case 0xaa: return c;
    case 0xab: return c | (~a & ~b);
    case 0xac: return (a | b) & (c | ~a);
    case 0xad: return (c | ~a) & (a | b | ~c);
    case 0xae: return c | (b & ~a);
    case 0xaf: return c | ~a;
    case 0xb0: return a & (c | ~b);
    case 0xb1: return (a | ~c) & (c | ~b);
    case 0xb2: return (a | c) & (a | ~b) & (c | ~b);
    case 0xb3: return ~b | (a & c);
    case 0xb4: return a ^ b ^ (b & c);
    case 0xb5: return (a | ~c) & (c | ~a | ~b);
    case 0xb6: return a ^ b ^ c ^ (a & c) ^ (a & b & c);
    case 0xb7: return ~b | (a & c) | (~a & ~c);
    case 0xb8: return (a | b) & (c | ~b);
    case 0xb9: return (c | ~b) & (a | b | ~c);
    case 0xba: return c | (a & ~b);
    case 0xbb: return c | ~b;
    case 0xbc: return a ^ b ^ (a & b & c);
    case 0xbd: return (a | b | ~c) & (c | ~a | ~b);
    case 0xbe: return c | (a & ~b) | (b & ~a);
    case 0xbf: return c | ~a | ~b;
    case 0xc0: return a & b;
    case 0xc1: return (a & b) | (~a & ~b & ~c);
    case 0xc2: return (a & b) | (c & ~a & ~b);
    case 0xc3: return a ^ b ^ ~T{};
    case 0xc4: return b & (a | ~c);
    case 0xc5: return (a | ~c) & (b | ~a);
    case 0xc6: return b ^ c ^ (a & c);
    case 0xc7: return (b | ~a) & (a | ~b | ~c);
    case 0xc8: return b & (a | c);
    case 0xc9: return a ^ b ^ c ^ ~T{} ^ (a & c);
    case 0xca: return (a | c) & (b | ~a);
    case 0xcb: return (b | ~a) & (a | c | ~b);
    case 0xcc: return b;
    case 0xcd: return b | (~a & ~c);
    case 0xce: return b | (c & ~a);
    case 0xcf: return b | ~a;
    case 0xd0: return a & (b | ~c);
    case 0xd1: return (a | ~b) & (b | ~c);
    case 0xd2: return a ^ c ^ (b & c);
    case 0xd3: return (a | ~b) & (b | ~a | ~c);
    case 0xd4: return (a | b) & (a | ~c) & (b | ~c);
    case 0xd5: return ~c | (a & b);
    case 0xd6: return a ^ b ^ c ^ (a & b) ^ (a & b & c);
    case 0xd7: return ~c | (a & b) | (~a & ~b);
    case 0xd8: return (a | c) & (b | ~c);
    case 0xd9: return (b | ~c) & (a | c | ~b);
    case 0xda: return a ^ c ^ (a & b & c);
    case 0xdb: return (a | c | ~b) & (b | ~a | ~c);
    case 0xdc: return b | (a & ~c);
    case 0xdd: return b | ~c;
    case 0xde: return b | (a & ~c) | (c & ~a);
    case 0xdf: return b | ~a | ~c;
    case 0xe0: return a & (b | c);
    case 0xe1: return a ^ b ^ c ^ ~T{} ^ (b & c);
    case 0xe2: return (b | c) & (a | ~b);
    case 0xe3: return (a | ~b) & (b | c | ~a);
    case 0xe4: return (b | c) & (a | ~c);
    case 0xe5: return (a | ~c) & (b | c | ~a);
    case 0xe6: return b ^ c ^ (a & b & c);
    case 0xe7: return (b | c | ~a) & (a | ~b | ~c);
    case 0xe8: return (a | b) & (a | c) & (b | c);
    case 0xe9: return a ^ b ^ c ^ ~T{} ^ (a & b & c);
    case 0xea: return c | (a & b);
    case 0xeb: return c | (a & b) | (~a & ~b);
    case 0xec: return b | (a & c);
    case 0xed: return b | (a & c) | (~a & ~c);
    case 0xee: return b | c;
    case 0xef: return b | c | ~a;
    case 0xf0: return a;
    case 0xf1: return a | (~b & ~c);
    case 0xf2: return a | (c & ~b);
    case 0xf3: return a | ~b;
    case 0xf4: return a | (b & ~c);
    case 0xf5: return a | ~c;
    case 0xf6: return a | (b & ~c) | (c & ~b);
    case 0xf7: return a | ~b | ~c;
    case 0xf8: return a | (b & c);
    case 0xf9: return a | (b & c) | (~b & ~c);
    case 0xfa: return a | c;
    case 0xfb: return a | c | ~b;
    case 0xfc: return a | b;
    case 0xfd: return a | b | ~c;
    case 0xfe: return a | b | c;
    case 0xff: return ~T{};
    }
}


template <size_t Begin, size_t End, typename CB>
void for_constexpr(CB cb) {
    if constexpr (Begin != End) {
        cb.template operator()<Begin>();
        for_constexpr<Begin+1, End>(cb);
    }
}

template <size_t N, typename T=std::bitset<N>>
static std::array<std::function<T(T const&, T const&, T const&)>, 256> lut_ternarylogic = []() {
    auto r = std::array<std::function<T(T const&, T const&, T const&)>, 256>{};
    for_constexpr<0, 256>([&]<size_t I>() {
        r[I] = &ternarylogic_impl<I, N, T>;
    });
    return r;
}();

template <size_t N, typename T=std::bitset<N>>
static std::array<T(*)(T const&, T const&, T const&), 256> lut_ternarylogic2 = []() {
    auto r = std::array<T(*)(T const&, T const&, T const&), 256>{};
    for_constexpr<0, 256>([&]<size_t I>() {
        r[I] = &ternarylogic_impl<I, N, T>;
    });
    return r;
}();



template <size_t N, typename T=std::bitset<N>>
auto ternarylogic(size_t R, T const& a, T const& b, T const& c) -> T {
//    return ternarylogic_impl2<N, T>(static_cast<uint8_t>(R), a, b, c);
    return lut_ternarylogic2<N, T>[R](a, b, c);
}

