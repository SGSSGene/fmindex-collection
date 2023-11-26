# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.14)

# Check every header compiles on its own.
# Validates that headers have all required includes.

set(target_name "test_header_fmc")

# Path to the library relative from this CMakeLists.txt file
set(library_path "../fmindex-collection")

file(GLOB_RECURSE header_files LIST_DIRECTORIES false RELATIVE ${CMAKE_CURRENT_LIST_DIR}/.. ${library_path}/*.h)
set(files)
foreach (header ${header_files})
    string (REPLACE "/" "__" test_name "${header}")
    set(cpp_file test_header_${target_name}_${test_name}.cpp)
    configure_file(test.cpp.in ${cpp_file})
    set(files ${files} ${cpp_file})
endforeach()

# Path to the library relative from this CMakeLists.txt file
set(library_path "../search_schemes")

file(GLOB_RECURSE header_files LIST_DIRECTORIES false RELATIVE ${CMAKE_CURRENT_LIST_DIR}/.. ${library_path}/*.h)
foreach (header ${header_files})
    string (REPLACE "/" "__" test_name "${header}")
    set(cpp_file test_header_${target_name}_${test_name}.cpp)
    configure_file(test.cpp.in ${cpp_file})
    set(files ${files} ${cpp_file})
endforeach()

add_library(${target_name} ${files})
target_link_libraries(${target_name} fmindex-collection::fmindex-collection)
