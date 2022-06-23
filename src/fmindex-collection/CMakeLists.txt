cmake_minimum_required (VERSION 3.8)

project(fmindex-collection)

# fmindex_collection library
add_library(${PROJECT_NAME} INTERFACE)

option(FMC_USE_SDSL "Activate occ tables that use SDSL as a backend" ${ROOT_PROJECT})

target_include_directories(${PROJECT_NAME}
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(${PROJECT_NAME}
    INTERFACE
)

if (FMC_USE_SDSL)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DFMC_USE_SDSL)
    target_include_directories(${PROJECT_NAME}
        SYSTEM INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/sdsl-lite/include
    )
endif()
