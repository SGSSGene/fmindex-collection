cmake_minimum_required (VERSION 3.8)

project(test_fmindex-collection)

# fmindex-collectionunittests
add_executable(${PROJECT_NAME}
    checkBiFMIndex.cpp
    checkDenseBiFMIndex.cpp
    checkFMIndex.cpp
    checkReverseFMIndex.cpp
    checkOccTables.cpp
    checkSearchPseudo.cpp
    checkSearchBacktracking.cpp
    main.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2
    fmindex-collection
)


target_include_directories(${PROJECT_NAME} PRIVATE ${DIVSUFSORT64_INCLUDE_DIRS})
target_compile_options(${PROJECT_NAME} PRIVATE ${DIVSUFSORT64_CFLAGS})
target_link_libraries(${PROJECT_NAME} ${DIVSUFSORT64_LIBRARIES})
