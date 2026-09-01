# BWA C++26 Core Library
add_library(bwa-core STATIC
    core/arena.cpp
)

target_include_directories(bwa-core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(bwa-core PUBLIC cxx_std_26)

# BWA Index Library
add_library(bwa-index STATIC
    index/fm_index.cpp
)

target_link_libraries(bwa-index PUBLIC bwa-core)
target_include_directories(bwa-index PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

# BWA Align Library
add_library(bwa-align STATIC
    align/sw.cpp
    align/mem.cpp
)

target_link_libraries(bwa-align PUBLIC bwa-core bwa-index)
target_include_directories(bwa-align PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

# BWA IO Library
add_library(bwa-io STATIC
    io/seq_io.cpp
)

target_link_libraries(bwa-io PUBLIC bwa-core)
target_include_directories(bwa-io PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)
find_package(ZLIB REQUIRED)
target_link_libraries(bwa-io PUBLIC ZLIB::ZLIB)

# BWA Pipeline Library
add_library(bwa-pipeline STATIC
    pipeline.cpp
)

target_link_libraries(bwa-pipeline PUBLIC bwa-core bwa-index bwa-align bwa-io)
target_include_directories(bwa-pipeline PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
)

# BWA Utils Library
add_library(bwa-utils STATIC
    utils/utils.cpp
)

target_link_libraries(bwa-utils PUBLIC bwa-core)

# Main executable
add_executable(bwa
    main.cpp
)

target_link_libraries(bwa PRIVATE
    bwa-core bwa-index bwa-align bwa-io bwa-pipeline bwa-utils
    Threads::Threads
    ZLIB::ZLIB
)

# Tests
add_executable(bwa-tests
    ../tests/unit/test_core.cpp
    ../tests/unit/test_index.cpp
    ../tests/unit/test_align.cpp
)

target_link_libraries(bwa-tests PRIVATE
    bwa-core bwa-index bwa-align bwa-io bwa-pipeline bwa-utils
    Threads::Threads
    ZLIB::ZLIB
)

add_test(NAME bwa_unit_tests COMMAND bwa-tests)

# Benchmarks
add_executable(bwa-bench
    ../benchmarks/bench_core.cpp
    ../benchmarks/bench_align.cpp
)

target_link_libraries(bwa-bench PRIVATE
    bwa-core bwa-index bwa-align bwa-io bwa-pipeline bwa-utils
    Threads::Threads
    ZLIB::ZLIB
)

# Install
install(TARGETS bwa
    RUNTIME DESTINATION bin
)

install(TARGETS bwa-core bwa-index bwa-align bwa-io bwa-pipeline bwa-utils
    EXPORT bwa-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY ../include/bwa
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)