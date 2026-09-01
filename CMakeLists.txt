cmake_minimum_required(VERSION 3.28)
project(bwa-cpp26 VERSION 0.1.0 LANGUAGES CXX)

# C++26 required
set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Strict warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow -Wnull-dereference -Wdouble-promotion -Wformat=2 -Wundef -Wcast-align -Wcast-qual -Wdisabled-optimization -Wfloat-equal -Wlogical-op -Wmissing-declarations -Wmissing-format-attribute -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wstrict-null-sentinel -Wsuggest-override -Wswitch-enum -Wswitch-default -Wuseless-cast -Wzero-as-null-pointer-constant")
    # Sanitizers (enable via -DSANITIZE=ON)
    option(SANITIZE "Enable sanitizers" ON)
    if(SANITIZE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined,leak -fno-omit-frame-pointer -fsanitize-address-use-after-scope -fsanitize=pointer-compare -fsanitize=pointer-subtract")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,undefined,leak")
    endif()
endif()

# Memory safety hardening
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORTIFY_SOURCE=3 -fstack-protector-strong -fstack-clash-protection -fcf-protection=full")

# PMR support
find_package(Threads REQUIRED)

# Optional: highway for portable SIMD
find_path(HIGHWAY_INCLUDE highway/hwy.h)
if(HIGHWAY_INCLUDE)
    message(STATUS "Found Highway SIMD at ${HIGHWAY_INCLUDE}")
    add_library(highway INTERFACE IMPORTED)
    set_target_properties(highway PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${HIGHWAY_INCLUDE}")
else()
    message(WARNING "Highway not found - SIMD kernels will use std::simd fallback")
endif()

# Optional: absl for flat_hash_map
find_package(absl CONFIG QUIET)
if(absl_FOUND)
    message(STATUS "Found Abseil")
else()
    message(WARNING "Abseil not found - using std::unordered_map fallback")
endif()

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(benchmarks)

# Install
install(DIRECTORY include/bwa DESTINATION include)
install(TARGETS bwa-core bwa-index bwa-align bwa-io bwa-utils
        EXPORT bwa-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)
install(EXPORT bwa-targets FILE bwa-config.cmake DESTINATION lib/cmake/bwa)