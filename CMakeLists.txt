cmake_minimum_required(VERSION 3.15)
project(availability_sim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -fsanitize=address")

# Find required packages
find_package(nlohmann_json 3.0 REQUIRED)
find_package(Threads REQUIRED)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/utils.cpp
    src/graph_structure.cpp
    src/disk.cpp
    src/erasure_coding.cpp
    src/simulation.cpp
    src/simulation_helpers.cpp
    src/simulation_core.cpp
    src/logger.cpp
    src/main.cpp
)

# Create executable
add_executable(availability_sim ${SOURCES})

# Link libraries
target_link_libraries(availability_sim
    PRIVATE
        nlohmann_json::nlohmann_json
        Threads::Threads
)

# Install target
install(TARGETS availability_sim DESTINATION bin)
