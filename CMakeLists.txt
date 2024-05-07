cmake_minimum_required(VERSION 3.5)
project(P_and_N_Equations)

# Add executable
add_executable(main main.cpp)

#Find and add python libraries
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Add source files from src directory
file(GLOB SRC_FILES src/*.cpp)
target_sources(main PRIVATE ${SRC_FILES} src/global.cpp)

# Include headers from src directory
target_include_directories(main PRIVATE include)

#Include necessary python libraries
target_link_libraries(main PRIVATE ${PYTHON_LIBRARIES})

# Manually specify the path to the FFTW headers
include_directories("/usr/include")

# Manually specify the path to the FFTW library and link it
find_library(FFTW_LIBRARIES NAMES fftw3 PATHS "/usr/lib/x86_64-linux-gnu/")
target_link_libraries(main PRIVATE ${FFTW_LIBRARIES})
