# Set minimum required CMake version
cmake_minimum_required(VERSION 3.10)

# Set project name
project(distribution_cpp)

# Enable C++17 standard (fallback to C++11 if necessary)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Find Python3 and include directories
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories(${Python3_INCLUDE_DIRS})

# Add the path to matplotlibcpp
include_directories(./matplotlibcpp)  # Replace with the actual path

# Add include directories for header files
include_directories(${CMAKE_SOURCE_DIR}/include)  # Add the include directory in the project root
include_directories(${CMAKE_SOURCE_DIR}/src)      # Add src directory if headers are also in src

# Optionally, link filesystem library for older GCC versions
if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lstdc++fs")
endif()

# Enable OpenMP for parallel computation
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Specify folder for compilation reports
set(REPORTS_DIR "${CMAKE_BINARY_DIR}/reports")
file(MAKE_DIRECTORY ${REPORTS_DIR})

# Enable optimizations and direct compilation reports to files
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fopt-info-vec-all=${REPORTS_DIR}/vectorization.txt -fopt-info-omp-all=${REPORTS_DIR}/openmp.txt")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /openmp /Qvec-report:2 /Qvec-report-file:${REPORTS_DIR}\\vectorization_msvc.txt")
endif()

# Specify your executable with correct source file path
add_executable(distribution_cpp src/main.cpp)  # Make sure this points to the correct path

# Link the Python3 library (for Python bindings) to your executable
target_link_libraries(distribution_cpp Python3::Python)

# Optionally, if using matplotlibcpp, link it to the executable
# target_link_libraries(distribution_cpp matplotlibcpp)

# Include additional libraries or dependencies here if necessary

