
cmake_minimum_required(VERSION 3.21)
project(agentic_cfd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_KOKKOS "Enable Kokkos portability layer" OFF)
option(BUILD_TESTING "Build unit tests" ON)

if(USE_KOKKOS)
  find_package(Kokkos QUIET)
  if(NOT Kokkos_FOUND)
    include(FetchContent)
    FetchContent_Declare(kokkos
      GIT_REPOSITORY https://github.com/kokkos/kokkos.git
      GIT_TAG 4.2.01)
    FetchContent_MakeAvailable(kokkos)
  endif()
endif()

add_library(cfd_core
  src/solver.cpp
  src/io.cpp
)
if(USE_KOKKOS)
  target_link_libraries(cfd_core PUBLIC Kokkos::kokkos)
  target_compile_definitions(cfd_core PUBLIC AGENTICCFD_USE_KOKKOS)
endif()
target_include_directories(cfd_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(agentic_cfd src/main.cpp)
target_link_libraries(agentic_cfd PRIVATE cfd_core)

include(CTest)
if(BUILD_TESTING)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip)
  FetchContent_MakeAvailable(googletest)
  enable_testing()
  add_executable(solver_tests tests/test_solver.cpp)
  target_link_libraries(solver_tests PRIVATE cfd_core GTest::gtest_main)
  include(GoogleTest)
  gtest_discover_tests(solver_tests)
endif()
