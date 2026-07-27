cmake_minimum_required(VERSION 3.23)

project(RayChasing LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

add_executable(RayChasing)

target_sources(RayChasing
  PRIVATE
    main.cc
)

add_executable(BraidPolynomialDemo
  braid_poly_demo.cc
)

target_include_directories(BraidPolynomialDemo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(BraidPolynomialTests
  braid_polynomial_tests.cc
)

target_include_directories(BraidPolynomialTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME braid_polynomial_demo COMMAND BraidPolynomialDemo 1,1,1)
add_test(NAME braid_polynomial_tests COMMAND BraidPolynomialTests)
