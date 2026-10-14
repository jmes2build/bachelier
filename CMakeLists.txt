cmake_minimum_required(VERSION 3.16)
project(bachelier LANGUAGES CXX VERSION 0.2.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Only the top-level build creates the CLI, tests and install rules. When this
# project is pulled in via add_subdirectory or FetchContent, a consumer wants
# the library and nothing else.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(BACHELIER_IS_TOP_LEVEL ON)
else()
    set(BACHELIER_IS_TOP_LEVEL OFF)
endif()

option(BACHELIER_BUILD_CLI "Build the bachelier command-line tool" ${BACHELIER_IS_TOP_LEVEL})
option(BACHELIER_BUILD_TESTS "Build the test suite" ${BACHELIER_IS_TOP_LEVEL})
option(BACHELIER_INSTALL "Generate install rules" ${BACHELIER_IS_TOP_LEVEL})

add_library(bachelier_core
    src/black_scholes.cpp
    src/implied_vol.cpp
    src/binomial.cpp
    src/normal_model.cpp
)
add_library(bachelier::bachelier ALIAS bachelier_core)

# The exported target takes its name from EXPORT_NAME, not from the target.
# Without this, installed consumers would have to link bachelier::bachelier_core
# while in-tree consumers used the bachelier::bachelier alias above - the same
# library under two names depending on how you got it.
set_target_properties(bachelier_core PROPERTIES EXPORT_NAME bachelier)

target_compile_features(bachelier_core PUBLIC cxx_std_20)

# Consumers of the installed package need the headers at their own prefix, so
# the include directory differs between building here and being installed.
target_include_directories(bachelier_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(MSVC)
    target_compile_options(bachelier_core PRIVATE /W4)
else()
    target_compile_options(bachelier_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(BACHELIER_BUILD_CLI)
    add_executable(bachelier src/main.cpp)
    target_link_libraries(bachelier PRIVATE bachelier_core)
endif()

if(BACHELIER_BUILD_TESTS)
    enable_testing()
    add_executable(bachelier_tests tests/tests.cpp)
    target_link_libraries(bachelier_tests PRIVATE bachelier_core)
    add_test(NAME bachelier_tests COMMAND bachelier_tests)
endif()

if(BACHELIER_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    install(TARGETS bachelier_core
        EXPORT bachelierTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
    if(BACHELIER_BUILD_CLI)
        install(TARGETS bachelier RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()

    install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(EXPORT bachelierTargets
        FILE bachelierTargets.cmake
        NAMESPACE bachelier::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bachelier
    )

    configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/bachelierConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/bachelierConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bachelier
    )
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/bachelierConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/bachelierConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/bachelierConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bachelier
    )
endif()
