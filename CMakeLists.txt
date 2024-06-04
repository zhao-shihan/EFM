cmake_minimum_required(VERSION 3.14)

project(EFM VERSION 0.0.1 LANGUAGES CXX)

message(STATUS "Looking for ROOT")
find_package(ROOT 6.30.00 REQUIRED Tree RIO)
message(STATUS "Looking for ROOT - found (version: ${ROOT_VERSION})")

add_library(EFM INTERFACE)
add_library(EFM::EFM ALIAS EFM)
target_link_libraries(EFM INTERFACE ROOT::Tree ROOT::RIO)

include(CMakeDependentOption)

# https://cmake.org/cmake/help/latest/variable/PROJECT_IS_TOP_LEVEL.html
string(COMPARE EQUAL ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR} PROJECT_IS_TOP_LEVEL)

option(EFM_INSTALL "Generate and install EFM target" ${PROJECT_IS_TOP_LEVEL})
option(EFM_TEST "Build and perform EFM tests" ${PROJECT_IS_TOP_LEVEL})
cmake_dependent_option(EFM_ENABLE_UBSAN_IN_DEBUG_BUILD "Enable UndefinedBehaviorSanitizer in debug build" ON EFM_TEST OFF)

# The implementation generally assumes a platform that implements at least C++17 support
target_compile_features(EFM INTERFACE "cxx_std_17")

# Setup include directory
add_subdirectory(include)

if (EFM_TEST)
    enable_testing()
    add_subdirectory(test)
endif()

if (EFM_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/EFM" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    set(export_name "EFMConfig")
    set(namespace "EFM::")
    set(cmake_files_install_dir ${CMAKE_INSTALL_DATADIR}/cmake/EFM)

    install(TARGETS EFM EXPORT ${export_name} INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(EXPORT ${export_name} NAMESPACE ${namespace} DESTINATION ${cmake_files_install_dir})
    export(TARGETS EFM NAMESPACE ${namespace} FILE ${export_name}.cmake)

    set(EFM_config_version "${CMAKE_CURRENT_BINARY_DIR}/EFMConfigVersion.cmake")

    write_basic_package_version_file(${EFM_config_version} COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)

    install(FILES ${EFM_config_version} DESTINATION ${cmake_files_install_dir})
endif()
