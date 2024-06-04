if(EFM_ENABLE_UBSAN_IN_DEBUG_BUILD)
    if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        add_compile_options($<$<CONFIG:Debug>:-fsanitize=undefined>
                            $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
                            $<$<CONFIG:RelWithDebInfo>:-fsanitize=undefined>
                            $<$<CONFIG:RelWithDebInfo>:-fno-omit-frame-pointer>)
        add_link_options($<$<CONFIG:Debug>:-fsanitize=undefined>
                         $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
                         $<$<CONFIG:RelWithDebInfo>:-fsanitize=undefined>
                         $<$<CONFIG:RelWithDebInfo>:-fno-omit-frame-pointer>)
        if(DEFINED CMAKE_CONFIGURATION_TYPES)
            # multi-config generator
            message(STATUS "UndefinedBehaviorSanitizer will be enabled in debug build (-fsanitize=undefined)")
        else()
            # single-config generator
            if(LOWERCASED_CMAKE_BUILD_TYPE STREQUAL "debug")
                message(STATUS "UndefinedBehaviorSanitizer enabled (-fsanitize=undefined)")
            endif()
        endif()
    endif()
endif()

# More warnings
if(CMAKE_COMPILER_IS_GNUCXX)
    list(APPEND EFM_COMPILE_OPTIONS -Wall -Wextra -Wduplicated-cond -Wnon-virtual-dtor -pedantic -Wundef -Wunused-macros)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    list(APPEND EFM_COMPILE_OPTIONS -WCL4 -Wmove -Wnon-virtual-dtor -pedantic -Wundef -Wunused-macros
                                     -Wno-gnu-zero-variadic-macro-arguments)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    list(APPEND EFM_COMPILE_OPTIONS /W4)
endif()

macro(add_executable_with_feature prefix feature)
    add_executable(${prefix}_${feature} ${prefix}.c++)
    target_link_libraries(${prefix}_${feature} PRIVATE EFM)
    target_compile_features(${prefix}_${feature} PRIVATE ${feature})
    target_compile_definitions(${prefix}_${feature} PRIVATE EFM_STATIC_TEST)
    target_compile_options(${prefix}_${feature} PRIVATE ${EFM_COMPILE_OPTIONS})
endmacro(add_executable_with_feature)

add_executable_with_feature(FieldMap3D cxx_std_17)

add_executable_with_feature(FieldMap3D cxx_std_20)

# add_executable_with_feature(FieldMap3D cxx_std_23)
