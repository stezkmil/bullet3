cmake_minimum_required(VERSION 3.24)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(BulletVbdGpu LANGUAGES CXX)
    option(INSTALL_LIBS "Install the GPU runtime backend" ON)
endif()

if(NOT WIN32)
    message(FATAL_ERROR "BUILD_VBD_GPU currently requires Windows (runtime loader and NVRTC compiler).")
endif()

if(VBD_PTX_COMPILER)
    # An explicit compiler (including Conan's pinned tool) must not silently fall back.
    find_package(CUDAToolkit 11.7 CONFIG REQUIRED)
    if(NOT EXISTS "${VBD_PTX_COMPILER}")
        message(FATAL_ERROR "VBD_PTX_COMPILER does not exist: ${VBD_PTX_COMPILER}")
    endif()
    set(vbd_ptx_command "${VBD_PTX_COMPILER}")
    set(vbd_ptx_dependency "${VBD_PTX_COMPILER}")
    message(STATUS "VBD GPU: using packaged PTX compiler ${VBD_PTX_COMPILER}")
else()
    # Standalone GitHub builds need only a local CUDA toolkit, not our Conan recipes.
    find_package(CUDAToolkit 11.7 MODULE REQUIRED)
    if(NOT TARGET CUDA::nvrtc)
        message(FATAL_ERROR "The CUDA toolkit must include NVRTC, or disable BUILD_VBD_GPU.")
    endif()
    if(CMAKE_CROSSCOMPILING)
        message(FATAL_ERROR "Standalone PTX compilation requires a native build or a supplied VBD_PTX_COMPILER.")
    endif()
    add_executable(BulletVbdCompilePtx compile_ptx.cpp)
    target_link_libraries(BulletVbdCompilePtx PRIVATE CUDA::nvrtc)
    target_compile_features(BulletVbdCompilePtx PRIVATE cxx_std_17)
    set_target_properties(BulletVbdCompilePtx PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tools/$<CONFIG>"
        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"
    )
    # CUDA 13 places DLLs in bin/x64; earlier toolkits use bin directly.
    file(GLOB vbd_nvrtc_dlls
        "${CUDAToolkit_BIN_DIR}/nvrtc64_*.dll" "${CUDAToolkit_BIN_DIR}/nvrtc-builtins64_*.dll"
        "${CUDAToolkit_BIN_DIR}/x64/nvrtc64_*.dll" "${CUDAToolkit_BIN_DIR}/x64/nvrtc-builtins64_*.dll")
    list(LENGTH vbd_nvrtc_dlls vbd_nvrtc_dll_count)
    if(vbd_nvrtc_dll_count LESS 2)
        message(FATAL_ERROR "NVRTC and its builtins DLL must be installed with the CUDA toolkit.")
    endif()
    add_custom_command(TARGET BulletVbdCompilePtx POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${vbd_nvrtc_dlls}
            "$<TARGET_FILE_DIR:BulletVbdCompilePtx>"
        VERBATIM)
    set(vbd_ptx_command "$<TARGET_FILE:BulletVbdCompilePtx>")
    set(vbd_ptx_dependency BulletVbdCompilePtx)
    message(STATUS "VBD GPU: building PTX compiler with CUDA toolkit ${CUDAToolkit_VERSION}")
endif()

set(vbd_ptx_header "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/guard_ptx.h")
add_custom_command(
    OUTPUT "${vbd_ptx_header}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>"
    COMMAND "${vbd_ptx_command}" "${CMAKE_CURRENT_SOURCE_DIR}/guard.cu" "${vbd_ptx_header}"
    DEPENDS ${vbd_ptx_dependency} guard.cu
    VERBATIM
)
add_library(BulletVbdGpu SHARED backend.cpp "${vbd_ptx_header}")
target_include_directories(BulletVbdGpu PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>" "${CMAKE_CURRENT_SOURCE_DIR}/../../src/BulletSoftBody")
target_link_libraries(BulletVbdGpu PRIVATE CUDA::cuda_driver)
target_compile_features(BulletVbdGpu PRIVATE cxx_std_17)

# The runtime loader uses a fixed DLL name and owns all allocations through the C API.
set_target_properties(BulletVbdGpu PROPERTIES
    DEBUG_POSTFIX "" RELEASE_POSTFIX "" RELWITHDEBINFO_POSTFIX "" MINSIZEREL_POSTFIX ""
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"
)
if(INSTALL_LIBS)
    install(TARGETS BulletVbdGpu RUNTIME DESTINATION bin COMPONENT VbdGpu)
endif()
