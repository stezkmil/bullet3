cmake_minimum_required(VERSION 3.24)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(BulletVbdGpu LANGUAGES CXX)
    option(INSTALL_LIBS "Install the GPU runtime backend" ON)
endif()

if(NOT WIN32)
    message(FATAL_ERROR "BUILD_VBD_GPU currently requires Windows (runtime loader and NVRTC compiler).")
endif()

# Both packages are provided by Conan; never select a machine-local CUDA toolkit.
find_package(CUDAToolkit 11.7 CONFIG REQUIRED)
find_package(CudaNvrtc CONFIG REQUIRED)
if(CMAKE_CROSSCOMPILING)
    message(FATAL_ERROR "The VBD PTX compiler currently requires a native Windows build.")
endif()
if(NOT IS_DIRECTORY "${VBD_NVRTC_BIN_DIR}")
    message(FATAL_ERROR "VBD_NVRTC_BIN_DIR must point to the Conan NVRTC package's bin directory.")
endif()

add_executable(BulletVbdCompilePtx compile_ptx.cpp)
target_link_libraries(BulletVbdCompilePtx PRIVATE CudaNvrtc::nvrtc)
target_compile_features(BulletVbdCompilePtx PRIVATE cxx_std_17)
set_target_properties(BulletVbdCompilePtx PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tools/$<CONFIG>"
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"
)
add_custom_command(TARGET BulletVbdCompilePtx POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${VBD_NVRTC_BIN_DIR}" "$<TARGET_FILE_DIR:BulletVbdCompilePtx>"
    VERBATIM)

set(vbd_ptx_header "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/guard_ptx.h")
add_custom_command(
    OUTPUT "${vbd_ptx_header}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>"
    COMMAND BulletVbdCompilePtx "${CMAKE_CURRENT_SOURCE_DIR}/guard.cu" "${vbd_ptx_header}"
    DEPENDS BulletVbdCompilePtx guard.cu
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
