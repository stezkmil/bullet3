cmake_minimum_required(VERSION 3.24)
project(BulletVbdGpu LANGUAGES CXX)
find_package(CUDAToolkit REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/guard_ptx.h
 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compile_ptx.py --cuda ${CUDAToolkit_LIBRARY_ROOT} --source ${CMAKE_CURRENT_SOURCE_DIR}/guard.cu --output ${CMAKE_CURRENT_BINARY_DIR}/guard_ptx.h
 DEPENDS guard.cu compile_ptx.py VERBATIM)
add_library(BulletVbdGpu SHARED backend.cpp ${CMAKE_CURRENT_BINARY_DIR}/guard_ptx.h)
target_include_directories(BulletVbdGpu PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ../../src/BulletSoftBody)
target_link_libraries(BulletVbdGpu PRIVATE CUDA::cuda_driver)
target_compile_features(BulletVbdGpu PRIVATE cxx_std_17)

set_property(TARGET BulletVbdGpu PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
