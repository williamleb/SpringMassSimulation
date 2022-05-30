cmake_minimum_required(VERSION 3.15)

project(labo-3)

find_package(OpenMP REQUIRED)

set(HEADERS ParticleSimApplication.h ParticleSimGLCanvas.h ParticleSystem.h Solvers.hpp Vector2d.h )
set(SOURCES main.cpp ParticleSimApplication.cpp ParticleSimGLCanvas.cpp ParticleSystem.cpp )
add_executable(labo-3 ${SOURCES} ${HEADERS})

target_link_libraries(labo-3 labo-1 nanogui ${NANOGUI_EXTRA_LIBS} OpenMP::OpenMP_CXX)