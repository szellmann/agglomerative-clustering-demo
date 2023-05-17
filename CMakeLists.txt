project(agglomerativeClustering LANGUAGES C CXX)
cmake_minimum_required(VERSION 3.12)

find_package(Qt6 COMPONENTS Widgets)
if(NOT Qt6_FOUND)
  find_package(Qt5 COMPONENTS Widgets REQUIRED)
  set(__extra_libs ${__extra_libs} Qt5::Widgets)
else()
  set(__extra_libs ${__extra_libs} Qt6::Widgets)
endif()
set(CMAKE_AUTOMOC ON)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME}
  PRIVATE
    clustering.cpp
    main.cpp
)

target_link_libraries(${PROJECT_NAME}
  PRIVATE
    ${__extra_libs}
)
