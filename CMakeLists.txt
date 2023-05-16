project(agglomerativeClustering LANGUAGES C CXX)
cmake_minimum_required(VERSION 3.12)

find_package(Qt6Widgets REQUIRED)
set(CMAKE_AUTOMOC ON)

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME}
  PRIVATE
    clustering.cpp
    main.cpp
)

target_link_libraries(${PROJECT_NAME}
  PRIVATE
    Qt6::Widgets
)
