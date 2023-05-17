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

if(WIN32)
  set_target_properties(${PROJECT_NAME}
    PROPERTIES
      WIN32_EXECUTABLE ON
  )
elseif(APPLE)
  set_target_properties(${PROJECT_NAME}
    PROPERTIES
      MACOSX_BUNDLE ON
  )

  ## if(Qt6_FOUND)
  ##   get_target_property(__qmake_executable Qt6::qmake IMPORTED_LOCATION)
  ## else()
  ##   get_target_property(__qmake_executable Qt5::qmake IMPORTED_LOCATION)
  ## endif()
  ## get_filename_component(__qt_bin_dir ${__qmake_executable} DIRECTORY)
  ## find_program(MACDEPLOYQT_EXECUTABLE macdeployqt HINTS ${__qt_bin_dir})

  ## add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
  ##   COMMAND rm
  ##   ARGS -rf "${PROJECT_NAME}.dmg"
  ##   COMMENT "Remove previously generated disk image"
  ## )
  ## add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
  ##   COMMAND ${MACDEPLOYQT_EXECUTABLE}
  ##   ARGS "${PROJECT_NAME}.app" "-dmg"
  ##   COMMENT "Execute macdeployqt to create macOS bundle"
  ## )
endif()
