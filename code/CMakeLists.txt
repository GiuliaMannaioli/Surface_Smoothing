cmake_minimum_required(VERSION 3.0)

SET(CMAKE_EXPORT_COMPILE_COMMANDS 1)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

project(IGRproject)

#add_definitions(-DSUPPORT_OPENGL_45)

add_executable(
  ${PROJECT_NAME}
  main.cpp
  #Error.cpp # Only if your system supports OpenGL 4.3 or later; don't forget to replace glad.
  Mesh.cpp
  ShaderProgram.cpp)

add_subdirectory(dep/glad)
target_link_libraries(${PROJECT_NAME} PRIVATE glad)

add_subdirectory(dep/glfw)
target_link_libraries(${PROJECT_NAME} PRIVATE glfw)

add_subdirectory(dep/glm)
target_link_libraries(${PROJECT_NAME} PRIVATE glm)

target_link_libraries(${PROJECT_NAME} PRIVATE ${CMAKE_DL_LIBS})

add_custom_command(TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_CURRENT_SOURCE_DIR})
