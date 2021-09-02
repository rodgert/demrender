cmake_minimum_required(VERSION 3.17)
project(demrender)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  add_compile_options(-Wall -fPIC)
endif()

# set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if (NOT CMAKE_THREAD_LIBS_INIT)
  set(CMAKE_THREAD_LIBS_INIT "")
endif()

macro(add_catch_test name)
  if (TEST_REPORT_FORMAT)
    add_test(NAME $(name) COMMAND $(name) -r $(TEST_REPORT_FORMAT) -o "${name}.test_out.xml")
  else()
    add_test(NAME ${name} COMMAND ${name})
  endif()
endmacro()

include_directories(${PROJECT_SOURCE_DIR})

# add_subdirectory(doc)
add_subdirectory(src)

enable_testing()
# add_subdirectory(testsuite)

