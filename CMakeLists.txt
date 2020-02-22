cmake_minimum_required(VERSION 3.10)

project(splat VERSION 2.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


if(WIN32)
  set(SHAREDIR ".")
  set(DOCDIR "docs")
elseif(APPLE)
  set(SHAREDIR "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}")
  set(DOCDIR "share/docs/${PROJECT_NAME}")
  set(MANDIR "share/man")
elseif(UNIX) # Linux, BSD etc
  set(SHAREDIR "share/${PROJECT_NAME}")
  set(DOCDIR "share/docs/${PROJECT_NAME}")
  set(MANDIR "share/man")
endif()


add_subdirectory(src)
add_subdirectory(utils)
add_subdirectory(docs)

install(DIRECTORY sample_data DESTINATION "${SHAREDIR}")
