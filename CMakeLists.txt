cmake_minimum_required(VERSION 2.6)
project(igampap)

enable_language(CXX)

##################################
# project directories
##################################


#PROJECT_BINARY_DIR
#PROJECT_SOURCE_DIR

##################################
# compiler options
##################################


set(COPTGNUOMP "-w -std=c++11 -fopenmp -O3 -Wall -Wextra -Wshadow -pedantic -Wno-uninitialized -Wno-sign-compare -Wno-write-strings  -Wno-deprecated  -Wno-unused  -Wno-format -Wno-narrowing  -Wno-reorder  -Wreturn-type -Wno-unused-parameter -fpermissive -frounding-math -DVTK_LEGACY_FORMAT")
set(COPTGNU "-w -std=c++11 -O3 -Wall -Wextra -Wshadow -pedantic -Wno-uninitialized -Wno-sign-compare -Wno-write-strings  -Wno-deprecated  -Wno-unused  -Wno-format -Wno-narrowing  -Wno-reorder  -Wreturn-type -Wno-unused-parameter -fpermissive -frounding-math -DVTK_LEGACY_FORMAT")
set(COPTGDB "-w -std=c++11 -g -O0 -Wall -Wextra -Wshadow -pedantic -Wno-uninitialized -Wno-sign-compare -Wno-write-strings  -Wno-deprecated  -Wno-unused  -Wno-format -Wno-narrowing  -Wno-reorder  -Wreturn-type -Wno-unused-parameter -fpermissive -frounding-math -DVTK_LEGACY_FORMAT")

set(COPTPGI "-w -std=c++11 -mp -fast -DVTK_LEGACY_FORMAT")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COPTGNU}")


##################################
# include directories
##################################
include_directories(
src/base
)

##################################
# source files
##################################


file(GLOB sources1 src/base/*.h)
file(GLOB sources2 src/base/*.cpp)

##################################
##################################

#add_executable(infsup                 src/infsup.cpp ${sources1} ${sources2})
#add_executable(incimplicit             src/incimplicit.cpp ${sources1} ${sources2})
add_executable(incexplicitSerial       src/incexplicitSerial.cpp ${sources1} ${sources2})
#add_executable(incexplicitParallel     src/incexplicitParallel.cpp ${sources1} ${sources2})
#add_executable(incexplicitdts          src/incexplicitdts.cpp ${sources1} ${sources2})



target_link_libraries(incexplicitSerial stdc++ m )


#install(TARGETS infsup                  RUNTIME DESTINATION /home/chenna/Documents/myCode/XCFD/bin)
#install(TARGETS incimplicit             RUNTIME DESTINATION /home/chenna/Documents/myCode/XCFD/bin)
install(TARGETS incexplicitSerial       RUNTIME DESTINATION bin)
#install(TARGETS incexplicitParallel     RUNTIME DESTINATION /home/chenna/Documents/myCode/XCFD/bin)
#install(TARGETS incexplicitdts          RUNTIME DESTINATION /home/chenna/Documents/myCode/XCFD/bin)


