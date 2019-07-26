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

# C++ compiler

set(CMAKE_CXX_COMPILER "/usr/bin/g++")
#set(CMAKE_CXX_COMPILER "/opt/mylib/ompp-0.8.5/bin/kinst-ompp /usr/bin/g++")
#set(CXX "kinst-ompp g++")
#set(CMAKE_CXX_COMPILER "kinst-ompp g++")

#if(CXX_COMPILER_NAME MATCHES "icc*")
#    message ("Using INTEL C++ compiler")
#    set(COPT "-O3 -Wall -DSELF_CHECK -DUNIX")
#endif()

#if(CXX_COMPILER_NAME MATCHES "gcc*")
#    message ("Using GNU C++ compiler")
#    set(COPT "-std=c++11 -O3 -w -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=1  -Wno-write-strings -Wno-uninitialized -Wno-deprecated  -Wno-format -fpermissive -frounding-math")
#endif()

#message ("${CXX_COMPILER_NAME}")
#message ("${CXX}")

set(COPT "-w -std=c++11 -O0 -mavx2 -ftree-vectorize -ftree-vectorizer-verbose=1 -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -Wno-uninitialized -Wno-sign-compare -Wno-write-strings  -Wno-deprecated  -Wno-unused  -Wno-format -Wno-narrowing  -Wno-reorder  -Wreturn-type -Wno-unused-parameter -fopenmp  -fpermissive -frounding-math")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COPT}")


##################################
# include directories
##################################
include_directories(
src/base
/usr/include/c++/7 
/usr/include/eigen3 
/usr/include/eigen3/unsupported 
/usr/include/vtk-6.3 
/usr/include/openmpi 
/usr/include/petsc 
/usr/include/superlu 
/usr/include/boost 
/opt/mylib/ompp-0.8.5/include
)



##################################
# lib directories
##################################
link_directories(
/usr/lib/gcc/x86_64-linux-gnu/7
/usr/lib/x86_64-linux-gnu 
/usr/lib/libblas
/usr/lib/lapack 
/usr/lib 
/usr/lib/x86_64-linux-gnu  
/usr/lib/gcc/x86_64-linux-gnu/openmpi/lib 
/usr/lib/petsc/lib 
/opt/mylib/ompp-0.8.5/lib
)



##################################
# source files
##################################


file(GLOB sources1 src/base/*.h)
file(GLOB sources2 src/base/*.cpp)

##################################
##################################

#add_executable(mpap  src/program/main.cpp ${sources1} ${sources2} ${sources3})

set(VTKLIB62 "vtkCommonCore-6.3 vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkRenderingFreeTypeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3")


#target_link_libraries(mpap stdc++  gmp  mpfr  stdc++ gfortran gomp pthread m  blas  lapack  metis parmetis  petsc  mpi  mpi_cxx  superlu  pardiso500-GNU481-X86-64  boost_thread boost_system  CGAL CGAL_Core CGAL_ImageIO vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkRenderingFreeTypeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)

#target_link_libraries(igampap stdc++  gmp  mpfr  "${VTKLIB62}"  boost_thread boost_system  CGAL  CGAL_Core  CGAL_ImageIO m)

#add_executable(infsup                 src/infsup.cpp ${sources1} ${sources2})
#add_executable(incimplicit             src/incimplicit.cpp ${sources1} ${sources2})
add_executable(incexplicitSerial       src/incexplicitSerial.cpp ${sources1} ${sources2})
#add_executable(incexplicitParallel     src/incexplicitParallel.cpp ${sources1} ${sources2})
#add_executable(incexplicitdts          src/incexplicitdts.cpp ${sources1} ${sources2})

#target_link_libraries(infsup  stdc++ m superlu vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkRenderingFreeTypeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)

#target_link_libraries(incimplicit  stdc++ m vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkRenderingFreeTypeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3 superlu)

target_link_libraries(incexplicitSerial stdc++ m ompp mpi mpi_cxx petsc_real metis vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)

#target_link_libraries(incexplicitParallel  stdc++ m mpi mpi_cxx petsc_real metis vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)

#target_link_libraries(incexplicitdts  stdc++ m vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkRenderingFreeTypeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)


#install(TARGETS infsup                  RUNTIME DESTINATION /home/chenna/Documents/myCode/cfdexplicit/bin)
#install(TARGETS incimplicit             RUNTIME DESTINATION /home/chenna/Documents/myCode/cfdexplicit/bin)
install(TARGETS incexplicitSerial       RUNTIME DESTINATION /home/chenna/Documents/myCode/cfdexplicit/bin)
#install(TARGETS incexplicitParallel     RUNTIME DESTINATION /home/chenna/Documents/myCode/cfdexplicit/bin)
#install(TARGETS incexplicitdts          RUNTIME DESTINATION /home/chenna/Documents/myCode/cfdexplicit/bin)


