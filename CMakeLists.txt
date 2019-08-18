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

set(CMAKE_CXX_COMPILER "g++")

set(COPT "-w -std=c++11 -g -O3 -Wall -Wextra -Wshadow -pedantic -Wno-uninitialized -Wno-sign-compare -Wno-write-strings  -Wno-deprecated  -Wno-unused  -Wno-format -Wno-narrowing  -Wno-reorder  -Wreturn-type -Wno-unused-parameter -fopenmp  -fpermissive -frounding-math -DVTK_LEGACY_FORMAT")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COPT}")


##################################
# include directories
##################################
include_directories(
src/base
/usr/include/eigen3 
/usr/include/eigen3/unsupported 
/usr/include/vtk-6.3 
)



##################################
# lib directories
##################################
link_directories(
/usr/lib/libblas
/usr/lib/lapack 
/usr/lib 
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

#target_link_libraries(infsup  stdc++ m superlu vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkRenderingFreeTypeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)

#target_link_libraries(incimplicit  stdc++ m vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkRenderingFreeTypeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3 superlu)

target_link_libraries(incexplicitSerial stdc++ m vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)

#target_link_libraries(incexplicitParallel  stdc++ m mpi mpi_cxx petsc_real metis vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)

#target_link_libraries(incexplicitdts  stdc++ m vtkCommonCore-6.3  vtkIOCore-6.3  vtkIOGeometry-6.3  vtkIOXML-6.3  vtkIOImage-6.3  vtkIOExport-6.3   vtkRenderingCore-6.3  vtkFiltersCore-6.3  vtkFiltersGeneral-6.3   vtkFiltersModeling-6.3  vtkFiltersSources-6.3  vtkRenderingAnnotation-6.3  vtkRenderingVolume-6.3  vtkRenderingVolumeOpenGL-6.3  vtkRenderingFreeTypeOpenGL-6.3  vtkInteractionStyle-6.3  vtkCommonExecutionModel-6.3  vtkCommonDataModel-6.3  vtkFiltersExtraction-6.3  vtkRenderingFreeType-6.3  vtkRenderingOpenGL-6.3  vtkIOParallel-6.3  vtkIOParallelXML-6.3)


#install(TARGETS infsup                  RUNTIME DESTINATION /home/chenna/Documents/myCode/XCFD/bin)
#install(TARGETS incimplicit             RUNTIME DESTINATION /home/chenna/Documents/myCode/XCFD/bin)
install(TARGETS incexplicitSerial       RUNTIME DESTINATION /home/chenna/Documents/myCode/XCFD/bin)
#install(TARGETS incexplicitParallel     RUNTIME DESTINATION /home/chenna/Documents/myCode/XCFD/bin)
#install(TARGETS incexplicitdts          RUNTIME DESTINATION /home/chenna/Documents/myCode/XCFD/bin)


