# XCFD
Explicit scheme for CFD using a new finite element scheme. At the moment limited to incompressible Navier-Stokes in 2D and 3D.

* Programming language: **C++**
* C++ standard: **C++11** (or above)
* Required third-party libraries:
  1.) [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  2.) [VTK](https://vtk.org/) (This is optional. Can output in VTK legacy format by defining the VTK_LEGACY_FORMAT variable in the compiler options.


## Compilation and building
* Create **build** and **bin** directories.
* Modify the CMakeLists.txt file accordingly.
  * Change the paths to the compilers, Eigen and VTK libraries.
  * Change the path to *bin*.
* Run the following commands (from the repository folder)
  * `cd build`
  * `cmake ..`
  * `make install`


## Execution
./incexplicitSerial \<inpfileprefix\>

One input file for the 2D lid-driven cavity benchmark is provided; this input is for the mesh B, Figure 13(b) in the paper. The control parameters in *control-parameters.dat* file are set to produce the results for Re=1000.

You can also test this by running the Python program *test-LDCT6-stru-meshB-Re1000.py*. This program makes a system call to the execution of the CFD program with the corresponding input file and compares the convergence data once the simulation is finished.

The simulation also produces 320 VTK files.
