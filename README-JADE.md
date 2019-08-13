# Installation instructions for the JADE cluster

* Create CMakeLists-\<yourname\>.txt from one of the existing files.

* Create a symbolic link to CMakeLists-\<yourname\>.txt.
    * `ln -s CMakeLists-<yourname>.txt CMakeLists.txt`

* Customise the CMakeLists-\<yourname\>.txt file accordingly.
  * Change the path to the compiler and the Eigen library.
  * Change the path to installation directory (*bin* in XCFD).
  * Add VTK_LEGACY_FORMAT preprocessor variable to the compiler options, if missing.

* Open a terminal and enter XCFD repository directory and enter the following commands
    * `module load gcc7/7.3.0`
    * `mkdir build`
    * `cd build`
    * `cmake ..`
    * `make install`

# Execution instructions

* Customise the `cfdjob.job` slurm script in the repository.