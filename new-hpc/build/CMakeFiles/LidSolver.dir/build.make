# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/li/Desktop/hpc/new-hpc

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/li/Desktop/hpc/new-hpc/build

# Include any dependencies generated for this target.
include CMakeFiles/LidSolver.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/LidSolver.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LidSolver.dir/flags.make

CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o: CMakeFiles/LidSolver.dir/flags.make
CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o: ../LidDrivenCavity.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/li/Desktop/hpc/new-hpc/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o -c /home/li/Desktop/hpc/new-hpc/LidDrivenCavity.cpp

CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/li/Desktop/hpc/new-hpc/LidDrivenCavity.cpp > CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.i

CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/li/Desktop/hpc/new-hpc/LidDrivenCavity.cpp -o CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.s

CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o.requires:

.PHONY : CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o.requires

CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o.provides: CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o.requires
	$(MAKE) -f CMakeFiles/LidSolver.dir/build.make CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o.provides.build
.PHONY : CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o.provides

CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o.provides.build: CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o


CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o: CMakeFiles/LidSolver.dir/flags.make
CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o: ../LidDrivenCavitySolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/li/Desktop/hpc/new-hpc/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o -c /home/li/Desktop/hpc/new-hpc/LidDrivenCavitySolver.cpp

CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/li/Desktop/hpc/new-hpc/LidDrivenCavitySolver.cpp > CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.i

CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/li/Desktop/hpc/new-hpc/LidDrivenCavitySolver.cpp -o CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.s

CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o.requires:

.PHONY : CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o.requires

CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o.provides: CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o.requires
	$(MAKE) -f CMakeFiles/LidSolver.dir/build.make CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o.provides.build
.PHONY : CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o.provides

CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o.provides.build: CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o


CMakeFiles/LidSolver.dir/Poisson.cpp.o: CMakeFiles/LidSolver.dir/flags.make
CMakeFiles/LidSolver.dir/Poisson.cpp.o: ../Poisson.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/li/Desktop/hpc/new-hpc/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/LidSolver.dir/Poisson.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LidSolver.dir/Poisson.cpp.o -c /home/li/Desktop/hpc/new-hpc/Poisson.cpp

CMakeFiles/LidSolver.dir/Poisson.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LidSolver.dir/Poisson.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/li/Desktop/hpc/new-hpc/Poisson.cpp > CMakeFiles/LidSolver.dir/Poisson.cpp.i

CMakeFiles/LidSolver.dir/Poisson.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LidSolver.dir/Poisson.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/li/Desktop/hpc/new-hpc/Poisson.cpp -o CMakeFiles/LidSolver.dir/Poisson.cpp.s

CMakeFiles/LidSolver.dir/Poisson.cpp.o.requires:

.PHONY : CMakeFiles/LidSolver.dir/Poisson.cpp.o.requires

CMakeFiles/LidSolver.dir/Poisson.cpp.o.provides: CMakeFiles/LidSolver.dir/Poisson.cpp.o.requires
	$(MAKE) -f CMakeFiles/LidSolver.dir/build.make CMakeFiles/LidSolver.dir/Poisson.cpp.o.provides.build
.PHONY : CMakeFiles/LidSolver.dir/Poisson.cpp.o.provides

CMakeFiles/LidSolver.dir/Poisson.cpp.o.provides.build: CMakeFiles/LidSolver.dir/Poisson.cpp.o


# Object files for target LidSolver
LidSolver_OBJECTS = \
"CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o" \
"CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o" \
"CMakeFiles/LidSolver.dir/Poisson.cpp.o"

# External object files for target LidSolver
LidSolver_EXTERNAL_OBJECTS =

LidSolver: CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o
LidSolver: CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o
LidSolver: CMakeFiles/LidSolver.dir/Poisson.cpp.o
LidSolver: CMakeFiles/LidSolver.dir/build.make
LidSolver: /usr/lib/x86_64-linux-gnu/libboost_system.so
LidSolver: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
LidSolver: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
LidSolver: /usr/lib/x86_64-linux-gnu/libboost_system.so
LidSolver: /usr/lib/x86_64-linux-gnu/libboost_program_options.so
LidSolver: /usr/lib/x86_64-linux-gnu/libboost_regex.so
LidSolver: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
LidSolver: /usr/lib/libmpi.so
LidSolver: /usr/lib/x86_64-linux-gnu/libf77blas.so
LidSolver: /usr/lib/x86_64-linux-gnu/libatlas.so
LidSolver: /usr/lib/x86_64-linux-gnu/liblapack.so
LidSolver: /usr/lib/x86_64-linux-gnu/libf77blas.so
LidSolver: /usr/lib/x86_64-linux-gnu/libatlas.so
LidSolver: /usr/lib/x86_64-linux-gnu/liblapack.so
LidSolver: CMakeFiles/LidSolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/li/Desktop/hpc/new-hpc/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable LidSolver"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LidSolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LidSolver.dir/build: LidSolver

.PHONY : CMakeFiles/LidSolver.dir/build

CMakeFiles/LidSolver.dir/requires: CMakeFiles/LidSolver.dir/LidDrivenCavity.cpp.o.requires
CMakeFiles/LidSolver.dir/requires: CMakeFiles/LidSolver.dir/LidDrivenCavitySolver.cpp.o.requires
CMakeFiles/LidSolver.dir/requires: CMakeFiles/LidSolver.dir/Poisson.cpp.o.requires

.PHONY : CMakeFiles/LidSolver.dir/requires

CMakeFiles/LidSolver.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LidSolver.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LidSolver.dir/clean

CMakeFiles/LidSolver.dir/depend:
	cd /home/li/Desktop/hpc/new-hpc/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/li/Desktop/hpc/new-hpc /home/li/Desktop/hpc/new-hpc /home/li/Desktop/hpc/new-hpc/build /home/li/Desktop/hpc/new-hpc/build /home/li/Desktop/hpc/new-hpc/build/CMakeFiles/LidSolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LidSolver.dir/depend

