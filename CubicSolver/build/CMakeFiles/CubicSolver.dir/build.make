# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/david/Code/curvature_constraints/CubicSolver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/david/Code/curvature_constraints/CubicSolver/build

# Include any dependencies generated for this target.
include CMakeFiles/CubicSolver.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CubicSolver.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CubicSolver.dir/flags.make

CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.o: CMakeFiles/CubicSolver.dir/flags.make
CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.o: ../CubicEquationSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/CubicSolver/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.o -c /home/david/Code/curvature_constraints/CubicSolver/CubicEquationSolver.cpp

CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/CubicSolver/CubicEquationSolver.cpp > CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.i

CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/CubicSolver/CubicEquationSolver.cpp -o CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.s

# Object files for target CubicSolver
CubicSolver_OBJECTS = \
"CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.o"

# External object files for target CubicSolver
CubicSolver_EXTERNAL_OBJECTS =

libCubicSolver.a: CMakeFiles/CubicSolver.dir/CubicEquationSolver.cpp.o
libCubicSolver.a: CMakeFiles/CubicSolver.dir/build.make
libCubicSolver.a: CMakeFiles/CubicSolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/david/Code/curvature_constraints/CubicSolver/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libCubicSolver.a"
	$(CMAKE_COMMAND) -P CMakeFiles/CubicSolver.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CubicSolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CubicSolver.dir/build: libCubicSolver.a

.PHONY : CMakeFiles/CubicSolver.dir/build

CMakeFiles/CubicSolver.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CubicSolver.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CubicSolver.dir/clean

CMakeFiles/CubicSolver.dir/depend:
	cd /home/david/Code/curvature_constraints/CubicSolver/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/david/Code/curvature_constraints/CubicSolver /home/david/Code/curvature_constraints/CubicSolver /home/david/Code/curvature_constraints/CubicSolver/build /home/david/Code/curvature_constraints/CubicSolver/build /home/david/Code/curvature_constraints/CubicSolver/build/CMakeFiles/CubicSolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CubicSolver.dir/depend

