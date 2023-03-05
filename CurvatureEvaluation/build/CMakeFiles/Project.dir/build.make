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
CMAKE_SOURCE_DIR = /home/david/Code/curvature_constraints/CurvatureEvaluation/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/david/Code/curvature_constraints/CurvatureEvaluation/build

# Include any dependencies generated for this target.
include CMakeFiles/Project.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Project.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Project.dir/flags.make

CMakeFiles/Project.dir/BsplineToBezier.o: CMakeFiles/Project.dir/flags.make
CMakeFiles/Project.dir/BsplineToBezier.o: /home/david/Code/curvature_constraints/CurvatureEvaluation/src/BsplineToBezier.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/CurvatureEvaluation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Project.dir/BsplineToBezier.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Project.dir/BsplineToBezier.o -c /home/david/Code/curvature_constraints/CurvatureEvaluation/src/BsplineToBezier.cpp

CMakeFiles/Project.dir/BsplineToBezier.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Project.dir/BsplineToBezier.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/CurvatureEvaluation/src/BsplineToBezier.cpp > CMakeFiles/Project.dir/BsplineToBezier.i

CMakeFiles/Project.dir/BsplineToBezier.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Project.dir/BsplineToBezier.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/CurvatureEvaluation/src/BsplineToBezier.cpp -o CMakeFiles/Project.dir/BsplineToBezier.s

CMakeFiles/Project.dir/CubicEquationSolver.o: CMakeFiles/Project.dir/flags.make
CMakeFiles/Project.dir/CubicEquationSolver.o: /home/david/Code/curvature_constraints/CurvatureEvaluation/src/CubicEquationSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/CurvatureEvaluation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Project.dir/CubicEquationSolver.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Project.dir/CubicEquationSolver.o -c /home/david/Code/curvature_constraints/CurvatureEvaluation/src/CubicEquationSolver.cpp

CMakeFiles/Project.dir/CubicEquationSolver.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Project.dir/CubicEquationSolver.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/CurvatureEvaluation/src/CubicEquationSolver.cpp > CMakeFiles/Project.dir/CubicEquationSolver.i

CMakeFiles/Project.dir/CubicEquationSolver.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Project.dir/CubicEquationSolver.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/CurvatureEvaluation/src/CubicEquationSolver.cpp -o CMakeFiles/Project.dir/CubicEquationSolver.s

CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.o: CMakeFiles/Project.dir/flags.make
CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.o: /home/david/Code/curvature_constraints/CurvatureEvaluation/src/SecondOrderCurvatureEvaluator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/CurvatureEvaluation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.o -c /home/david/Code/curvature_constraints/CurvatureEvaluation/src/SecondOrderCurvatureEvaluator.cpp

CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/CurvatureEvaluation/src/SecondOrderCurvatureEvaluator.cpp > CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.i

CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/CurvatureEvaluation/src/SecondOrderCurvatureEvaluator.cpp -o CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.s

CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.o: CMakeFiles/Project.dir/flags.make
CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.o: /home/david/Code/curvature_constraints/CurvatureEvaluation/src/ThirdOrderCurvatureEvaluator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/CurvatureEvaluation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.o -c /home/david/Code/curvature_constraints/CurvatureEvaluation/src/ThirdOrderCurvatureEvaluator.cpp

CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/CurvatureEvaluation/src/ThirdOrderCurvatureEvaluator.cpp > CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.i

CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/CurvatureEvaluation/src/ThirdOrderCurvatureEvaluator.cpp -o CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.s

CMakeFiles/Project.dir/MDMAlgorithm.o: CMakeFiles/Project.dir/flags.make
CMakeFiles/Project.dir/MDMAlgorithm.o: /home/david/Code/curvature_constraints/CurvatureEvaluation/src/MDMAlgorithm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/CurvatureEvaluation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/Project.dir/MDMAlgorithm.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Project.dir/MDMAlgorithm.o -c /home/david/Code/curvature_constraints/CurvatureEvaluation/src/MDMAlgorithm.cpp

CMakeFiles/Project.dir/MDMAlgorithm.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Project.dir/MDMAlgorithm.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/CurvatureEvaluation/src/MDMAlgorithm.cpp > CMakeFiles/Project.dir/MDMAlgorithm.i

CMakeFiles/Project.dir/MDMAlgorithm.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Project.dir/MDMAlgorithm.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/CurvatureEvaluation/src/MDMAlgorithm.cpp -o CMakeFiles/Project.dir/MDMAlgorithm.s

CMakeFiles/Project.dir/MDMAlgorithmClass.o: CMakeFiles/Project.dir/flags.make
CMakeFiles/Project.dir/MDMAlgorithmClass.o: /home/david/Code/curvature_constraints/CurvatureEvaluation/src/MDMAlgorithmClass.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/CurvatureEvaluation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/Project.dir/MDMAlgorithmClass.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Project.dir/MDMAlgorithmClass.o -c /home/david/Code/curvature_constraints/CurvatureEvaluation/src/MDMAlgorithmClass.cpp

CMakeFiles/Project.dir/MDMAlgorithmClass.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Project.dir/MDMAlgorithmClass.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/CurvatureEvaluation/src/MDMAlgorithmClass.cpp > CMakeFiles/Project.dir/MDMAlgorithmClass.i

CMakeFiles/Project.dir/MDMAlgorithmClass.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Project.dir/MDMAlgorithmClass.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/CurvatureEvaluation/src/MDMAlgorithmClass.cpp -o CMakeFiles/Project.dir/MDMAlgorithmClass.s

CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.o: CMakeFiles/Project.dir/flags.make
CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.o: /home/david/Code/curvature_constraints/CurvatureEvaluation/src/FifthOrderCurvatureEvaluator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/CurvatureEvaluation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.o -c /home/david/Code/curvature_constraints/CurvatureEvaluation/src/FifthOrderCurvatureEvaluator.cpp

CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/CurvatureEvaluation/src/FifthOrderCurvatureEvaluator.cpp > CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.i

CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/CurvatureEvaluation/src/FifthOrderCurvatureEvaluator.cpp -o CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.s

# Object files for target Project
Project_OBJECTS = \
"CMakeFiles/Project.dir/BsplineToBezier.o" \
"CMakeFiles/Project.dir/CubicEquationSolver.o" \
"CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.o" \
"CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.o" \
"CMakeFiles/Project.dir/MDMAlgorithm.o" \
"CMakeFiles/Project.dir/MDMAlgorithmClass.o" \
"CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.o"

# External object files for target Project
Project_EXTERNAL_OBJECTS =

libProject.so: CMakeFiles/Project.dir/BsplineToBezier.o
libProject.so: CMakeFiles/Project.dir/CubicEquationSolver.o
libProject.so: CMakeFiles/Project.dir/SecondOrderCurvatureEvaluator.o
libProject.so: CMakeFiles/Project.dir/ThirdOrderCurvatureEvaluator.o
libProject.so: CMakeFiles/Project.dir/MDMAlgorithm.o
libProject.so: CMakeFiles/Project.dir/MDMAlgorithmClass.o
libProject.so: CMakeFiles/Project.dir/FifthOrderCurvatureEvaluator.o
libProject.so: CMakeFiles/Project.dir/build.make
libProject.so: CMakeFiles/Project.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/david/Code/curvature_constraints/CurvatureEvaluation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX shared library libProject.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Project.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Project.dir/build: libProject.so

.PHONY : CMakeFiles/Project.dir/build

CMakeFiles/Project.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Project.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Project.dir/clean

CMakeFiles/Project.dir/depend:
	cd /home/david/Code/curvature_constraints/CurvatureEvaluation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/david/Code/curvature_constraints/CurvatureEvaluation/src /home/david/Code/curvature_constraints/CurvatureEvaluation/src /home/david/Code/curvature_constraints/CurvatureEvaluation/build /home/david/Code/curvature_constraints/CurvatureEvaluation/build /home/david/Code/curvature_constraints/CurvatureEvaluation/build/CMakeFiles/Project.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Project.dir/depend

