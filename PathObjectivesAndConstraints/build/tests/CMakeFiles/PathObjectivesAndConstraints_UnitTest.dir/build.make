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
CMAKE_SOURCE_DIR = /home/david/Code/curvature_constraints/PathObjectivesAndConstraints

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build

# Include any dependencies generated for this target.
include tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.o: ../tests/UnitTestBsplineToMinvo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestBsplineToMinvo.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestBsplineToMinvo.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestBsplineToMinvo.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.s

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.o: ../tests/UnitTestBsplineToBezier.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestBsplineToBezier.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestBsplineToBezier.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestBsplineToBezier.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.s

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.o: ../tests/UnitTestCubicEquationSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestCubicEquationSolver.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestCubicEquationSolver.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestCubicEquationSolver.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.s

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.o: ../tests/UnitTestThirdOrderCurvatureEvaluator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestThirdOrderCurvatureEvaluator.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestThirdOrderCurvatureEvaluator.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestThirdOrderCurvatureEvaluator.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.s

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.o: ../tests/UnitTestMDMAlgorithmClass.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestMDMAlgorithmClass.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestMDMAlgorithmClass.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestMDMAlgorithmClass.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.s

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.o: ../tests/UnitTestDerivativeEvaluator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestDerivativeEvaluator.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestDerivativeEvaluator.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestDerivativeEvaluator.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.s

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.o: ../tests/UnitTestObjectiveFunctions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestObjectiveFunctions.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestObjectiveFunctions.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestObjectiveFunctions.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.s

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.o: ../tests/UnitTestWaypointConstraints.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestWaypointConstraints.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestWaypointConstraints.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestWaypointConstraints.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.s

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.o: ../tests/UnitTestRectPrismCollisionChecker.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestRectPrismCollisionChecker.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestRectPrismCollisionChecker.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestRectPrismCollisionChecker.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.s

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.o: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/flags.make
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.o: ../tests/UnitTestConvexHullCollisionChecker.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.o"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.o -c /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestConvexHullCollisionChecker.cpp

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.i"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestConvexHullCollisionChecker.cpp > CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.i

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.s"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests/UnitTestConvexHullCollisionChecker.cpp -o CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.s

# Object files for target PathObjectivesAndConstraints_UnitTest
PathObjectivesAndConstraints_UnitTest_OBJECTS = \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.o" \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.o" \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.o" \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.o" \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.o" \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.o" \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.o" \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.o" \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.o" \
"CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.o"

# External object files for target PathObjectivesAndConstraints_UnitTest
PathObjectivesAndConstraints_UnitTest_EXTERNAL_OBJECTS =

tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToMinvo.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestBsplineToBezier.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestCubicEquationSolver.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestThirdOrderCurvatureEvaluator.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestMDMAlgorithmClass.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestDerivativeEvaluator.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestObjectiveFunctions.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestWaypointConstraints.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestRectPrismCollisionChecker.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/UnitTestConvexHullCollisionChecker.cpp.o
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/build.make
tests/PathObjectivesAndConstraints_UnitTest: /usr/lib/x86_64-linux-gnu/libgtest.a
tests/PathObjectivesAndConstraints_UnitTest: /usr/lib/x86_64-linux-gnu/libgtest_main.a
tests/PathObjectivesAndConstraints_UnitTest: src/libPathObjectivesAndConstraints.so
tests/PathObjectivesAndConstraints_UnitTest: tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking CXX executable PathObjectivesAndConstraints_UnitTest"
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/build: tests/PathObjectivesAndConstraints_UnitTest

.PHONY : tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/build

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/clean:
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/clean

tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/depend:
	cd /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/david/Code/curvature_constraints/PathObjectivesAndConstraints /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/tests /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests /home/david/Code/curvature_constraints/PathObjectivesAndConstraints/build/tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/PathObjectivesAndConstraints_UnitTest.dir/depend

