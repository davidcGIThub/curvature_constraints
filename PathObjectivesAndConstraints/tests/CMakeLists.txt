set(TEST_NAME ${PROJECT_NAME}_UnitTest)

add_executable(${TEST_NAME}
    UnitTestBsplineToBezier.cpp
    UnitTestCubicEquationSolver.cpp
    UnitTestThirdOrderCurvatureEvaluator.cpp
    UnitTestMDMAlgorithm.cpp
    UnitTestMDMAlgorithmClass.cpp
    UnitTestDerivativeEvaluator.cpp
    UnitTestObjectiveFunctions.cpp
    UnitTestWaypointConstraints.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_SOURCE_DIR}/src")

target_link_libraries(${TEST_NAME} ${GTEST_LIBRARIES}
    pthread ${GTEST_MAIN_LIBRARIES} ${PROJECT_NAME}
    ${Eigen_LIBRARIES})