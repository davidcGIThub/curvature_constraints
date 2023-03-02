add_library(${PROJECT_NAME} SHARED 
    BsplineToBezier.hpp
    CubicEquationSolver.hpp
    SecondOrderCurvatureEvaluator.hpp
    ThirdOrderCurvatureEvaluator.hpp
    BsplineToBezier.cpp
    CubicEquationSolver.cpp
    SecondOrderCurvatureEvaluator.cpp
    ThirdOrderCurvatureEvaluator.cpp
    MDMAlgorithm.cpp
    MDMAlgorithm.hpp
    MDMAlgorithmClass.cpp
    MDMAlgorithmClass.hpp
    FifthOrderCurvatureEvaluator.cpp
    FifthOrderCurvatureEvaluator.hpp)

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")