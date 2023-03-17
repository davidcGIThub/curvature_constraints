add_library(${PROJECT_NAME} SHARED 
    BsplineToBezier.cpp
    BsplineToBezier.hpp
    BsplineToMinvo.cpp
    BsplineToMinvo.hpp
    CubicEquationSolver.cpp
    CubicEquationSolver.hpp
    DerivativeEvaluator.cpp 
    DerivativeEvaluator.hpp
    ThirdOrderCurvatureEvaluator.cpp
    ThirdOrderCurvatureEvaluator.hpp
    MDMAlgorithmClass.cpp
    MDMAlgorithmClass.hpp
    WaypointConstraints.cpp
    WaypointConstraints.hpp
    ObjectiveFunctions.cpp
    ObjectiveFunctions.hpp
    CBindingHelper.cpp
    CBindingHelper.hpp
    RectPrismCollisionChecker.cpp
    RectPrismCollisionChecker.hpp
    ConvexHullCollisionChecker.cpp
    ConvexHullCollisionChecker.hpp)

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")