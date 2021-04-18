CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
SET(CMAKE_TOOLCHAIN_FILE $ENV{VCPKG_CMAKE})
PROJECT(Isconna)

SET(CMAKE_CXX_STANDARD 11)

INCLUDE_DIRECTORIES(
	include
	dep/AUROC
	dep/mio/include
	$ENV{VCPKG_INCLUDE_DIR}
)
FIND_PACKAGE(TBB)
FIND_PACKAGE(SQLite3)

IF(TBB_FOUND AND SQLite3_FOUND)
	LINK_LIBRARIES(
		TBB::tbb
		SQLite::SQLite3
	)
ENDIF()

ADD_COMPILE_DEFINITIONS(
	_CRT_SECURE_NO_WARNINGS
	NOMINMAX
	SOLUTION_DIR="${CMAKE_SOURCE_DIR}/"
)

# In case you put datasets elsewhere, e.g., using WSL
IF(NOT DEFINED ENV{DATASET_DIR})
	SET(ENV{DATASET_DIR} "${CMAKE_SOURCE_DIR}/data")
ENDIF()
FILE(TO_CMAKE_PATH $ENV{DATASET_DIR} DATASET_DIR)
ADD_COMPILE_DEFINITIONS(DATASET_DIR="${DATASET_DIR}/")

# Commonly used optimization flags
IF(NOT CMAKE_BUILD_TYPE STREQUAL Debug)
	IF(MSVC)
		ADD_COMPILE_OPTIONS(/arch:AVX2 /fp:fast /GL /Gw)
		ADD_LINK_OPTIONS(/LTCG /OPT:REF,ICF) # 4075: /LTCG overwrites /INCREMENTAL
	ELSE()
		ADD_COMPILE_OPTIONS(-Ofast -mavx2 -march=native -fwhole-program -funroll-loops)
	ENDIF()
ENDIF()

ADD_EXECUTABLE(Demo example/Demo.cpp)
ADD_EXECUTABLE(Experiment.Accuracy example/Experiment.Accuracy.cpp)
ADD_EXECUTABLE(Experiment.Efficiency example/Experiment.Efficiency.cpp)
ADD_EXECUTABLE(Experiment.Scalability.Edge example/Experiment.Scalability.Edge.cpp)
ADD_EXECUTABLE(Experiment.Scalability.Row example/Experiment.Scalability.Row.cpp)
ADD_EXECUTABLE(Experiment.Scalability.Column example/Experiment.Scalability.Column.cpp)
