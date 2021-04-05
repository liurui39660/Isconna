CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
SET(CMAKE_TOOLCHAIN_FILE $ENV{VCPKG_CMAKE})
PROJECT(Isconna)

SET(CMAKE_CXX_STANDARD 11)

INCLUDE_DIRECTORIES(
	src
	dep/AUROC
	$ENV{VCPKG_INCLUDE_DIR}
)
FIND_PACKAGE(mio REQUIRED)
FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(SQLite3 REQUIRED)

LINK_LIBRARIES(
	mio::mio mio::mio_base
	TBB::tbb
	SQLite::SQLite3
)

ADD_COMPILE_DEFINITIONS(
	_CRT_SECURE_NO_WARNINGS
	NOMINMAX
	SOLUTION_DIR="${CMAKE_SOURCE_DIR}/"
)

IF(DEFINED ENV{DATASET_DIR}) # Because I don't want a copy in all projects
	FILE(TO_CMAKE_PATH $ENV{DATASET_DIR} DATASET_DIR)
	ADD_COMPILE_DEFINITIONS(DATASET_DIR="${DATASET_DIR}/")
ELSE()
	ADD_COMPILE_DEFINITIONS(DATASET_DIR="${CMAKE_SOURCE_DIR}/data/")
ENDIF()

ADD_EXECUTABLE(Demo example/Demo.cpp)
ADD_EXECUTABLE(Experiment.Accuracy example/Experiment.Accuracy.cpp)
ADD_EXECUTABLE(Experiment.Efficiency example/Experiment.Efficiency.cpp)
ADD_EXECUTABLE(Experiment.Scalability.Edge example/Experiment.Scalability.Edge.cpp)
ADD_EXECUTABLE(Experiment.Scalability.Row example/Experiment.Scalability.Row.cpp)
ADD_EXECUTABLE(Experiment.Scalability.Column example/Experiment.Scalability.Column.cpp)
