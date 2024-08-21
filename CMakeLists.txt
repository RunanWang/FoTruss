cmake_minimum_required(VERSION 3.22)
project(FoTruss)

set(CMAKE_CXX_STANDARD 17)

# Unify using clang as compiler
if (APPLE)
    SET(CMAKE_C_COMPILER /opt/homebrew/opt/llvm/bin/clang++)
    SET(CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++)
else ()
    # SET(CMAKE_C_COMPILER /usr/bin/clang++)
    # SET(CMAKE_CXX_COMPILER /usr/bin/clang++)
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -stdlib=libstdc++")
endif ()

# Debug/Release
#set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Find gflags
find_package(gflags REQUIRED)
set(FLAG_LIBRARIES gflags)

# Find glog
find_package(glog REQUIRED)
set(LOG_LIBRARIES glog::glog)

find_package(fmt REQUIRED)
set(FMT_LIBRARIES fmt::fmt)

add_executable(FoTruss
        main.cpp constant.h
        Utils/Timer.cpp Utils/Timer.h Utils/Memory.h Utils/logging.h
        DataStructure/Graph.h DataStructure/Graph.cpp DataStructure/MLGWithSchema.cpp DataStructure/MLGWithSchema.h
        DataStructure/Bucket.cpp DataStructure/Bucket.h DataStructure/Combo.cpp DataStructure/Combo.h
        DataStructure/EqualTree.cpp DataStructure/EqualTree.h DataStructure/UnionFind.cpp DataStructure/UnionFind.h
        Algorithm/Algo.cpp Algorithm/Algo.h Algorithm/TrussOpt.cpp Algorithm/TrussOpt.h
        test/testHeader.h test/SOTATest.cpp test/DecompositionTest.cpp test/SearchTimeTest.cpp test/IndexConstructionTest.cpp)

target_link_libraries(FoTruss
        ${FLAG_LIBRARIES}
        ${LOG_LIBRARIES}
        ${FMT_LIBRARIES}
        )