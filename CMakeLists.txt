PROJECT(obstacle LANGUAGES CUDA CXX)

CMAKE_MINIMUM_REQUIRED(VERSION 3.4)

include(CheckLanguage)
check_language(CUDA)


set (CMAKE_CXX_STANDARD 17)
# set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_52,code=sm_52;-gencode arch=compute_61,code=sm_61;-gencode arch=compute_35,code=sm_35;-std=c++11;-O3)
set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ;-I../../include --expt-extended-lambda;-Xcompiler -fPIC  -Xcompiler -fabi-version=2 -Xcompiler -fno-omit-frame-pointer ;--generate-code arch=compute_35,code=sm_35;-gencode arch=compute_52,code=sm_52;-std=c++11;-O3)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-sign-compare -Wno-int-in-bool-context -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# find_package(pthread)

#mysql
execute_process(COMMAND mysql_config --cflags
OUTPUT_VARIABLE MYSQL_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND mysql_config --libs
    OUTPUT_VARIABLE MYSQL_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)


#boost
find_package(Boost COMPONENTS program_options serialization wserialization iostreams REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

#cuda
set(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda)
SET(CUDA_SEPARABLE_COMPILATION ON)
find_package(CUDA REQUIRED)
SET(CUDA_SEPARABLE_COMPILATION ON)
include_directories(${CUDA_INCLUDE_DIRS})
link_directories(${CUDA_LIBRARY_DIRS})

#function_ref
include_directories(function_ref/include)

#hnswlib
include_directories(hnswlib)

#genie
include_directories(genie-dev-master/src)

#mongodb
find_package(libmongocxx REQUIRED)
find_package(libbsoncxx REQUIRED)
include_directories(${LIBMONGOCXX_INCLUDE_DIR} "/usr/local/include/mongocxx/v_noabi")
include_directories(${LIBBSONCXX_INCLUDE_DIR} "/usr/local/include/bsoncxx/v_noabi")
include_directories("json/single_include/nlohmann")

message("LIBBSONCXX_INCLUDE_DIR" ${LIBBSONCXX_INCLUDE_DIR})

AUX_SOURCE_DIRECTORY(. DIR_SRCS)
ADD_EXECUTABLE(interpolate_main "DataPreprocess/interpolate_main.cpp" "DataPreprocess/interpolation.cpp")
ADD_EXECUTABLE(heatmap_main "DataPreprocess/heatmap_main.cpp" "DataPreprocess/interpolation.cpp")
ADD_EXECUTABLE(index_test "index/index_test.cpp" "index/index.cpp")
ADD_EXECUTABLE(obst_test "main.cpp" "index/index.cpp" "obstacle/obstacle_detector.cpp")
ADD_EXECUTABLE(density_test "obstacle/obstacle_density_test.cpp" "obstacle/obstacle_detector.cpp" "lsh_density/index.cpp")
ADD_EXECUTABLE(kde_anomaly_test "obstacle/kde_anomaly_detection.cpp" "obstacle/obstacle_detector.cpp" "lsh_density/index.cpp")

target_link_libraries(interpolate_main  ${LIBMONGOCXX_LIBRARIES} ${LIBBSONCXX_LIBRARIES})
target_link_libraries(heatmap_main  ${LIBMONGOCXX_LIBRARIES} ${LIBBSONCXX_LIBRARIES})
target_link_libraries(index_test  ${LIBMONGOCXX_LIBRARIES} ${LIBBSONCXX_LIBRARIES})
target_link_libraries(density_test  ${LIBMONGOCXX_LIBRARIES} ${LIBBSONCXX_LIBRARIES}  ${Boost_LIBRARIES})
target_link_libraries(obst_test  ${LIBMONGOCXX_LIBRARIES} ${LIBBSONCXX_LIBRARIES}  ${Boost_LIBRARIES} )
target_link_libraries(kde_anomaly_test  ${LIBMONGOCXX_LIBRARIES} ${LIBBSONCXX_LIBRARIES}  ${Boost_LIBRARIES})



ADD_EXECUTABLE(density_test_gen_data "lsh_density/generate_datasets.cpp" "lsh_density/index.cpp")
target_link_libraries(density_test_gen_data  ${LIBMONGOCXX_LIBRARIES} ${LIBBSONCXX_LIBRARIES}  ${Boost_LIBRARIES})

ADD_EXECUTABLE(density_benchmarker "lsh_density/lsh_density_benchmarker.cpp" "lsh_density/index.cpp")
target_link_libraries(density_benchmarker  ${LIBMONGOCXX_LIBRARIES} ${LIBBSONCXX_LIBRARIES}  ${Boost_LIBRARIES})

# set(CUDA_ALL_LIBRARIES ${CUDA_LIBRARIES} ${CUDA_cudadevrt_LIBRARY} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cupti_LIBRARY} ${CUDA_curand_LIBRARY} ${CUDA_cusparse_LIBRARY} ${CUDA_npp_LIBRARY} ${CUDA_nppc_LIBRARY} ${CUDA_nppi_LIBRARY} ${CUDA_npps_LIBRARY} ${CUDA_nvcuvenc_LIBRARY} ${CUDA_nvcuvid_LIBRARY})
