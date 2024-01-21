cmake_minimum_required (VERSION 3.1)
project (scsnv C CXX)

if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()

set (CMAKE_CXX_STANDARD 11)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
#set(HDF5_USE_STATIC_LIBRARIES ON)

SET( EXTERNAL_DIR "${PROJECT_SOURCE_DIR}/external" )
include(ExternalProject)

ExternalProject_Add(htslib
    DOWNLOAD_DIR ${EXTERNAL_DIR}
    DOWNLOAD_COMMAND curl -k -L https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 -o htslib-1.19.tar.bz2 &&
        rm -rf htslib && mkdir htslib && tar -xf htslib-1.19.tar.bz2 --strip-components=1 -C htslib && rm htslib-1.19.tar.bz2
    SOURCE_DIR ${EXTERNAL_DIR}/htslib
    CONFIGURE_COMMAND 
         ./configure --disable-libcurl  --without-libdeflate
            CC=${CMAKE_C_COMPILER}
            CFLAGS=-g\ -Wall\ -O2\ -fvisibility=hidden\ -fPIC ${ARMFLAGS}
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} CC=${CMAKE_C_COMPILER}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND cmake -E echo "Skipping install step."
)

SET( HTS_LIBRARY "${PROJECT_SOURCE_DIR}/external/htslib/libhts.a" )
SET( HTS_INCLUDE "${PROJECT_SOURCE_DIR}/external/htslib" )

ExternalProject_Add(bwa
    SOURCE_DIR ${EXTERNAL_DIR}/bwa
    GIT_REPOSITORY https://github.com/lh3/bwa.git#139f68fc4c3747813783a488aef2adc86626b01b
    CONFIGURE_COMMAND cd ${EXTERNAL_DIR}/bwa
        COMMAND cmake -E echo "BWA config"
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} CC=${CMAKE_C_COMPILER} CLAGS=${ARMFLAGS}
    INSTALL_COMMAND cmake -E echo "Skipping bwa install step."
    BUILD_IN_SOURCE 1
)

SET( BWA_LIBRARY "${PROJECT_SOURCE_DIR}/external/bwa/libbwa.a" )
SET( BWA_INCLUDE "${PROJECT_SOURCE_DIR}/external/bwa" )


ExternalProject_Add(argagg
    SOURCE_DIR ${EXTERNAL_DIR}/argagg
    GIT_REPOSITORY https://github.com/vietjtnguyen/argagg.git#6cd4a64ae09795c178d9ccba6d763c5fc9dcb72a
    CONFIGURE_COMMAND cmake -E echo "Skipping config step"
    BUILD_COMMAND cmake -E echo "Skipping build step."
    INSTALL_COMMAND cmake -E echo "Skipping install step."
    BUILD_IN_SOURCE 1
)

ExternalProject_Add(parallel_hashmap
    SOURCE_DIR ${EXTERNAL_DIR}/parallel_hashmap
    GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git#946ebad67a21212d11a0dd4deb7cdedc297d47bc
    CONFIGURE_COMMAND cmake -E echo "Skipping config step"
    BUILD_COMMAND cmake -E echo "Skipping build step."
    INSTALL_COMMAND cmake -E echo "Skipping install step."
    BUILD_IN_SOURCE 1
)
INCLUDE_DIRECTORIES("${EXTERNAL_DIR}/parallel_hashmap")

find_package( HDF5 REQUIRED COMPONENTS C CXX HL)
INCLUDE_DIRECTORIES(SYSTEM ${HDF5_INCLUDE_DIRS})
MESSAGE( "HDF5 LIBRARY: ${HDF5_LIBRARIES}" )
MESSAGE( "HDF5 LIBRARY CXX: ${HDF5_CXX_LIBRARIES}" )
MESSAGE( "HDF5 INCLUDE: ${HDF5_INCLUDE_DIRS}" )

MESSAGE( "BUILD TYPE: ${CMAKE_BUILD_TYPE}" )
MESSAGE( "INSTALL PREFIX: ${CMAKE_INSTALL_PREFIX}" )

SET( EXTERNAL_DIR "${PROJECT_SOURCE_DIR}/external" )

#add_subdirectory("bwa")
add_subdirectory("src")
