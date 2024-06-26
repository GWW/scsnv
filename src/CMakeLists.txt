#SET(BWA_LIBRARY ${EXTERNAL_DIR}/bwa/libbwa.a)
#add_custom_command(OUTPUT ${BWA_LIBRARY}
#                   COMMAND make
#                   WORKING_DIRECTORY ${EXTERNAL_DIR}/bwa
#)

#add_custom_target(libbwa DEPENDS ${BWA_LIBRARY})

#MESSAGE( "BWA LIBRARY: ${BWA_LIBRARY}" )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/)

add_library(scsnvlib
    "aux.cpp"
    "barcodes.cpp"
    "build.cpp"
    "index.cpp"
    "map_worker.cpp"
    "quant_worker.cpp"
    "reader.cpp"
    "transcript_align.cpp"
    "genome_align.cpp"

    "annotation.cpp"
    "iprocessor.cpp"
    "gtf.cpp"
    "gzstream.cpp"
    "sequence.cpp"
    "sbam_writer.cpp"

    "collapse_worker.cpp"
    "collapse_aux.cpp"
    "pileup.cpp"
    "pileup_worker.cpp"

    "pbase.cpp"
    "pbarcodes.cpp"
    "pmap.cpp"
    "pcollapse.cpp"
    "ppileup.cpp"
    "pmixture.cpp"
    "psnvcounts.cpp"
    "paccuracy.cpp"
    "ptrim.cpp"
)

target_link_libraries(scsnvlib ${HTS_LIBRARY} ${BWA_LIBRARY} ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES} ${HDF5_C_HL_LIBRARIES} z pthread lzma bz2 )

add_executable(scsnv scsnv.cpp)
target_link_libraries(scsnv scsnvlib z)

add_executable(rfilter rfilter.cpp)
target_link_libraries(rfilter scsnvlib z)
