add_library(scsnvlib
    "scsnv/aux.cpp"
    "scsnv/barcodes.cpp"
    "scsnv/build.cpp"
    "scsnv/index.cpp"
    "scsnv/map_worker.cpp"
    "scsnv/quant_worker.cpp"
    "scsnv/reader.cpp"
    "scsnv/tag_reader.cpp"
    "scsnv/transcript_align.cpp"
    "scsnv/genome_align.cpp"
    "scsnv/tag_reader.cpp"

    "util/annotation.cpp"
    "util/gtf.cpp"
    "util/gzstream.cpp"
    "util/sequence.cpp"
    "util/sbam_writer.cpp"

    "muts/collapse_worker.cpp"
    "muts/collapse_aux.cpp"
    "muts/pileup.cpp"
    "muts/pileup_worker.cpp"

    "progs/base.cpp"
    "progs/barcodes.cpp"
    "progs/map.cpp"
    "progs/quant.cpp"
    "progs/collapse.cpp"
    "progs/pileup.cpp"
)

SET(BWA_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/bwa/libbwa.a)
add_custom_command(OUTPUT ${BWA_LIBRARY}
                   COMMAND make
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bwa
)

add_custom_target(libbwa DEPENDS ${BWA_LIBRARY})

MESSAGE( "BWA LIBRARY: ${BWA_LIBRARY}" )


target_link_libraries(scsnvlib z ${HTS_LIBRARY} ${BWA_LIBRARY} ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES} ${HDF5_C_HL_LIBRARIES} pthread lzma bz2 curl )
add_dependencies(scsnvlib libbwa)

add_executable(scsnv scsnv.cpp)
target_link_libraries(scsnv scsnvlib)
