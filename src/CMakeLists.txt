#SET(BWA_LIBRARY ${EXTERNAL_DIR}/bwa/libbwa.a)
#add_custom_command(OUTPUT ${BWA_LIBRARY}
#                   COMMAND make
#                   WORKING_DIRECTORY ${EXTERNAL_DIR}/bwa
#)

#add_custom_target(libbwa DEPENDS ${BWA_LIBRARY})

#MESSAGE( "BWA LIBRARY: ${BWA_LIBRARY}" )

SET(HTS_LIBRARY ${EXTERNAL_DIR}/htslib/libhts.a)
add_custom_command(OUTPUT ${HTS_LIBRARY}
    COMMAND autoreconf -i && ./configure --without-libdeflate CC=${CMAKE_C_COMPILER} && make
                   WORKING_DIRECTORY ${EXTERNAL_DIR}/htslib
)

add_custom_target(libhts DEPENDS ${HTS_LIBRARY})

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/bwa/)
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
    "psnvcounts.cpp"
    "paccuracy.cpp"
    "ptrim.cpp"

    "bwa/utils.c"
    "bwa/kthread.c"
    #    "bwa/kstring.c"
    "bwa/ksw.c"
    "bwa/bwt.c"
    "bwa/bntseq.c"
    "bwa/bwa.c"
    "bwa/bwamem.c"
    "bwa/bwamem_pair.c"
    "bwa/bwamem_extra.c"
    "bwa/malloc_wrap.c"
    "bwa/QSufSort.c"
    "bwa/bwt_gen.c"
    "bwa/rope.c"
    "bwa/rle.c"
    "bwa/is.c"
    "bwa/bwtindex.c"


)

target_link_libraries(scsnvlib ${HTS_LIBRARY} ${BWA_LIBRARY} ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES} ${HDF5_C_HL_LIBRARIES} z pthread lzma bz2 curl crypto )
add_dependencies(scsnvlib libhts)

add_executable(scsnv scsnv.cpp)
target_link_libraries(scsnv scsnvlib z)

add_executable(rfilter rfilter.cpp)
target_link_libraries(rfilter scsnvlib z)
