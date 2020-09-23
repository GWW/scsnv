#cython: c_string_type=str, c_string_encoding=ascii, language_level=3str
from libc.stdint cimport int32_t, uint64_t, uint32_t
from libcpp cimport bool
from libcpp.string cimport string

cdef extern from './_align_mats.hpp':
    cdef cppclass CSRout:
        CSRout(Py_ssize_t N, Py_ssize_t M, Py_ssize_t D, int32_t * indptr, int32_t * indices, int32_t * ref, int32_t * alt) 

    Py_ssize_t sizeMatrices(size_t N, size_t M, int32_t * ref, int32_t * alt)
    Py_ssize_t mergeDups(int32_t * m, size_t N)
    void alignMatrices(CSRout * out, int32_t * ref, int32_t * alt)
