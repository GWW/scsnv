#cython: c_string_type=str, c_string_encoding=ascii, language_level=3str

import numpy as NP
cimport numpy as NP
import cython
cimport _amats as amats
from scipy.sparse import csr_matrix

def align_matrices(NP.ndarray[NP.int32_t, ndim=2] ref, NP.ndarray[NP.int32_t, ndim=2] alt):
    cdef Py_ssize_t N = ref.shape[0]
    cdef Py_ssize_t M = ref.shape[1]
    cdef Py_ssize_t D = amats.sizeMatrices(N, M, &ref[0, 0], &alt[0, 0])
    print(f'Nonzero {D} out of {N * M}')

    cdef NP.ndarray[NP.int32_t] indptr = NP.zeros(N + 1, dtype=NP.int32)
    cdef NP.ndarray[NP.int32_t] indices = NP.zeros(D, dtype=NP.int32)
    cdef NP.ndarray[NP.int32_t] refd = NP.zeros(D, dtype=NP.int32)
    cdef NP.ndarray[NP.int32_t] altd = NP.zeros(D, dtype=NP.int32)

    cdef amats.CSRout * ret = new amats.CSRout(N, M, D, &indptr[0], &indices[0], &refd[0], &altd[0])
    amats.alignMatrices(ret, &ref[0, 0], &alt[0, 0])

    return indptr, indices, refd, altd
