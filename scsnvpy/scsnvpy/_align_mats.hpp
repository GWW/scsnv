#pragma once
#include <stdint.h>
#include <iostream>
#include <fstream>

struct CSRout{
    CSRout(){

    }
    CSRout(size_t N, size_t M, size_t D, int32_t * indptr, int32_t * indices, int32_t * ref, int32_t * alt)
        : N(N), M(M), D(D), indptr(indptr), indices(indices), ref(ref), alt(alt)
    {

    }
    size_t N = 0;
    size_t M = 0;
    size_t D = 0; 
    int32_t * indptr = nullptr;
    int32_t * indices = nullptr;
    int32_t * ref = nullptr;
    int32_t * alt = nullptr;
};


inline size_t sizeMatrices(size_t N, size_t M, int32_t * ref, int32_t * alt){
    size_t D = 0;
    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j < M; j++){
            size_t idx = i * M + j;
            if(ref[idx] > 0 || alt[idx] > 0){
                D++;
            }
        }
    }
    return D;
}

inline void alignMatrices(CSRout * out, int32_t * ref, int32_t * alt){
    auto & o = *out;
    size_t D = 0;
    for(size_t i = 0; i < o.N; i++){
        o.indptr[i] = D;
        for(size_t j = 0; j < o.M; j++){
            size_t idx = i * o.M + j;
            if(ref[idx] > 0 || alt[idx] > 0){
                o.indices[D] = j;
                o.ref[D] = ref[idx];
                o.alt[D] = alt[idx];
                D++;
            }
        }
    }
    o.indptr[o.N] = D;

}
