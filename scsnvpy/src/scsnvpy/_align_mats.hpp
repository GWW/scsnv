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

struct COOout{
    COOout(){

    }
    COOout(size_t N, size_t M, int32_t * xp, int32_t * yp, int32_t * ref, int32_t * alt)
        : N(N), M(M), xp(xp), yp(yp), ref(ref), alt(alt)
    {

    }
    size_t N = 0;
    size_t M = 0;
    size_t D = 0; 
    int32_t * xp = nullptr;
    int32_t * yp = nullptr;
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

inline size_t mergeDups(int32_t * m, size_t N){
    int32_t * out = m;
    for(int32_t * r = (m + 3); r != (m + N * 3); r+=3){
        //std::cout << "curr = " << (r - m) << " oidx = " << (out - m) << " curr = (" << r[0] << ", " << r[1] << ", " << r[2] << ") out = (" << out[0] << ", " << out[1] << ", " << out[2] << ") ";
        if(r[0] == out[0] && r[1] == out[1]){
            //std::cout << " dup\n";
            out[2] = r[2];
        }else{
            out+=3;
            //std::cout << " unique\n";
            out[0] = r[0];
            out[1] = r[1];
            out[2] = r[2];
        }
    }
    return ((out + 3) - m) / 3;
}
