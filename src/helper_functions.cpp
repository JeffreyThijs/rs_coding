#include "helper_functions.h"

galois::GaloisFieldPolynomial 
vector_product(galois::GaloisFieldPolynomial x,
               galois::GaloisFieldPolynomial y){

    auto elements = std::min(x.deg(), y.deg()) + 1;
    galois::GaloisFieldPolynomial z(x.field(), elements-1);
    for(unsigned int i=0; i < elements; i++){
        z[i] = x[i] * y[i];
    }

    return z;
}

galois::GaloisFieldElement
dot_product(galois::GaloisFieldPolynomial x,
               galois::GaloisFieldPolynomial y){

    auto elements = std::min(x.deg(), y.deg()) + 1;
    galois::GaloisFieldElement z(x.field(), 0);
    for(unsigned int i=0; i < elements; i++){
        z += (x[i] * y[i]);
    }
    return z;
}