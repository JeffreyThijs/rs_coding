#pragma once

#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"

galois::GaloisFieldPolynomial 
vector_product(galois::GaloisFieldPolynomial x,
               galois::GaloisFieldPolynomial y);

galois::GaloisFieldElement
dot_product(galois::GaloisFieldPolynomial x,
               galois::GaloisFieldPolynomial y);