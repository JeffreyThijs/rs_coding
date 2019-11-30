#pragma once

#include <memory>
#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"
#include "rs_code.h"

class ReedSolomonDecoder {

    public:
        ReedSolomonDecoder(std::shared_ptr<rs_code_t> rs_code);
        ~ReedSolomonDecoder();
        void decode(std::string r_x, std::string& c_x);
        void decode(galois::GaloisFieldPolynomial r_x, std::string& c_x);
    private:
        std::string _decode(galois::GaloisFieldPolynomial r_x);
        std::string poly_to_string(galois::GaloisFieldPolynomial message);
        galois::GaloisFieldPolynomial string_to_poly(std::string message);
        galois::GaloisFieldPolynomial idft(galois::GaloisFieldPolynomial x);
        galois::GaloisFieldPolynomial calc_syndrome(galois::GaloisFieldPolynomial r_x);
        galois::GaloisFieldPolynomial calc_lambda_vector(galois::GaloisFieldPolynomial syndrome);
        galois::GaloisFieldPolynomial calc_error_vector(galois::GaloisFieldPolynomial syndrome,
                                                        galois::GaloisFieldPolynomial lambda_vector);
        std::shared_ptr<galois::GaloisField> gf;
        std::shared_ptr<galois::GaloisFieldElement> alpha;
        std::shared_ptr<rs_code_t> rs_code;
};