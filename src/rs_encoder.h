#pragma once

#include <memory>
#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"
#include "rs_code.h"

class ReedSolomonEncoder {

    public:
        ReedSolomonEncoder(std::shared_ptr<rs_code_t> rs_code);
        ~ReedSolomonEncoder();
        void encode(std::string data, std::string& enc_data);
        void encode(std::string data, galois::GaloisFieldPolynomial& enc_data);
    private:
        galois::GaloisFieldPolynomial _encode(std::string data);
        galois::GaloisFieldPolynomial string_to_poly(std::string data);
        std::string poly_to_string(galois::GaloisFieldPolynomial message);
        galois::GaloisFieldPolynomial g_x;
        std::shared_ptr<galois::GaloisField> gf;
        std::shared_ptr<galois::GaloisFieldElement> alpha;
        std::shared_ptr<rs_code_t> rs_code;
};