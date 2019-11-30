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
        std::string decode(galois::GaloisFieldPolynomial enc_data);
    private:
        std::string poly_to_string(galois::GaloisFieldPolynomial message);
        std::shared_ptr<galois::GaloisField> gf;
        std::shared_ptr<galois::GaloisFieldElement> alpha;
        std::shared_ptr<rs_code_t> rs_code;
};