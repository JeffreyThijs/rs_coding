#pragma once

#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
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
        void encode_file(std::string src, std::string dst);
    private:
        unsigned int read_from_file(std::string src, std::vector<std::string>& blocks);
        void write_to_file(std::string dst, std::string blocks);
        galois::GaloisFieldPolynomial _encode(std::string data);
        galois::GaloisFieldPolynomial string_to_poly(std::string data);
        std::string poly_to_string(galois::GaloisFieldPolynomial message);
        galois::GaloisFieldPolynomial g_x;
        galois::GaloisField * gf;
        galois::GaloisFieldElement alpha;
        std::shared_ptr<rs_code_t> rs_code;
        
};