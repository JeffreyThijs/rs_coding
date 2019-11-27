#include <memory>
#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"
#include "rs_code.h"

class ReedSolomonEncoder {

    public:
        ReedSolomonEncoder(std::shared_ptr<rs_code_t> rs_code);
        ~ReedSolomonEncoder();
        galois::GaloisFieldPolynomial encode(std::string data);
        int get_data_block_size();
    private:
        galois::GaloisFieldPolynomial conv(galois::GaloisFieldPolynomial x, 
                                           const int nf,
                                           galois::GaloisFieldPolynomial y,
                                           const int ng);
        galois::GaloisFieldPolynomial get_generator();
        galois::GaloisFieldPolynomial string_to_poly(std::string data);
        galois::GaloisFieldPolynomial g_x;
        galois::GaloisField * gf;
        galois::GaloisFieldElement alpha;
        std::shared_ptr<rs_code_t> rs_code;
};