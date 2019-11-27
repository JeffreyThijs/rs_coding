#include <memory>
#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"

// default prime_poly: D^8 + D^4 + D^3 + D^2 + 1

struct rs_code_t {
    int q = 37;
    int n = 36; // degree: q - 1
    int t = 6;  // t-error correcting code
    int k = 24; // for rs k = 2 * t

    rs_code_t(int n=36, int t=6)
    : n(n)
    , t(t)
    {
        this->q = n + 1;
        this->k = n - 2 * t; 
    }
};

class ReedSolomonEncoder {

    public:
        ReedSolomonEncoder(const int power=8, unsigned int * prime_polynomial=nullptr);
        ~ReedSolomonEncoder();
    private:
        galois::GaloisFieldPolynomial conv(galois::GaloisFieldPolynomial x, 
                                           const int nf,
                                           galois::GaloisFieldPolynomial y,
                                           const int ng);
        galois::GaloisFieldPolynomial get_generator();
        const int power;
        galois::GaloisFieldPolynomial g_x;
        galois::GaloisField * gf;
        galois::GaloisFieldElement alpha;
        std::unique_ptr<rs_code_t> rs_code;
        unsigned int default_prime_polynomial[9] = {1,0,1,1,1,0,0,0,1};
};