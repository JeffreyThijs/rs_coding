#include "rs_encoder.h"

ReedSolomonEncoder::ReedSolomonEncoder(const int power,
                                       unsigned int * prime_polynomial)
: power(power)
{

    if(!prime_polynomial){
        prime_polynomial = default_prime_polynomial;
    }

    gf =  new galois::GaloisField(power, prime_polynomial);
    alpha = galois::GaloisFieldElement(gf, 2);
    rs_code = std::make_unique<rs_code_t>();
    g_x = get_generator();
}

ReedSolomonEncoder::~ReedSolomonEncoder(){
    delete gf;
}

galois::GaloisFieldPolynomial ReedSolomonEncoder::get_generator(){

    auto size = 2 * rs_code->t;
    galois::GaloisFieldPolynomial g_x(gf, size);
    g_x[0] = 1;

    galois::GaloisFieldElement alpha_vector[2] = { galois::GaloisFieldElement(gf, 1),
                                                   galois::GaloisFieldElement(gf, 1) };

    for(auto i=0; i < size; ++i){
        alpha_vector[0] = alpha ^ (i+1);
        galois::GaloisFieldPolynomial alpha_polynomial(gf, 1, alpha_vector);
        g_x = g_x * alpha_polynomial;
    }

    std::cout << "g(x)  = " << g_x              << std::endl;

    return g_x;
}