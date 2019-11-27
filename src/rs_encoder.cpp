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

// systematic encoder
galois::GaloisFieldPolynomial ReedSolomonEncoder::encode(std::string data){

    // convert to poly
    galois::GaloisFieldPolynomial enc_data = string_to_poly(data);

    // multiple with x^(n-k) = x^tt
    galois::GaloisFieldPolynomial x_tt(gf, rs_code->tt);
    x_tt[rs_code->tt] = 1;
    enc_data = enc_data * x_tt;

    // calculate the remainder of enc_data / g_x
    galois::GaloisFieldPolynomial r_x = enc_data % g_x;

    // the resulting code word is equal to enc_data + rx
    enc_data = enc_data + r_x;

    return enc_data;
}

galois::GaloisFieldPolynomial ReedSolomonEncoder::get_generator(){

    auto size = rs_code->tt;
    galois::GaloisFieldPolynomial g_x(gf, size);
    g_x[0] = 1;

    galois::GaloisFieldElement alpha_vector[2] = { galois::GaloisFieldElement(gf, 1),
                                                   galois::GaloisFieldElement(gf, 1) };

    for(auto i=0; i < size; ++i){
        alpha_vector[0] = alpha ^ (i+1);
        galois::GaloisFieldPolynomial alpha_polynomial(gf, 1, alpha_vector);
        g_x = g_x * alpha_polynomial;
    }

    // std::cout << "g(x) = " << g_x << std::endl;

    return g_x;
}

galois::GaloisFieldPolynomial ReedSolomonEncoder::string_to_poly(std::string data){
    
    galois::GaloisFieldPolynomial message(gf, rs_code->n);

    for(int i = 0; i < rs_code->k; ++i){
        message[i] = data[rs_code->k - i - 1];
    }

    return message;
}

int ReedSolomonEncoder::get_data_block_size(){
    return rs_code->k;
}