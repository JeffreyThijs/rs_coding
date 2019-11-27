#include "rs_encoder.h"

ReedSolomonEncoder::ReedSolomonEncoder(std::shared_ptr<rs_code_t> rs_code){

    this->rs_code = rs_code;

    gf =  new galois::GaloisField(rs_code->power, 
                                  rs_code->prime_poly);
    alpha = galois::GaloisFieldElement(gf, 2);
    g_x = rs_code->get_generator(gf);
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
galois::GaloisFieldPolynomial ReedSolomonEncoder::string_to_poly(std::string data){
    
    galois::GaloisFieldPolynomial message(gf, rs_code->n);

    for(int i = 0; i < rs_code->k; ++i){
        message[i] = data[rs_code->k - i - 1];
    }

    return message;
}