#include "rs_encoder.h"

ReedSolomonEncoder::ReedSolomonEncoder(std::shared_ptr<rs_code_t> rs_code){

    this->rs_code = rs_code;
    this->gf = this->rs_code->gf;
    this->alpha = this->rs_code->alpha;
    g_x = rs_code->get_generator();
}

ReedSolomonEncoder::~ReedSolomonEncoder(){
}

// systematic encoder
galois::GaloisFieldPolynomial ReedSolomonEncoder::_encode(std::string data){
     // convert to poly
    galois::GaloisFieldPolynomial enc_data = string_to_poly(data);

    // multiple with x^(n-k) = x^tt
    galois::GaloisFieldPolynomial x_tt(gf.get(), rs_code->tt);
    x_tt[rs_code->tt] = 1;
    enc_data = enc_data * x_tt;

    // calculate the remainder of enc_data / g_x
    galois::GaloisFieldPolynomial r_x = enc_data % g_x;

    // the resulting code word is equal to enc_data + rx
    enc_data = enc_data + r_x;

    return enc_data;
}

void ReedSolomonEncoder::encode(std::string data, std::string& enc_data){
    enc_data = poly_to_string(_encode(data));
}

void ReedSolomonEncoder::encode(std::string data, galois::GaloisFieldPolynomial& enc_data){
    enc_data = _encode(data);
}

galois::GaloisFieldPolynomial ReedSolomonEncoder::string_to_poly(std::string data){
    
    galois::GaloisFieldPolynomial message(gf.get(), rs_code->n - 1);

    for(int i = 0; i < rs_code->k; ++i){
        message[i] = data[rs_code->k - i - 1];
    }
    return message;
}

std::string ReedSolomonEncoder::poly_to_string(galois::GaloisFieldPolynomial message){
    
    std::string data(rs_code->n, 0x0);
    message.set_degree(rs_code->n + 1);

    for(int i = 0; i < rs_code->n; ++i){
        data[i] = message[i].poly();
    }

    return data;
}