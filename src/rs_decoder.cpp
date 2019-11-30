#include "rs_decoder.h"

ReedSolomonDecoder::ReedSolomonDecoder(std::shared_ptr<rs_code_t> rs_code){
    this->rs_code = rs_code;
    this->gf = this->rs_code->gf;
    this->alpha = this->rs_code->alpha;
}

ReedSolomonDecoder::~ReedSolomonDecoder(){
}

static galois::GaloisFieldPolynomial 
vector_product(galois::GaloisFieldPolynomial x,
               galois::GaloisFieldPolynomial y){

    auto elements = std::min(x.deg(), y.deg()) + 1;
    galois::GaloisFieldPolynomial z(x.field(), elements-1);
    for(auto i=0; i < elements; i++){
        z[i] = x[i] * y[i];
    }

    return z;
}

static galois::GaloisFieldElement
dot_product(galois::GaloisFieldPolynomial x,
               galois::GaloisFieldPolynomial y){

    auto elements = std::min(x.deg(), y.deg()) + 1;
    galois::GaloisFieldElement z(x.field(), 0);
    for(auto i=0; i < elements; i++){
        z += (x[i] * y[i]);
    }
    return z;
}

galois::GaloisFieldPolynomial ReedSolomonDecoder::idft(galois::GaloisFieldPolynomial x){

    galois::GaloisFieldPolynomial y(gf.get(), x.deg());
    galois::GaloisFieldElement scale_factor(gf.get(), ((rs_code->q_m - 1) % rs_code->q));

    // inverse fourier transform
    for(int i=0; i < rs_code->q_m - 1; i++){
        for(int j=0; j < rs_code->q_m - 1; j++){
            y[i] += x[j] * (*alpha^(-j*i));
        }
         y[i] /= scale_factor;
    }

    return y;
}

galois::GaloisFieldPolynomial ReedSolomonDecoder::calc_syndrome(galois::GaloisFieldPolynomial r_x){

    galois::GaloisFieldPolynomial syndrome(gf.get(), rs_code->tt-1);
    galois::GaloisFieldPolynomial helper_1(gf.get(), r_x.deg());
    galois::GaloisFieldPolynomial helper_2;

    for(int j=0; j < (r_x.deg()+1); j++){
        helper_1[j] = (*alpha)^j;
    }
    
    helper_2 = helper_1;

    for(int i=0; i < rs_code->tt; i++){
        syndrome[i] = dot_product(helper_2, r_x);
        helper_2 = vector_product(helper_2, helper_1);
    }

    // galois::GaloisFieldElement _alpha;
    // for(int i=0; i < rs_code->tt; i++){
    //     syndrome[i] = 0;
    //     for(int j=0; j < r_x.deg() + 1; j++){
    //         _alpha = ((*alpha)^(i+rs_code->m0)) ^ j;
    //         syndrome[i] += (r_x[j] * _alpha);
    //     }
    //     std::cout << std::endl;
    // }

    syndrome.simplify();

    return syndrome;
}


galois::GaloisFieldPolynomial ReedSolomonDecoder::calc_lambda_vector(galois::GaloisFieldPolynomial syndrome){
    galois::GaloisFieldPolynomial lambda_prev(gf.get(), 0);
    galois::GaloisFieldPolynomial lambda_current(gf.get(), 0);
    galois::GaloisFieldPolynomial omega_prev(gf.get(), rs_code->tt);
    galois::GaloisFieldPolynomial omega_current = syndrome;
    galois::GaloisFieldPolynomial q(gf.get(), 0);
    galois::GaloisFieldPolynomial tmp(gf.get(), 0);
    lambda_current[0] = 1;
    omega_prev[rs_code->tt] = 1;

     while(omega_current.deg() > rs_code->t){

        // std::cout << "lambda_prev: " << lambda_prev << std::endl;
        // std::cout << "lambda_current: " << lambda_current << std::endl;
        // std::cout << "omega_prev: " << omega_prev << std::endl;
        // std::cout << "omega_current: " << omega_current << std::endl;
        // std::cout << std::endl;

        // calc new omega
        q = omega_prev / omega_current;
        tmp = omega_prev - (q * omega_current);
        omega_prev = omega_current;
        omega_current = tmp;

        // calc new lambba
        tmp = lambda_prev - (q * lambda_current);
        lambda_prev = lambda_current;
        lambda_current = tmp;
    }

    // final lambda_current needs to be used further
    // make last term 1 so E_vector[0] = 0 
    lambda_current = lambda_current * (*alpha * (rs_code->q_m - 1 - lambda_current[0].index()));

    return lambda_current;
}

galois::GaloisFieldPolynomial ReedSolomonDecoder::calc_error_vector(galois::GaloisFieldPolynomial syndrome,
                                                                    galois::GaloisFieldPolynomial lambda_vector){

    galois::GaloisFieldPolynomial E_vector(gf.get(), rs_code->q_m-2);

    // fill known error vectors in error vector (= syndrome)
    for(int i=rs_code->m0; i < (rs_code->m0 + rs_code->tt); i++){
        E_vector[i] = syndrome[i - rs_code->m0];
    }

    // calculate the remaining the lambda equation
    for(int i=(rs_code->m0 + rs_code->tt); i < rs_code->q_m - 1; i++){
        for(int j=0; j < lambda_vector.deg(); j++){
            E_vector[i] += lambda_vector[j+1] * E_vector[i-j-1];
        }
    }

    return E_vector;
}

// systematic decoder
std::string ReedSolomonDecoder::_decode(galois::GaloisFieldPolynomial r_x){

    galois::GaloisFieldPolynomial c_x;
    galois::GaloisFieldPolynomial syndrome;
    galois::GaloisFieldPolynomial lambda_vector;
    galois::GaloisFieldPolynomial E_vector;
    galois::GaloisFieldPolynomial e_vector;

    // copy r_x into c_x
    c_x = r_x;
    
    // calc syndrome
    syndrome = calc_syndrome(r_x);
    // std::cout << "syndrome: " << syndrome << std::endl;

    if(syndrome.valid()){
        // calc lambda vector
        lambda_vector = calc_lambda_vector(syndrome);
        // std::cout << "lambda_vector: " << lambda_vector << std::endl;

        // calc error vector (in fourier domain)
        E_vector = calc_error_vector(syndrome, lambda_vector);
        // std::cout << "E_vector: " << E_vector << std::endl;

        // transform the error vector to the time domain
        e_vector = idft(E_vector);
        // std::cout << "e_vector: " << e_vector << std::endl;

        // correct r_x with e_vector
        c_x = r_x + e_vector;
        c_x.set_degree(rs_code->n + 1);
        std::cout << "c_x: " << c_x << std::endl;
    } else {
        // std::cout << "not valid!" << std::endl;
    }

    return poly_to_string(c_x);
}

void ReedSolomonDecoder::decode(galois::GaloisFieldPolynomial r_x, std::string& c_x){
     c_x = _decode(r_x);
}

void ReedSolomonDecoder::decode(std::string r_x, std::string& c_x){
    c_x = _decode(string_to_poly(r_x));
}

galois::GaloisFieldPolynomial ReedSolomonDecoder::string_to_poly(std::string message){
    
    galois::GaloisFieldPolynomial data(gf.get(), rs_code->n);

    for(int i = 0; i < rs_code->n; ++i){
        data[i] = message[i];
    }

    return data;
}

std::string ReedSolomonDecoder::poly_to_string(galois::GaloisFieldPolynomial message){
    
    std::string data(rs_code->k, 0x0);
    message.set_degree(rs_code->n + 1);

    for(int i = 0; i < rs_code->k; ++i){
        data[i] = message[rs_code->k - i - 1 + rs_code->tt].poly();
    }

    return data;
}