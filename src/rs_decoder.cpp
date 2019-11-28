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

// systematic decoder
std::string ReedSolomonDecoder::decode(galois::GaloisFieldPolynomial enc_data){

    galois::GaloisFieldPolynomial lambda_prev(gf.get(), 0);
    galois::GaloisFieldPolynomial lambda_current(gf.get(), 0);
    galois::GaloisFieldPolynomial omega_prev(gf.get(), rs_code->tt);
    galois::GaloisFieldPolynomial omega_current(gf.get(), rs_code->tt-1);
    galois::GaloisFieldPolynomial g_x = rs_code->get_generator();
    galois::GaloisFieldPolynomial q(gf.get(), 0);
    galois::GaloisFieldPolynomial tmp(gf.get(), 0);
    galois::GaloisFieldPolynomial E_vector(gf.get(), rs_code->q_m-1);
    galois::GaloisFieldPolynomial e_vector(gf.get(), rs_code->q_m-1);
    galois::GaloisFieldPolynomial syndrome(gf.get(), rs_code->tt-1);
    lambda_current[0] = 1;
    omega_prev[rs_code->tt] = 1;


    galois::GaloisFieldPolynomial helper_1(gf.get(), enc_data.deg());
    for(int j=0; j < enc_data.deg()+1; j++){
        helper_1[j] = (*alpha)^j;
    }
    galois::GaloisFieldPolynomial helper_2 = helper_1;

    for(int i=0; i < rs_code->tt; i++){
        omega_current[i] = dot_product(helper_2, enc_data);
        helper_2 = vector_product(helper_2, helper_1);
        syndrome[i] = omega_current[i];
    }

    std::cout << "syndrome: " << omega_current << std::endl;

    while(omega_current.deg() > rs_code->t){

        // std::cout << "lambda_prev: " << lambda_prev << std::endl;
        // std::cout << "lambda_current: " << lambda_current << std::endl;
        // std::cout << "omega_prev: " << omega_prev << std::endl;
        // std::cout << "omega_current: " << omega_current << std::endl;

        // calc new omega
        q = omega_prev / omega_current;
        tmp = omega_prev - (q * omega_current);
        omega_prev = omega_current;
        omega_current = tmp;

        // calc new lambba
        tmp = lambda_prev - (q * lambda_current);
        lambda_current = lambda_prev;
        lambda_current = tmp;
    }

    // final lambda_current needs to be used further
    // make last term 1 so E_vector[0] = 0 
    lambda_current = lambda_current * (*alpha * (rs_code->q_m - 1 - lambda_current[0].index()));

    // fill known error vectors in error vector (= syndrome)
    for(int i=rs_code->m0; i < (rs_code->m0 + rs_code->tt); i++){
        E_vector[i] = syndrome[i - rs_code->m0];
    }

    // calculate the remaining the lambda equation
    std::cout << "lambda_current: " << lambda_current << std::endl;
    std::cout << "lambda_current: " << lambda_current.deg() << std::endl;
    for(int i=(rs_code->m0 + rs_code->tt); i < rs_code->q_m - 1; i++){
        // std::cout << "E[" << i << "] = ";
        for(int j=0; j < lambda_current.deg(); j++){
            E_vector[i] += lambda_current[j+1] * E_vector[i-j-1];
            //  std::cout << "L[" << j+1 << "] * E[" << i-j-1 << "] + ";
        }
        // std::cout << std::endl;
    }

    // E_vector.simplify();

    // inverse fourier transform
    for(int i=1; i < rs_code->q_m - 1; i++){
        for(int j=1; j < rs_code->q_m - 1; j++){
            e_vector[i] = E_vector[j] * (*alpha^(-j*i));
        }
    }
    std::cout << "e_vector: " << e_vector << std::endl;

    std::string a(5, 0x0);
    return a;
}

 std::string ReedSolomonDecoder::poly_to_string(galois::GaloisFieldPolynomial message){
    
    std::string data(rs_code->k, 0x0);

    for(int i = 0; i < rs_code->k; ++i){
        message[i] = data[rs_code->k - i - 1];
    }

    return data;
}