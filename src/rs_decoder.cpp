#include "rs_decoder.h"

ReedSolomonDecoder::ReedSolomonDecoder(std::shared_ptr<rs_code_t> rs_code){
    this->rs_code = rs_code;
    this->gf = this->rs_code->gf.get();
    this->alpha = *this->rs_code->alpha.get();
}

ReedSolomonDecoder::~ReedSolomonDecoder(){
}

galois::GaloisFieldPolynomial ReedSolomonDecoder::idft(galois::GaloisFieldPolynomial x){

    galois::GaloisFieldPolynomial y(gf, x.deg());
    galois::GaloisFieldElement scale_factor(gf, ((rs_code->q_m - 1) % rs_code->q));

    // inverse fourier transform
    for(int i=0; i < rs_code->q_m - 1; i++){
        for(int j=0; j < rs_code->q_m - 1; j++){
            y[i] += x[j] * (alpha^(-j*i));
        }
        y[i] /= scale_factor;
    }

    return y;
}

galois::GaloisFieldPolynomial ReedSolomonDecoder::calc_syndrome(galois::GaloisFieldPolynomial r_x){

    galois::GaloisFieldPolynomial syndrome(gf, rs_code->tt-1);
    galois::GaloisFieldPolynomial helper_1(gf, r_x.deg());
    galois::GaloisFieldPolynomial helper_2;

    for(int j=0; j < (r_x.deg()+1); j++){
        helper_1[j] = (alpha)^j;
    }
    
    helper_2 = helper_1;

    for(int i=0; i < rs_code->tt; i++){
        syndrome[i] = dot_product(helper_2, r_x);
        helper_2 = vector_product(helper_2, helper_1);
    }

    syndrome.simplify();

    return syndrome;
}


galois::GaloisFieldPolynomial ReedSolomonDecoder::calc_lambda_vector(galois::GaloisFieldPolynomial syndrome){
    galois::GaloisFieldPolynomial lambda_prev(gf, 0);
    galois::GaloisFieldPolynomial lambda_current(gf, 0);
    galois::GaloisFieldPolynomial omega_prev(gf, rs_code->tt);
    galois::GaloisFieldPolynomial omega_current = syndrome;
    galois::GaloisFieldPolynomial q(gf, 0);
    galois::GaloisFieldPolynomial tmp(gf, 0);
    lambda_current[0] = 1;
    omega_prev[rs_code->tt] = 1;

    // std::cout << "lambda_prev: " << lambda_prev << std::endl;
    // std::cout << "lambda_current: " << lambda_current << std::endl;
    // std::cout << "omega_prev: " << omega_prev << std::endl;
    // std::cout << "omega_current: " << omega_current << std::endl;
    // std::cout << std::endl;

     while(omega_current.deg() > rs_code->t){
        // calc new omega
        q = omega_prev / omega_current;
        
        tmp = omega_prev - (q * omega_current);
        omega_prev = omega_current;
        omega_current = tmp;

        // calc new lambba
        tmp = lambda_prev - (q * lambda_current);
        lambda_prev = lambda_current;
        lambda_current = tmp;

        // std::cout << "q: " << q << std::endl;
        // std::cout << "lambda_prev: " << lambda_prev << std::endl;
        // std::cout << "lambda_current: " << lambda_current << std::endl;
        // std::cout << "omega_prev: " << omega_prev << std::endl;
        // std::cout << "omega_current: " << omega_current << std::endl;
        // std::cout << std::endl;
    }

    // final lambda_current needs to be used further
    // make last term 1 so E_vector[0] = 0 
    auto scale_factor = lambda_current[0].index();
    for(int i=0; i < lambda_current.deg()+1; i++){
        lambda_current[i] *= (alpha)^(-scale_factor);
    }

    return lambda_current;
}

galois::GaloisFieldPolynomial ReedSolomonDecoder::calc_error_vector(galois::GaloisFieldPolynomial syndrome,
                                                                    galois::GaloisFieldPolynomial lambda_vector){

    galois::GaloisFieldPolynomial E_vector(gf, rs_code->q_m-2);

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

    // only the values for E smaller than rs_code->m0 need to be found
    // for reed solomon m0 = 1 so only E[0] remains to be found
    galois::GaloisFieldPolynomial helper(gf, lambda_vector.deg()-1);
    auto scale_factor = lambda_vector[lambda_vector.deg()].index();
    for(int j=0; j < lambda_vector.deg(); j++){
        helper[j] = lambda_vector[j] / ((alpha)^scale_factor);
    }

    for(int i=0; i < rs_code->m0; i++){
        E_vector[i] = 0;
        for(int j=0; j < helper.deg()+1; j++){
            E_vector[i] += helper[j] * E_vector[lambda_vector.deg() - j];
        }

    }

    return E_vector;
}

// systematic frequency domain non-cyclic decoder using Shuhong Gao Algorithm:
// http://www.math.clemson.edu/~sgao/papers/RS.pdf
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
        // std::cout << "c_x: " << c_x << std::endl;
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
    
    galois::GaloisFieldPolynomial data(gf, rs_code->n);

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