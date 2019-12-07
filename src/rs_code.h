#pragma once

#include <assert.h>

const unsigned int pp_power_16[17] = {1,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1};
const unsigned int pp_power_15[16] = {1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1};
const unsigned int pp_power_14[15] = {1,1,0,0,0,0,1,0,0,0,1,0,0,0,1};
const unsigned int pp_power_13[14] = {1,1,0,1,1,0,0,0,0,0,0,0,0,1};
const unsigned int pp_power_12[13] = {1,1,0,0,1,0,1,0,0,0,0,0,1};
const unsigned int pp_power_11[12] = {1,0,1,0,0,0,0,0,0,0,0,1};
const unsigned int pp_power_10[11] = {1,0,0,1,0,0,0,0,0,0,1};
const unsigned int pp_power_9[ 10] = {1,0,0,0,1,0,0,0,0,1};
const unsigned int pp_power_8[ 9 ] = {1,0,1,1,1,0,0,0,1};
const unsigned int pp_power_7[ 8 ] = {1,0,0,1,0,0,0,1};
const unsigned int pp_power_6[ 7 ] = {1,1,0,0,0,0,1};
const unsigned int pp_power_5[ 6 ] = {1,0,1,0,0,1};
const unsigned int pp_power_4[ 5 ] = {1,1,0,0,1};
const unsigned int pp_power_3[ 4 ] = {1,1,0,1};
const unsigned int pp_power_2[ 3 ] = {1,1,1};
const unsigned int pp_power_1[ 2 ] = {1,1};

struct rs_code_t {
    unsigned int n = 36; // degree: q - 1
    unsigned int t = 6;  // t-error correcting code
    unsigned int q = 2;
    unsigned int power = 8;
    unsigned int tt = 12;
    unsigned int k = 24; // for rs k = 2 * t
    unsigned int q_m = 255;
    unsigned int m0 = 1; // always for reed solomun
    const unsigned int * prime_poly = pp_power_8;
    std::shared_ptr<galois::GaloisField> gf;
    std::shared_ptr<galois::GaloisFieldElement> alpha;
    
    rs_code_t(unsigned int n=36, 
              unsigned int t=6, 
              unsigned int power=8, 
              unsigned int q=2)
    : n(n)
    , t(t)
    , q(q)
    {
        this->power = power;
        this->tt = 2 * t;
        this->k = n - tt;   
        this->q_m = q << (power - 1);
        this->prime_poly = get_default_prime_poly(power);
        this->gf = std::make_shared<galois::GaloisField>(this->power, 
                                                         this->prime_poly);
        this->alpha = std::make_shared<galois::GaloisFieldElement>(this->gf.get(), this->q);

        // some safety checks
        assert(tt > 0);
        assert(k > 0);
        assert(q_m > 0);
        assert(n > 0);
        assert(power > 0);
        assert(q > 0);
    }

    const unsigned int * get_default_prime_poly(int choice){
        switch(choice){
            case 16: return pp_power_16;
            case 15: return pp_power_15;
            case 14: return pp_power_14;
            case 13: return pp_power_13;
            case 12: return pp_power_12;
            case 11: return pp_power_11;
            case 10: return pp_power_10;
            case  9: return pp_power_9;
            case  8: return pp_power_8;
            case  7: return pp_power_7;
            case  6: return pp_power_6;
            case  5: return pp_power_5;
            case  4: return pp_power_4;
            case  3: return pp_power_3;
            case  2: return pp_power_2;
            case  1: return pp_power_1;
            default: return nullptr;
        }

    }

    galois::GaloisFieldPolynomial get_generator(){

        unsigned int size = this->tt;
        galois::GaloisField * gf_ptr = this->gf.get(); 
        galois::GaloisFieldElement alpha = *this->alpha;
        galois::GaloisFieldPolynomial g_x(gf_ptr, size);
        g_x[0] = 1;

        galois::GaloisFieldElement alpha_vector[2] = { galois::GaloisFieldElement(gf_ptr, 1),
                                                       galois::GaloisFieldElement(gf_ptr, 1) };

        for(unsigned int i=0; i < size; ++i){
            alpha_vector[0] = alpha ^ (i+1);
            galois::GaloisFieldPolynomial alpha_polynomial(gf_ptr, 1, alpha_vector);
            g_x = g_x * alpha_polynomial;
        }

        return g_x;
    }

};