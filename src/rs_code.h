#pragma once

const unsigned int pp_power_8[9] = {1,0,1,1,1,0,0,0,1};
const unsigned int pp_power_4[5] = {1,1,0,0,1};

struct rs_code_t {
    int m0 = 1; // always for reed solomun
    int q = 2;
    int n = 36; // degree: q - 1
    int t = 6;  // t-error correcting code
    int k = 24; // for rs k = 2 * t
    int tt = 12;
    int power = 8;
    int q_m = 255;
    const unsigned int * prime_poly = pp_power_8;
    std::shared_ptr<galois::GaloisField> gf;
    std::shared_ptr<galois::GaloisFieldElement> alpha;
    
    rs_code_t(int n=36, int t=6, int power=8, int q=2)
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
    }

    const unsigned int * get_default_prime_poly(int choice){
        switch(choice){
            case 8:
                return pp_power_8;
            case 4:
                return pp_power_4;
            default:
                return nullptr;
        }
    }

    galois::GaloisFieldPolynomial get_generator(){

        auto size = this->tt;
        galois::GaloisField * gf_ptr = this->gf.get(); 
        galois::GaloisFieldElement alpha = *this->alpha;
        galois::GaloisFieldPolynomial g_x(gf_ptr, size);
        g_x[0] = 1;

        galois::GaloisFieldElement alpha_vector[2] = { galois::GaloisFieldElement(gf_ptr, 1),
                                                       galois::GaloisFieldElement(gf_ptr, 1) };

        for(auto i=0; i < size; ++i){
            alpha_vector[0] = alpha ^ (i+1);
            galois::GaloisFieldPolynomial alpha_polynomial(gf_ptr, 1, alpha_vector);
            g_x = g_x * alpha_polynomial;
        }

        return g_x;
    }

};