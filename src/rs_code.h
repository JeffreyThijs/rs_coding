const unsigned int pp_power_8[9] = {1,0,1,1,1,0,0,0,1};

struct rs_code_t {
    int q = 37;
    int n = 36; // degree: q - 1
    int t = 6;  // t-error correcting code
    int k = 24; // for rs k = 2 * t
    int tt = 12;
    int power = 8;
    const unsigned int * prime_poly = pp_power_8;
    
    rs_code_t(int n=36, int t=6, int power=8)
    : n(n)
    , t(t)
    {
        this->power = power;
        this->q = n + 1;
        this->tt = 2 * t;
        this->k = n - tt;   
        this->prime_poly = get_default_prime_poly(power);
    }

    const unsigned int * get_default_prime_poly(int choice){
        switch(choice){
            case 8:
                return pp_power_8;
            default:
                return nullptr;
        }
    }
};