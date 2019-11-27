#include <memory>
#include "rs_encoder.h"


int main(){

    std::unique_ptr<ReedSolomonEncoder> encoder(new ReedSolomonEncoder());
    std::string data (encoder->get_data_block_size(), 0xa5);
    
    galois::GaloisFieldPolynomial enc_data = encoder->encode(data);

    std::cout << "enc_data = " << enc_data << std::endl;

    return 0;
}