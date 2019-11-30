#include <memory>
#include "rs_encoder.h"
#include "rs_decoder.h"

int main(){

    std::shared_ptr<rs_code_t> rs_code(new rs_code_t());
    std::unique_ptr<ReedSolomonEncoder> encoder(new ReedSolomonEncoder(rs_code));
    std::unique_ptr<ReedSolomonDecoder> decoder(new ReedSolomonDecoder(rs_code));
    std::string data(rs_code->k, 0xa5);

    galois::GaloisFieldPolynomial enc_data = encoder->encode(data);
    std::cout << "enc_data = " << enc_data << std::endl;
    enc_data[0] = enc_data[0] + enc_data[1]; 

    decoder->decode(enc_data);

    
    return 0;
}