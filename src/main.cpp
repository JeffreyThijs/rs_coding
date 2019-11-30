#include <memory>
#include "rs_encoder.h"
#include "rs_decoder.h"

int main(){

    std::shared_ptr<rs_code_t> rs_code(new rs_code_t());
    std::unique_ptr<ReedSolomonEncoder> encoder(new ReedSolomonEncoder(rs_code));
    std::unique_ptr<ReedSolomonDecoder> decoder(new ReedSolomonDecoder(rs_code));
    std::string data = "abcdefghijklmnopqrstuvwx";
    std::string decoded_message = "";
    std::string enc_data;

    std::cout << "original message = " << data << std::endl;
    encoder->encode(data, enc_data);

    // introduce some errors
    enc_data[35] += 1;
    enc_data[33] += 1;
    enc_data[7] += 1;
    enc_data[2] += 1;
    std::cout << "enc_data = " << enc_data << std::endl;

    decoder->decode(enc_data, decoded_message);
    std::cout << "decoded_message: " << decoded_message << std::endl;

    
    return 0;
}