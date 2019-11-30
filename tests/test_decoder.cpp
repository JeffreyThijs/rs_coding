#include <memory>
#include "rs_encoder.h"
#include "rs_decoder.h"

int main(){
    std::shared_ptr<rs_code_t> rs_code(new rs_code_t(15, 3, 4, 2));
    std::unique_ptr<ReedSolomonDecoder> decoder(new ReedSolomonDecoder(rs_code));
    galois::GaloisFieldElement x[11] = {
        galois::GaloisFieldElement(rs_code->gf.get(), 0),
        galois::GaloisFieldElement(rs_code->gf.get(), 0),
        galois::GaloisFieldElement(rs_code->gf.get(), 0),
        galois::GaloisFieldElement(rs_code->gf.get(), 1),
        galois::GaloisFieldElement(rs_code->gf.get(), 0),
        galois::GaloisFieldElement(rs_code->gf.get(), 0),
        galois::GaloisFieldElement(rs_code->gf.get(), 0),
        galois::GaloisFieldElement(rs_code->gf.get(), 0),
        galois::GaloisFieldElement(rs_code->gf.get(), 0),
        galois::GaloisFieldElement(rs_code->gf.get(), 0),
        galois::GaloisFieldElement(rs_code->gf.get(), 1)
    };

    galois::GaloisFieldPolynomial enc_data(rs_code->gf.get(), 10, x);
    std::cout << "enc_data: " << enc_data << std::endl;

    decoder->decode(enc_data);

    return 0;
}