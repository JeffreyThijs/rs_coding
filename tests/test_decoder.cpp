#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <memory>
#include "rs_encoder.h"
#include "rs_decoder.h"

BOOST_AUTO_TEST_CASE(TestDecoder_case_1){

    
    std::string decoded_message = "";
    std::shared_ptr<rs_code_t> rs_code(new rs_code_t(15, 3, 4, 2));
    std::unique_ptr<ReedSolomonDecoder> decoder(new ReedSolomonDecoder(rs_code));
    std::string original_message(rs_code->k, 0x0);

    galois::GaloisFieldElement x[2] = {
        galois::GaloisFieldElement(rs_code->gf.get(), 2),
        galois::GaloisFieldElement(rs_code->gf.get(), 11),
    };

    galois::GaloisFieldPolynomial enc_data(rs_code->gf.get(), 1, x);
    decoder->decode(enc_data, decoded_message);

    BOOST_CHECK_EQUAL(original_message, decoded_message);
}

BOOST_AUTO_TEST_CASE(TestDecoder_case_2){
    std::string decoded_message = "";
    std::shared_ptr<rs_code_t> rs_code(new rs_code_t(15, 3, 4, 2));
    std::unique_ptr<ReedSolomonDecoder> decoder(new ReedSolomonDecoder(rs_code));
    std::string original_message(rs_code->k, 0x0);

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

    galois::GaloisFieldPolynomial enc_data(rs_code->gf.get(), 1, x);
    decoder->decode(enc_data, decoded_message);

    BOOST_CHECK_EQUAL(original_message, decoded_message);
}