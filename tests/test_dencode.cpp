#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <memory>
#include <iostream>
#include <fstream>
#include "rs_encoder.h"
#include "rs_decoder.h"

BOOST_AUTO_TEST_CASE(TestDEncoder){

    std::shared_ptr<rs_code_t> rs_code(new rs_code_t());
    std::unique_ptr<ReedSolomonEncoder> encoder(new ReedSolomonEncoder(rs_code));
    std::unique_ptr<ReedSolomonDecoder> decoder(new ReedSolomonDecoder(rs_code));
    std::string original_message = "abcdefghijklmnopqrstuvwx";
    std::string decoded_message = "";
    std::string encoded_message;

    for(int i=0; i < 100; i++){

        encoder->encode(original_message, encoded_message);

        // introduce random number of bit errors
        auto n_errors = rand() % rs_code->t;
        for(unsigned int i=0; i < n_errors; i++){
            encoded_message[rand() % rs_code->n] += 1;
        }

        decoder->decode(encoded_message, decoded_message);
        BOOST_CHECK_EQUAL(original_message, decoded_message);
    }
}

BOOST_AUTO_TEST_CASE(TestFileEncoder){

    std::shared_ptr<rs_code_t> rs_code(new rs_code_t());
    std::unique_ptr<ReedSolomonEncoder> encoder(new ReedSolomonEncoder(rs_code));
    std::unique_ptr<ReedSolomonDecoder> decoder(new ReedSolomonDecoder(rs_code));

    std::string input_file = "input.raw";
    std::string output_file = "output.raw";
    std::string encoded_file = "data.enc";
    std::string original_message = "";
    std::string received_message = "";
    int blocks = 100;

    // write to file
    std::ofstream x(input_file, std::ifstream::binary);
    for(auto i=0; i < blocks; i++){
        std::string block(rs_code->k, 0x32 + rand() % 26);
        x << block;
        original_message.append(block);
    }
    x.close();

    // file encoding/decoding
    encoder->encode_file(input_file, encoded_file);
    decoder->decode_file(encoded_file, output_file);

    // read from file
    std::ifstream y(output_file, std::ifstream::binary);
    y.seekg(0, y.end);
    int length = y.tellg();
    y.seekg(0, y.beg);

    // assert file length is a multiple of the block lengt to ensure correct decoding
    BOOST_CHECK_EQUAL(length % rs_code->k, 0);

    auto runs = length / rs_code->k;
    std::string block(rs_code->k, 0x0);

    for(unsigned int i=0; i < runs; i++){
        y.read(&block.front(), rs_code->k);
        received_message.append(block);
    }
    y.close();

    BOOST_CHECK_EQUAL(original_message, received_message);
}

BOOST_AUTO_TEST_CASE(test_non_div_block_size){

    std::shared_ptr<rs_code_t> rs_code(new rs_code_t());
    std::unique_ptr<ReedSolomonEncoder> encoder(new ReedSolomonEncoder(rs_code));
    std::unique_ptr<ReedSolomonDecoder> decoder(new ReedSolomonDecoder(rs_code));

    std::string input_file = "input.raw";
    std::string output_file = "output.raw";
    std::string encoded_file = "data.enc";
    std::string original_message = "";
    std::string received_message = "";
    int blocks = 100;

    // write to file
    std::ofstream x(input_file, std::ifstream::binary);
    for(auto i=0; i < blocks; i++){
        std::string block(rs_code->k, 0x32 + rand() % 26);
        x << block;
        original_message.append(block);
    }

    // add random block size at the end
    std::string _block(rand() % rs_code->k, 0x32 + rand() % 26);
    x << _block;
    original_message.append(_block);
    x.close();

    // file encoding/decoding
    encoder->encode_file(input_file, encoded_file);
    decoder->decode_file(encoded_file, output_file);

    // read from file
    std::ifstream y(output_file, std::ifstream::binary);
    y.seekg(0, y.end);
    int length = y.tellg();
    y.seekg(0, y.beg);

    std::string block(length, 0x0);
    y.read(&block.front(), length);
    received_message.append(block);
    y.close();

    BOOST_CHECK_EQUAL(original_message, received_message);
}