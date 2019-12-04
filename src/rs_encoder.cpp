#include "rs_encoder.h"

ReedSolomonEncoder::ReedSolomonEncoder(std::shared_ptr<rs_code_t> rs_code){

    this->rs_code = rs_code;
    this->gf = this->rs_code->gf.get();
    this->alpha = *this->rs_code->alpha.get();
    g_x = rs_code->get_generator();
}

ReedSolomonEncoder::~ReedSolomonEncoder(){
}

// systematic encoder
galois::GaloisFieldPolynomial ReedSolomonEncoder::_encode(std::string data){
     // convert to poly
    galois::GaloisFieldPolynomial enc_data = string_to_poly(data);

    // multiple with x^(n-k) = x^tt
    galois::GaloisFieldPolynomial x_tt(gf, rs_code->tt);
    x_tt[rs_code->tt] = 1;
    enc_data = enc_data * x_tt;

    // calculate the remainder of enc_data / g_x
    galois::GaloisFieldPolynomial r_x = enc_data % g_x;

    // the resulting code word is equal to enc_data + rx
    enc_data = enc_data + r_x;

    return enc_data;
}

void ReedSolomonEncoder::encode(std::string data, std::string& enc_data){
    enc_data = poly_to_string(_encode(data));
}

void ReedSolomonEncoder::encode(std::string data, galois::GaloisFieldPolynomial& enc_data){
    enc_data = _encode(data);
}

galois::GaloisFieldPolynomial ReedSolomonEncoder::string_to_poly(std::string data){
    
    galois::GaloisFieldPolynomial message(gf, rs_code->n - 1);

    for(int i = 0; i < rs_code->k; ++i){
        message[i] = data[rs_code->k - i - 1];
    }
    return message;
}

std::string ReedSolomonEncoder::poly_to_string(galois::GaloisFieldPolynomial message){
    
    std::string data(rs_code->n, 0x0);
    message.set_degree(rs_code->n + 1);

    for(int i = 0; i < rs_code->n; ++i){
        data[i] = message[i].poly();
    }

    return data;
}

int ReedSolomonEncoder::read_from_file(std::string src, std::vector<std::string>& blocks){

    int runs = 0;

    std::ifstream src_file(src, std::ifstream::binary);
    if (src_file) {
        src_file.seekg(0, src_file.end);
        int length = src_file.tellg();
        src_file.seekg(0, src_file.beg);

        runs = static_cast<int>(ceil(static_cast<float>(length) / static_cast<float>(rs_code->k)));
        std::string block(rs_code->k, 0x0);
    
        for(int i=0; i < runs; i++){
            src_file.read(&block.front(), rs_code->k);
            blocks.push_back(block);
            block.clear();
            block.resize(rs_code->k);
        }

        src_file.close();

    } else {
        std::cout << "Could not open " << src << std::endl;
    }

    return runs;
}

void ReedSolomonEncoder::write_to_file(std::string dst, std::string blocks){
    std::ofstream dst_file(dst, std::ifstream::binary);
    dst_file << blocks;
    dst_file.close();
}

void ReedSolomonEncoder::encode_file(std::string src, std::string dst){

    std::vector<std::string> src_blocks;
    int block_size = read_from_file(src, src_blocks);

    if(block_size == 0){
        std::cout << "Encoding failed!" << std::endl;
        return;
    }

    std::string block(rs_code->n, 0x0);
    std::string dst_blocks = "";
    
    for(int i=0; i < block_size; i++){
        encode(src_blocks[i], block);
        dst_blocks.append(block);
    }

    write_to_file(dst, dst_blocks);
}