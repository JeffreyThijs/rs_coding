#include <boost/program_options.hpp>
#include <memory>
#include "rs_encoder.h"
#include "rs_decoder.h"

bool argparse(int argc, const char *argv[], boost::program_options::variables_map& vm){
   
    try{
        boost::program_options::options_description desc{"Options"};

        desc.add_options()
        ("help,h", "Help screen")
        ("message,m", boost::program_options::value<std::string>()->default_value(""), "Input message to encode & decode")
        ("code-length,n", boost::program_options::value<int>()->default_value(36), "Code length")
        ("errors,t", boost::program_options::value<int>()->default_value(6), "Error correction capability")
        ("power,q", boost::program_options::value<int>()->default_value(8), "Galois power GF(2^power)")
        ("input-file,i", boost::program_options::value<std::string>()->default_value(""), "Path to input file")
        ("output-file,o", boost::program_options::value<std::string>()->default_value(""), "Path to output file")
        ("decode,d", "Add this flag to decode")
        ;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        if (vm.count("help")){
            std::cout << desc << std::endl;
        }
        
    } catch (const boost::program_options::error &ex){
        std::cerr << ex.what() << '\n';
        return false;
    }

    return true;
}


int main(int argc, const char *argv[])
{
    boost::program_options::variables_map vm;
    if(!argparse(argc, argv, vm)){
        return -1;
    }

    std::shared_ptr<rs_code_t> rs_code(new rs_code_t(
        vm["code-length"].as<int>(),
        vm["errors"].as<int>(),
        vm["power"].as<int>()
    ));

    bool decode = vm.count("decode") ? true : false;
    std::unique_ptr<ReedSolomonEncoder> encoder(new ReedSolomonEncoder(rs_code));
    std::unique_ptr<ReedSolomonDecoder> decoder(new ReedSolomonDecoder(rs_code));

    // message encoding/decoding
    if(!vm["message"].as<std::string>().empty()){

        std::string original_message = vm["message"].as<std::string>();
        std::string decoded_message = "";
        std::string encoded_message = "";

        std::cout << "Original message: " << original_message << std::endl;

        encoder->encode(original_message, encoded_message);
        std::cout << "Encoded message:" << encoded_message << std::endl;

        decoder->decode(encoded_message, decoded_message);
        std::cout << "Decoded message: " << decoded_message << std::endl;

        return 0;
    }

    // file encoding/decoding   
    if(!vm["input-file"].as<std::string>().empty()){

        std::string input_file = vm["input-file"].as<std::string>();
        std::string output_file = (!vm["output-file"].as<std::string>().empty()) ? 
                                  vm["output-file"].as<std::string>() :
                                  (decode ? std::string("data.dec") : std::string("data.enc"));

        if(!decode){
            encoder->encode_file(input_file, output_file);
        } else {
            decoder->decode_file(input_file, output_file);    
        }
    }

    return 0;
}