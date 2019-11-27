#include <memory>
#include "rs_encoder.h"


int main(){

    std::unique_ptr<ReedSolomonEncoder> encoder(new ReedSolomonEncoder());

    return 0;
}