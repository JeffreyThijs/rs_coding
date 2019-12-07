# Reed Solomon Coding

A basic software implementation of systematic Reed Solomon encoder and frequency domain decoder written in C++. Due to using the [Galois Field Arithmetic Library](http://www.partow.net/projects/galois/index.html) only Galois fields GF(q<sup>m</sup>) with the prime q = 2 and m is and integer power are supported. The decoder uses the Shuhong Gao algorithm which has advantages over classical decoding methods such as the Berlekamp-Massey algorithm, details on this can be found [here](http://www.math.clemson.edu/~sgao/papers/RS.pdf).

# Todo's

- Support for variable message size encoding/decoding
- Optimization of used functions: optimize cyclic IDFT, etc.
- Support for other BCH codes, RS codes are a special case of BCH codes where m<sub>0</sub> = 1.
- Improve rs_code config
- Add boost program_options for basic cli usage
- Add load/save to/from file

# Usage

```
$ ./rs --help
Options:
  -h [ --help ]                  Help screen
  -m [ --message ] arg           Input message to encode & decode
  -n [ --code-length ] arg (=36) Code length
  -t [ --errors ] arg (=6)       Error correction capability
  -q [ --power ] arg (=8)        Galois power GF(2^power)
  -i [ --input-file ] arg        Path to input file
  -o [ --output-file ] arg       Path to output file
  -d [ --decode ]                Add this flag to decode
```