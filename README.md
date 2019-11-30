# Reed Solomon Coding

A basic software implementation of systematic Reed Solomon encoder and frequency domain decoder written in C++. Due to using the [Galois Field Arithmetic Library](http://www.partow.net/projects/galois/index.html) only Galois fields GF(q<sup>m</sup>) with the prime q = 2 and m is and integer power are supported.

# Todo's

- Support for variable message size encoding/decoding
- Optimization of used functions: optimize cyclic IDFT, etc.
- Support for other BCH codes, RS codes are a special case of BCH codes where m<sub>0</sub> = 1.
- Improve rs_code config
- Add boost program_options for basic cli usage
- Add load/save to/from file