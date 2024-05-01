#ifndef XMATH_H
#define XMATH_H
#include <stdint.h>

uint64_t xmath_residue_add(uint64_t x, uint64_t y);
uint64_t xmath_residue_sub(uint64_t x, uint64_t y);
uint64_t xmath_umul128(uint64_t multiplier, uint64_t multiplicand, uint64_t* product_hi);
uint64_t xmath_ui128_quotient_mod(uint64_t lo, uint64_t hi);
uint64_t xmath_residue_mul(uint64_t x, uint64_t y);
uint64_t xmath_residue_pow(uint64_t x, uint64_t y);
uint64_t xmath_inverse(uint64_t u, uint64_t v);
uint64_t xmath_residue_inv(uint64_t x);
uint64_t xmath_residue_sqrt(uint64_t what);

#endif
