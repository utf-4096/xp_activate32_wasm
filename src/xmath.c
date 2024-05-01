#include "defs.h"
#include <stdint.h>

static inline uint64_t _xmath_emulu(unsigned int a, unsigned int b) {
    return (uint32_t) ((uint64_t) a) * ((uint64_t) b);
}

uint64_t xmath_residue_add(uint64_t x, uint64_t y) {
    uint64_t z = x + y;
    // z = z - (z >= MOD ? MOD : 0);
    if (z >= MOD)
        z -= MOD;
    return z;
}

uint64_t xmath_residue_sub(uint64_t x, uint64_t y) {
    uint64_t z = x - y;
    // z += (x < y ? MOD : 0);
    if (x < y)
        z += MOD;
    return z;
}

// copypasted from https://stackoverflow.com/questions/46870373/umul128-on-windows-32-bits
uint64_t xmath_umul128(uint64_t multiplier, uint64_t multiplicand, uint64_t* product_hi) {
    // multiplier   = ab = a * 2^32 + b
    // multiplicand = cd = c * 2^32 + d
    // ab * cd = a * c * 2^64 + (a * d + b * c) * 2^32 + b * d
    uint64_t a = multiplier >> 32;
    uint64_t b = (uint32_t) multiplier; // & 0xFFFFFFFF;
    uint64_t c = multiplicand >> 32;
    uint64_t d = (uint32_t) multiplicand; // & 0xFFFFFFFF;

    // uint64_t ac = _xmath_emulu(a, c);
    uint64_t ad = _xmath_emulu(a, d);
    // uint64_t bc = _xmath_emulu(b, c);
    uint64_t bd = _xmath_emulu(b, d);

    uint64_t adbc       = ad + _xmath_emulu(b, c);
    uint64_t adbc_carry = (adbc < ad); // ? 1 : 0;
    // MSVC gets confused by the ternary and makes worse code than using a boolean in an
    // integer context for 1 : 0

    // multiplier * multiplicand = product_hi * 2^64 + product_lo
    uint64_t product_lo       = bd + (adbc << 32);
    uint64_t product_lo_carry = (product_lo < bd); // ? 1 : 0;
    *product_hi =
        _xmath_emulu(a, c) + (adbc >> 32) + (adbc_carry << 32) + product_lo_carry;

    return product_lo;
}

uint64_t xmath_ui128_quotient_mod(uint64_t lo, uint64_t hi) {
    // hi:lo * ceil(2**170/MOD) >> (64 + 64 + 42)
    uint64_t prod1;
    xmath_umul128(lo, 0x604fa6a1c6346a87, &prod1);
    uint64_t part1hi;
    uint64_t part1lo = xmath_umul128(lo, 0x2d351c6d04f8b, &part1hi);
    uint64_t part2hi;
    uint64_t part2lo   = xmath_umul128(hi, 0x604fa6a1c6346a87, &part2hi);
    uint64_t sum1      = part1lo + part2lo;
    unsigned sum1carry = (sum1 < part1lo);
    sum1 += prod1;
    sum1carry += (sum1 < prod1);
    uint64_t prod2 = part1hi + part2hi + sum1carry;
    uint64_t prod3hi;
    uint64_t prod3lo = xmath_umul128(hi, 0x2d351c6d04f8b, &prod3hi);
    prod3lo += prod2;
    prod3hi += (prod3lo < prod2);
    return (prod3lo >> 42) | (prod3hi << 22);
}

uint64_t xmath_residue_mul(uint64_t x, uint64_t y) {
    // * ceil(2**170/MOD) = 0x2d351 c6d04f8b|604fa6a1 c6346a87 for (p-1)*(p-1) max
    uint64_t hi;
    uint64_t lo       = xmath_umul128(x, y, &hi);
    uint64_t quotient = xmath_ui128_quotient_mod(lo, hi);
    return lo - quotient * MOD;
}

uint64_t xmath_residue_pow(uint64_t x, uint64_t y) {
    if (y == 0)
        return 1;
    uint64_t cur = x;
    while (!(y & 1)) {
        cur = xmath_residue_mul(cur, cur);
        y >>= 1;
    }
    uint64_t res = cur;
    while ((y >>= 1) != 0) {
        cur = xmath_residue_mul(cur, cur);
        if (y & 1)
            res = xmath_residue_mul(res, cur);
    }
    return res;
}

uint64_t xmath_inverse(uint64_t u, uint64_t v) {
    int64_t  tmp;
    int64_t  xu = 1, xv = 0;
    uint64_t v0 = v;
    while (u > 1) {
        uint64_t d         = v / u;
        uint64_t remainder = v % u;
        tmp                = u;
        u                  = remainder;
        v                  = tmp;
        tmp                = xu;
        xu                 = xv - d * xu;
        xv                 = tmp;
    }
    xu += (xu < 0 ? v0 : 0);
    return xu;
}

uint64_t xmath_residue_inv(uint64_t x) {
    return xmath_inverse(x, MOD);
}
//{ return xmath_residue_pow(x, MOD - 2); }

uint64_t xmath_residue_sqrt(uint64_t what) {
    if (!what)
        return 0;
    uint64_t g = NON_RESIDUE, z, y, r, x, b, t;
    uint64_t e = 0, q = MOD - 1;
    while (!(q & 1))
        e++, q >>= 1;
    z = xmath_residue_pow(g, q);
    y = z;
    r = e;
    x = xmath_residue_pow(what, (q - 1) / 2);
    b = xmath_residue_mul(xmath_residue_mul(what, x), x);
    x = xmath_residue_mul(what, x);
    while (b != 1) {
        uint64_t m = 0, b2 = b;
        do {
            m++;
            b2 = xmath_residue_mul(b2, b2);
        } while (b2 != 1);
        if (m == r)
            return BAD;
        t = xmath_residue_pow(y, 1 << (r - m - 1));
        y = xmath_residue_mul(t, t);
        r = m;
        x = xmath_residue_mul(x, t);
        b = xmath_residue_mul(b, y);
    }
    if (xmath_residue_mul(x, x) != what) {
        // printf("internal error in sqrt\n");
        return BAD;
    }
    return x;
}
