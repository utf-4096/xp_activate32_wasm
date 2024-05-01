#include "defs.h"
#include "mem.h"
#include "xmath.h"

#include <stdint.h>

static const uint64_t f[6] = {0,
                              0x21840136C85381ULL,
                              0x44197B83892AD0ULL,
                              0x1400606322B3B04ULL,
                              0x1400606322B3B04ULL,
                              1};

#define divisor_double(src, dst) \
    divisor_add(src, src, dst)

static int find_divisor_v(TDivisor* d) {
    // u | v^2 - f
    // u = u0 + u1*x + x^2
    // f%u = f0 + f1*x
    uint64_t v1;
    uint64_t f2[6];
    int      i, j;

    for (i = 0; i < 6; i++)
        f2[i] = f[i];

    const uint64_t u0 = d->u[0];
    const uint64_t u1 = d->u[1];
    for (j = 4; j--;) {
        f2[j]     = xmath_residue_sub(f2[j], xmath_residue_mul(u0, f2[j + 2]));
        f2[j + 1] = xmath_residue_sub(f2[j + 1], xmath_residue_mul(u1, f2[j + 2]));
        f2[j + 2] = 0;
    }

    // v = v0 + v1*x
    // u | (v0^2 - f0) + (2*v0*v1 - f1)*x + v1^2*x^2 = u0*v1^2 + u1*v1^2*x + v1^2*x^2
    // v0^2 - f0 = u0*v1^2
    // 2*v0*v1 - f1 = u1*v1^2
    // v0^2 = f0 + u0*v1^2 = (f1 + u1*v1^2)^2 / (2*v1)^2
    // (f1^2) + 2*(f1*u1-2*f0) * v1^2 + (u1^2-4*u0) * v1^4 = 0
    // v1^2 = ((2*f0-f1*u1) +- 2*sqrt(-f0*f1*u1 + f0^2 + f1^2*u0))) / (u1^2-4*u0)
    const uint64_t f0       = f2[0];
    const uint64_t f1       = f2[1];
    const uint64_t u0double = xmath_residue_add(u0, u0);
    const uint64_t coeff2   = xmath_residue_sub(xmath_residue_mul(u1, u1),
                                              xmath_residue_add(u0double, u0double));
    const uint64_t coeff1 =
        xmath_residue_sub(xmath_residue_add(f0, f0), xmath_residue_mul(f1, u1));

    if (coeff2 == 0) {
        if (coeff1 == 0) {
            if (f1 == 0) {
                // impossible
                // printf("bad f(), double root detected\n");
            }
            return 0;
        }

        uint64_t sqr =
            xmath_residue_mul(xmath_residue_mul(f1, f1),
                              xmath_residue_inv(xmath_residue_add(coeff1, coeff1)));

        v1 = xmath_residue_sqrt(sqr);
        if (v1 == BAD)
            return 0;
    } else {
        uint64_t d = xmath_residue_add(
            xmath_residue_mul(f0, f0),
            xmath_residue_mul(
                f1,
                xmath_residue_sub(xmath_residue_mul(f1, u0), xmath_residue_mul(f0, u1))));

        d = xmath_residue_sqrt(d);
        if (d == BAD)
            return 0;

        d             = xmath_residue_add(d, d);
        uint64_t inv  = xmath_residue_inv(coeff2);
        uint64_t root = xmath_residue_mul(xmath_residue_add(coeff1, d), inv);
        v1            = xmath_residue_sqrt(root);

        if (v1 == BAD) {
            root = xmath_residue_mul(xmath_residue_sub(coeff1, d), inv);
            v1   = xmath_residue_sqrt(root);
            if (v1 == BAD)
                return 0;
        }
    }
    uint64_t v0 = xmath_residue_mul(
        xmath_residue_add(f1, xmath_residue_mul(u1, xmath_residue_mul(v1, v1))),
        xmath_residue_inv(xmath_residue_add(v1, v1)));

    d->v[0] = v0;
    d->v[1] = v1;
    return 1;
}

// generic short slow code
static int polynomial_mul(int            adeg,
                          const uint64_t a[],
                          int            bdeg,
                          const uint64_t b[],
                          int            resultprevdeg,
                          uint64_t       result[]) {
    if (adeg < 0 || bdeg < 0)
        return resultprevdeg;

    int i, j;
    for (i = resultprevdeg + 1; i <= adeg + bdeg; i++)
        result[i] = 0;

    resultprevdeg = i - 1;
    for (i = 0; i <= adeg; i++)
        for (j = 0; j <= bdeg; j++)
            result[i + j] =
                xmath_residue_add(result[i + j], xmath_residue_mul(a[i], b[j]));

    while (resultprevdeg >= 0 && result[resultprevdeg] == 0)
        --resultprevdeg;

    return resultprevdeg;
}

static int polynomial_div_monic(int            adeg,
                                uint64_t       a[],
                                int            bdeg,
                                const uint64_t b[],
                                uint64_t*      quotient) {
    int i, j;
    for (i = adeg - bdeg; i >= 0; i--) {
        uint64_t q = a[i + bdeg];
        if (quotient)
            quotient[i] = q;

        for (j = 0; j < bdeg; j++)
            a[i + j] = xmath_residue_sub(a[i + j], xmath_residue_mul(q, b[j]));

        a[i + j] = 0;
    }

    i += bdeg;
    while (i >= 0 && a[i] == 0)
        i--;

    return i;
}

static void polynomial_xgcd(int            adeg,
                            const uint64_t a[3],
                            int            bdeg,
                            const uint64_t b[3],
                            int*           pgcddeg,
                            uint64_t       gcd[3],
                            int*           pmult1deg,
                            uint64_t       mult1[3],
                            int*           pmult2deg,
                            uint64_t       mult2[3]) {
    int      sdeg     = -1;
    uint64_t s[3]     = {0, 0, 0};
    int      mult1deg = 0;
    mult1[0]          = 1;
    mult1[1]          = 0;
    mult1[2]          = 0;
    int      tdeg     = 0;
    uint64_t t[3]     = {1, 0, 0};
    int      mult2deg = -1;
    mult2[0]          = 0;
    mult2[1]          = 0;
    mult2[2]          = 0;
    int      rdeg     = bdeg;
    uint64_t r[3]     = {b[0], b[1], b[2]};
    int      gcddeg   = adeg;
    gcd[0]            = a[0];
    gcd[1]            = a[1];
    gcd[2]            = a[2];

    // s*u1 + t*u2 = r
    // mult1*u1 + mult2*u2 = gcd
    while (rdeg >= 0) {
        if (rdeg > gcddeg) {
            unsigned tmp;
            int      tmpi;
            tmp      = rdeg;
            rdeg     = gcddeg;
            gcddeg   = tmp;
            tmpi     = sdeg;
            sdeg     = mult1deg;
            mult1deg = tmpi;
            tmpi     = tdeg;
            tdeg     = mult2deg;
            mult2deg = tmpi;
            uint64_t tmp2;
            tmp2     = r[0];
            r[0]     = gcd[0];
            gcd[0]   = tmp2;
            tmp2     = r[1];
            r[1]     = gcd[1];
            gcd[1]   = tmp2;
            tmp2     = r[2];
            r[2]     = gcd[2];
            gcd[2]   = tmp2;
            tmp2     = s[0];
            s[0]     = mult1[0];
            mult1[0] = tmp2;
            tmp2     = s[1];
            s[1]     = mult1[1];
            mult1[1] = tmp2;
            tmp2     = s[2];
            s[2]     = mult1[2];
            mult1[2] = tmp2;
            tmp2     = t[0];
            t[0]     = mult2[0];
            mult2[0] = tmp2;
            tmp2     = t[1];
            t[1]     = mult2[1];
            mult2[1] = tmp2;
            tmp2     = t[2];
            t[2]     = mult2[2];
            mult2[2] = tmp2;
            continue;
        }

        int      delta = gcddeg - rdeg;
        uint64_t mult  = xmath_residue_mul(gcd[gcddeg], xmath_residue_inv(r[rdeg]));

        // quotient = mult * x**delta
        for (int i = 0; i <= rdeg; i++)
            gcd[i + delta] =
                xmath_residue_sub(gcd[i + delta], xmath_residue_mul(mult, r[i]));

        while (gcddeg >= 0 && gcd[gcddeg] == 0)
            gcddeg--;

        for (int i = 0; i <= sdeg; i++)
            mult1[i + delta] =
                xmath_residue_sub(mult1[i + delta], xmath_residue_mul(mult, s[i]));

        if (mult1deg < sdeg + delta)
            mult1deg = sdeg + delta;

        while (mult1deg >= 0 && mult1[mult1deg] == 0)
            mult1deg--;

        for (int i = 0; i <= tdeg; i++)
            mult2[i + delta] =
                xmath_residue_sub(mult2[i + delta], xmath_residue_mul(mult, t[i]));

        if (mult2deg < tdeg + delta)
            mult2deg = tdeg + delta;

        while (mult2deg >= 0 && mult2[mult2deg] == 0)
            mult2deg--;
    }

    // d1 = gcd, e1 = mult1, e2 = mult2
    *pgcddeg   = gcddeg;
    *pmult1deg = mult1deg;
    *pmult2deg = mult2deg;
}

static int u2poly(const TDivisor* src, uint64_t polyu[3], uint64_t polyv[2]) {
    if (src->u[1] != BAD) {
        polyu[0] = src->u[0];
        polyu[1] = src->u[1];
        polyu[2] = 1;
        polyv[0] = src->v[0];
        polyv[1] = src->v[1];
        return 2;
    }

    if (src->u[0] != BAD) {
        polyu[0] = src->u[0];
        polyu[1] = 1;
        polyv[0] = src->v[0];
        polyv[1] = 0;
        return 1;
    }

    polyu[0] = 1;
    polyv[0] = 0;
    polyv[1] = 0;
    return 0;
}

static void divisor_add(const TDivisor* src1, const TDivisor* src2, TDivisor* dst) {
    uint64_t u1[3], u2[3], v1[2], v2[2];
    int      u1deg = u2poly(src1, u1, v1);
    int      u2deg = u2poly(src2, u2, v2);
    // extended gcd: d1 = gcd(u1, u2) = e1*u1 + e2*u2
    int      d1deg, e1deg, e2deg;
    uint64_t d1[3], e1[3], e2[3];
    polynomial_xgcd(u1deg, u1, u2deg, u2, &d1deg, d1, &e1deg, e1, &e2deg, e2);

    // extended gcd again: d = gcd(d1, v1+v2) = c1*d1 + c2*(v1+v2)
    uint64_t b[3] = {xmath_residue_add(v1[0], v2[0]), xmath_residue_add(v1[1], v2[1]), 0};
    int      bdeg = (b[1] == 0 ? (b[0] == 0 ? -1 : 0) : 1);
    int      ddeg, c1deg, c2deg;
    uint64_t d[3], c1[3], c2[3];
    polynomial_xgcd(d1deg, d1, bdeg, b, &ddeg, d, &c1deg, c1, &c2deg, c2);

    uint64_t dmult = xmath_residue_inv(d[ddeg]);
    int      i;
    for (i = 0; i < ddeg; i++)
        d[i] = xmath_residue_mul(d[i], dmult);

    d[i] = 1;
    for (i = 0; i <= c1deg; i++)
        c1[i] = xmath_residue_mul(c1[i], dmult);

    for (i = 0; i <= c2deg; i++)
        c2[i] = xmath_residue_mul(c2[i], dmult);

    uint64_t u[5];
    int      udeg = polynomial_mul(u1deg, u1, u2deg, u2, -1, u);
    // u is monic
    uint64_t v[7], tmp[7];
    int      vdeg, tmpdeg;
    // c1*(e1*u1*v2 + e2*u2*v1) + c2*(v1*v2 + f)
    // c1*(e1*u1*(v2-v1) + d1*v1) + c2*(v1*v2 + f)
    v[0]   = xmath_residue_sub(v2[0], v1[0]);
    v[1]   = xmath_residue_sub(v2[1], v1[1]);
    tmpdeg = polynomial_mul(e1deg, e1, 1, v, -1, tmp);
    vdeg   = polynomial_mul(u1deg, u1, tmpdeg, tmp, -1, v);
    vdeg   = polynomial_mul(d1deg, d1, 1, v1, vdeg, v);

    for (i = 0; i <= vdeg; i++)
        v[i] = xmath_residue_mul(v[i], c1[0]);

    memcpy(tmp, f, 6 * sizeof(f[0]));
    tmpdeg = 5;
    tmpdeg = polynomial_mul(1, v1, 1, v2, tmpdeg, tmp);
    vdeg   = polynomial_mul(c2deg, c2, tmpdeg, tmp, vdeg, v);

    if (ddeg > 0) {
        uint64_t udiv[5];
        polynomial_div_monic(udeg, u, ddeg, d, udiv);
        udeg -= ddeg;
        polynomial_div_monic(udeg, udiv, ddeg, d, u);
        udeg -= ddeg;
        if (vdeg >= 0) {
            polynomial_div_monic(vdeg, v, ddeg, d, udiv);
            vdeg -= ddeg;
            memcpy(v, udiv, (vdeg + 1) * sizeof(v[0]));
        }
    }

    vdeg = polynomial_div_monic(vdeg, v, udeg, u, NULL);
    while (udeg > 2) {
        // u' = monic((f-v^2)/u), v'=-v mod u'
        tmpdeg = polynomial_mul(vdeg, v, vdeg, v, -1, tmp);
        for (i = 0; i <= tmpdeg && i <= 5; i++)
            tmp[i] = xmath_residue_sub(f[i], tmp[i]);
        for (; i <= tmpdeg; i++)
            tmp[i] = xmath_residue_sub(0, tmp[i]);
        for (; i <= 5; i++)
            tmp[i] = f[i];
        tmpdeg = i - 1;
        uint64_t udiv[5];
        polynomial_div_monic(tmpdeg, tmp, udeg, u, udiv);
        udeg          = tmpdeg - udeg;
        uint64_t mult = xmath_residue_inv(udiv[udeg]);
        for (i = 0; i < udeg; i++)
            u[i] = xmath_residue_mul(udiv[i], mult);
        u[i] = 1;
        for (i = 0; i <= vdeg; i++)
            v[i] = xmath_residue_sub(0, v[i]);
        vdeg = polynomial_div_monic(vdeg, v, udeg, u, NULL);
    }

    if (udeg == 2) {
        dst->u[0] = u[0];
        dst->u[1] = u[1];
        dst->v[0] = (vdeg >= 0 ? v[0] : 0);
        dst->v[1] = (vdeg >= 1 ? v[1] : 0);
    } else if (udeg == 1) {
        dst->u[0] = u[0];
        dst->u[1] = BAD;
        dst->v[0] = (vdeg >= 0 ? v[0] : 0);
        dst->v[1] = BAD;
    } else {
        dst->u[0] = BAD;
        dst->u[1] = BAD;
        dst->v[0] = BAD;
        dst->v[1] = BAD;
    }
}

static void divisor_mul(const TDivisor* src, uint64_t mult, TDivisor* dst) {
    if (mult == 0) {
        dst->u[0] = BAD;
        dst->u[1] = BAD;
        dst->v[0] = BAD;
        dst->v[1] = BAD;
        return;
    }

    TDivisor cur = *src;
    while (!(mult & 1)) {
        divisor_double(&cur, &cur);
        mult >>= 1;
    }

    *dst = cur;
    while ((mult >>= 1) != 0) {
        divisor_double(&cur, &cur);
        if (mult & 1)
            divisor_add(dst, &cur, dst);
    }
}

static void
divisor_mul128(const TDivisor* src, uint64_t mult_lo, uint64_t mult_hi, TDivisor* dst) {
    if (mult_lo == 0 && mult_hi == 0) {
        dst->u[0] = BAD;
        dst->u[1] = BAD;
        dst->v[0] = BAD;
        dst->v[1] = BAD;
        return;
    }

    TDivisor cur = *src;
    while (!(mult_lo & 1)) {
        divisor_double(&cur, &cur);
        mult_lo >>= 1;
        if (mult_hi & 1)
            mult_lo |= (1ULL << 63);
        mult_hi >>= 1;
    }

    *dst = cur;
    for (;;) {
        mult_lo >>= 1;
        if (mult_hi & 1)
            mult_lo |= (1ULL << 63);
        mult_hi >>= 1;
        if (mult_lo == 0 && mult_hi == 0)
            break;
        divisor_double(&cur, &cur);
        if (mult_lo & 1)
            divisor_add(dst, &cur, dst);
    }
}

static unsigned rol(unsigned x, int shift) {
    return (x << shift) | (x >> (32 - shift));
}

static void sha1_single_block(unsigned char input[64], unsigned char output[20]) {
    unsigned a, b, c, d, e;
    a = 0x67452301;
    b = 0xEFCDAB89;
    c = 0x98BADCFE;
    d = 0x10325476;
    e = 0xC3D2E1F0;

    unsigned w[80];
    size_t   i;
    for (i = 0; i < 16; i++)
        w[i] = input[4 * i] << 24 | input[4 * i + 1] << 16 | input[4 * i + 2] << 8 |
               input[4 * i + 3];

    for (i = 16; i < 80; i++)
        w[i] = rol(w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16], 1);

    for (i = 0; i < 20; i++) {
        unsigned tmp = rol(a, 5) + ((b & c) | (~b & d)) + e + w[i] + 0x5A827999;
        e            = d;
        d            = c;
        c            = rol(b, 30);
        b            = a;
        a            = tmp;
    }

    for (i = 20; i < 40; i++) {
        unsigned tmp = rol(a, 5) + (b ^ c ^ d) + e + w[i] + 0x6ED9EBA1;
        e            = d;
        d            = c;
        c            = rol(b, 30);
        b            = a;
        a            = tmp;
    }

    for (i = 40; i < 60; i++) {
        unsigned tmp = rol(a, 5) + ((b & c) | (b & d) | (c & d)) + e + w[i] + 0x8F1BBCDC;
        e            = d;
        d            = c;
        c            = rol(b, 30);
        b            = a;
        a            = tmp;
    }

    for (i = 60; i < 80; i++) {
        unsigned tmp = rol(a, 5) + (b ^ c ^ d) + e + w[i] + 0xCA62C1D6;
        e            = d;
        d            = c;
        c            = rol(b, 30);
        b            = a;
        a            = tmp;
    }

    a += 0x67452301;
    b += 0xEFCDAB89;
    c += 0x98BADCFE;
    d += 0x10325476;
    e += 0xC3D2E1F0;

    output[0]  = a >> 24;
    output[1]  = a >> 16;
    output[2]  = a >> 8;
    output[3]  = a;
    output[4]  = b >> 24;
    output[5]  = b >> 16;
    output[6]  = b >> 8;
    output[7]  = b;
    output[8]  = c >> 24;
    output[9]  = c >> 16;
    output[10] = c >> 8;
    output[11] = c;
    output[12] = d >> 24;
    output[13] = d >> 16;
    output[14] = d >> 8;
    output[15] = d;
    output[16] = e >> 24;
    output[17] = e >> 16;
    output[18] = e >> 8;
    output[19] = e;
}

static void
mix(unsigned char* buffer, size_t bufSize, const unsigned char* key, size_t keySize) {
    unsigned char sha1_input[64];
    unsigned char sha1_result[20];
    size_t        half = bufSize / 2;
    int           external_counter;
    for (external_counter = 0; external_counter < 4; external_counter++) {
        memset(sha1_input, 0, sizeof(sha1_input));
        memcpy(sha1_input, buffer + half, half);
        memcpy(sha1_input + half, key, keySize);
        sha1_input[half + keySize]         = 0x80;
        sha1_input[sizeof(sha1_input) - 1] = (half + keySize) * 8;
        sha1_input[sizeof(sha1_input) - 2] = (half + keySize) * 8 / 0x100;
        sha1_single_block(sha1_input, sha1_result);

        size_t i;
        for (i = half & ~3; i < half; i++)
            sha1_result[i] = sha1_result[i + 4 - (half & 3)];

        for (i = 0; i < half; i++) {
            unsigned char tmp = buffer[i + half];
            buffer[i + half]  = buffer[i] ^ sha1_result[i];
            buffer[i]         = tmp;
        }
    }
}

static void
unmix(unsigned char* buffer, size_t bufSize, const unsigned char* key, size_t keySize) {
    unsigned char sha1_input[64];
    unsigned char sha1_result[20];
    size_t        half = bufSize / 2;
    int           external_counter;

    for (external_counter = 0; external_counter < 4; external_counter++) {
        memset(sha1_input, 0, sizeof(sha1_input));
        memcpy(sha1_input, buffer, half);
        memcpy(sha1_input + half, key, keySize);
        sha1_input[half + keySize]         = 0x80;
        sha1_input[sizeof(sha1_input) - 1] = (half + keySize) * 8;
        sha1_input[sizeof(sha1_input) - 2] = (half + keySize) * 8 / 0x100;
        sha1_single_block(sha1_input, sha1_result);
        size_t i;

        for (i = half & ~3; i < half; i++)
            sha1_result[i] = sha1_result[i + 4 - (half & 3)];

        for (i = 0; i < half; i++) {
            unsigned char tmp = buffer[i];
            buffer[i]         = buffer[i + half] ^ sha1_result[i];
            buffer[i + half]  = tmp;
        }
    }
}

int generate(const char* installation_id_str, char confirmation_id[49]) {
    unsigned char installation_id[19]; // 10**45 < 256**19
    size_t        installation_id_len = 0;
    const char*   p                   = installation_id_str;
    size_t        count = 0, totalCount = 0;
    unsigned      check = 0;
    size_t        i;

    for (; *p; p++) {
        if (*p == ' ' || *p == '-')
            continue;

        int d = *p - '0';
        if (d < 0 || d > 9)
            return ERR_INVALID_CHARACTER;

        if (count == 5 || p[1] == 0) {
            if (!count)
                return (totalCount == 45) ? ERR_TOO_LARGE : ERR_TOO_SHORT;
            if (d != check % 7)
                return (count < 5) ? ERR_TOO_SHORT : ERR_INVALID_CHECK_DIGIT;
            check = 0;
            count = 0;
            continue;
        }

        check += (count % 2 ? d * 2 : d);
        count++;
        totalCount++;

        if (totalCount > 45)
            return ERR_TOO_LARGE;

        unsigned char carry = d;
        for (i = 0; i < installation_id_len; i++) {
            unsigned x         = installation_id[i] * 10 + carry;
            installation_id[i] = x & 0xFF;
            carry              = x >> 8;
        }

        if (carry) {
            installation_id[installation_id_len++] = carry;
        }
    }

    if (totalCount != 41 && totalCount < 45)
        return ERR_TOO_SHORT;

    for (; installation_id_len < sizeof(installation_id); installation_id_len++)
        installation_id[installation_id_len] = 0;

    static const unsigned char iid_key[4] = {0x6A, 0xC8, 0x5E, 0xD4};
    unmix(installation_id, totalCount == 41 ? 17 : 19, iid_key, 4);

    if (installation_id[18] >= 0x10)
        return ERR_UNKNOWN_VERSION;

    key_t parsed;
    memcpy(&parsed, installation_id, sizeof(parsed));
    unsigned productId1 = parsed.ProductIDLow & ((1 << 17) - 1);
    unsigned productId2 = (parsed.ProductIDLow >> 17) & ((1 << 10) - 1);
    unsigned productId3 = (parsed.ProductIDLow >> 27) & ((1 << 25) - 1);
    unsigned version    = (parsed.ProductIDLow >> 52) & 7;
    unsigned productId4 = (parsed.ProductIDLow >> 55) | (parsed.ProductIDHigh << 9);

    if (version != (totalCount == 41 ? 4 : 5))
        return ERR_UNKNOWN_VERSION;

    // printf("Product ID: %05u-%03u-%07u-%05u\n", productId1, productId2, productId3,
    // productId4);

    unsigned char keybuf[16];
    memcpy(keybuf, &parsed.HardwareID, 8);
    uint64_t productIdMixed = (uint64_t) productId1 << 41 | (uint64_t) productId2 << 58 |
                              (uint64_t) productId3 << 17 | productId4;
    memcpy(keybuf + 8, &productIdMixed, 8);

    TDivisor      d;
    unsigned char attempt;
    for (attempt = 0; attempt <= 0x80; attempt++) {
        union
        {
            unsigned char buffer[14];
            struct
            {
                uint64_t lo;
                uint64_t hi;
            };
        } u;
        u.lo        = 0;
        u.hi        = 0;
        u.buffer[7] = attempt;

        mix(u.buffer, 14, keybuf, 16);

        uint64_t x2 = xmath_ui128_quotient_mod(u.lo, u.hi);
        uint64_t x1 = u.lo - x2 * MOD;
        x2++;
        d.u[0] =
            xmath_residue_sub(xmath_residue_mul(x1, x1),
                              xmath_residue_mul(NON_RESIDUE, xmath_residue_mul(x2, x2)));
        d.u[1] = xmath_residue_add(x1, x1);

        if (find_divisor_v(&d))
            break;
    }

    if (attempt > 0x80)
        return ERR_UNLUCKY;

    divisor_mul128(&d, 0x04e21b9d10f127c1, 0x40da7c36d44c, &d);
    union
    {
        struct
        {
            uint64_t encoded_lo, encoded_hi;
        };
        struct
        {
            uint32_t encoded[4];
        };
    } e;

    if (d.u[0] == BAD) {
        // we can not get the zero divisor, actually...
        e.encoded_lo = xmath_umul128(MOD + 2, MOD, &e.encoded_hi);
    } else if (d.u[1] == BAD) {
        // O(1/MOD) chance
        // encoded = (unsigned __int128)(MOD + 1) * d.u[0] + MOD; // * MOD + d.u[0] is
        // fine too
        e.encoded_lo = xmath_umul128(MOD + 1, d.u[0], &e.encoded_hi);
        e.encoded_lo += MOD;
        e.encoded_hi += (e.encoded_lo < MOD);
    } else {
        uint64_t x1    = (d.u[1] % 2 ? d.u[1] + MOD : d.u[1]) / 2;
        uint64_t x2sqr = xmath_residue_sub(xmath_residue_mul(x1, x1), d.u[0]);
        uint64_t x2    = xmath_residue_sqrt(x2sqr);

        if (x2 == BAD) {
            x2 = xmath_residue_sqrt(
                xmath_residue_mul(x2sqr, xmath_residue_inv(NON_RESIDUE)));
            e.encoded_lo = xmath_umul128(MOD + 1, MOD + x2, &e.encoded_hi);
            e.encoded_lo += x1;
            e.encoded_hi += (e.encoded_lo < x1);
        } else {
            // points (-x1+x2, v(-x1+x2)) and (-x1-x2, v(-x1-x2))
            uint64_t x1a = xmath_residue_sub(x1, x2);
            uint64_t y1  = xmath_residue_sub(d.v[0], xmath_residue_mul(d.v[1], x1a));
            uint64_t x2a = xmath_residue_add(x1, x2);
            uint64_t y2  = xmath_residue_sub(d.v[0], xmath_residue_mul(d.v[1], x2a));

            if (x1a > x2a) {
                uint64_t tmp = x1a;
                x1a          = x2a;
                x2a          = tmp;
            }

            if ((y1 ^ y2) & 1) {
                uint64_t tmp = x1a;
                x1a          = x2a;
                x2a          = tmp;
            }

            e.encoded_lo = xmath_umul128(MOD + 1, x1a, &e.encoded_hi);
            e.encoded_lo += x2a;
            e.encoded_hi += (e.encoded_lo < x2a);
        }
    }

    unsigned char decimal[35];
    for (i = 0; i < 35; i++) {
        unsigned c = e.encoded[3] % 10;
        e.encoded[3] /= 10;
        unsigned c2     = ((uint64_t) c << 32 | e.encoded[2]) % 10;
        e.encoded[2]    = ((uint64_t) c << 32 | e.encoded[2]) / 10;
        unsigned c3     = ((uint64_t) c2 << 32 | e.encoded[1]) % 10;
        e.encoded[1]    = ((uint64_t) c2 << 32 | e.encoded[1]) / 10;
        unsigned c4     = ((uint64_t) c3 << 32 | e.encoded[0]) % 10;
        e.encoded[0]    = ((uint64_t) c3 << 32 | e.encoded[0]) / 10;
        decimal[34 - i] = c4;
    }

    char* q = confirmation_id;
    for (i = 0; i < 7; i++) {
        if (i)
            *q++ = '-';
        unsigned char* p = decimal + i * 5;
        q[0]             = p[0] + '0';
        q[1]             = p[1] + '0';
        q[2]             = p[2] + '0';
        q[3]             = p[3] + '0';
        q[4]             = p[4] + '0';
        q[5]             = ((p[0] + p[1] * 2 + p[2] + p[3] * 2 + p[4]) % 7) + '0';
        q += 6;
    }

    *q++ = 0;
    return ERR_NONE;
}
