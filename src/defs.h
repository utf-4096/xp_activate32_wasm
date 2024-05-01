#ifndef DEFS_H
#define DEFS_H
#include <stdint.h>

typedef uintmax_t size_t;

#define MOD         0x16A6B036D7F2A79ULL
#define NON_RESIDUE 43
#define BAD         0xFFFFFFFFFFFFFFFFull

typedef enum {
    ERR_NONE                = 0,
    ERR_TOO_SHORT           = 1,
    ERR_TOO_LARGE           = 2,
    ERR_INVALID_CHARACTER   = 3,
    ERR_INVALID_CHECK_DIGIT = 4,
    ERR_UNKNOWN_VERSION     = 5,
    ERR_UNLUCKY             = 6,
} error_t;

#pragma pack(push, 1)
typedef struct
{
    uint64_t       HardwareID;
    uint64_t       ProductIDLow;
    unsigned char  ProductIDHigh;
    unsigned short KeySHA1;
} key_t;
#pragma pack(pop)

typedef struct
{
    uint64_t u[2];
    uint64_t v[2];
} TDivisor;

#endif
