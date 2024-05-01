#ifndef MEM_H
#define MEM_H
#include "defs.h"

#define NULL ((void*) 0)
void* memcpy(void* dest, const void* src, unsigned long size);
void* memset(void* dest, int data, unsigned long size);

#endif
