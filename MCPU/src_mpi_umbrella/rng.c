#include "rng.h"


static threefry4x64_ctr_t ctr = {{}};
static threefry4x64_key_t key = {{}};

static unsigned char randomNumberIndex = 4;
static threefry4x64_ctr_t result = {{}};


void increment_counter()
{
    static unsigned char index = 0;
    if (index == 4) {
        index = 0;
    }
    ctr.v[index++]++;
}


/* Functions implementations for rng.h */

void set_threefry_array(uint64_t user_key)
{
    key.v[0] = user_key;
    key.v[1] = user_key;
    key.v[2] = user_key;
    key.v[3] = user_key;
}

double threefryrand()
{
    if (randomNumberIndex == 4)
    {
        result = threefry4x64(ctr, key);
        randomNumberIndex = 0;
        increment_counter();
    }
    return u01fixedpt_closed_open_64_53(result.v[randomNumberIndex++]);
}

uint64_t threefryrand_int()
{
   if (randomNumberIndex == 4)
   {
        result = threefry4x64(ctr, key);
        randomNumberIndex = 0;
        increment_counter();
    }
    return result.v[randomNumberIndex++];
}
