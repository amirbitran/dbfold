#include "Random123/threefry.h"
#include "Random123/u01fixedpt.h"

/* 
 * threefryrand() generates random doubles in the range [0.0, 1.0]
 * using an algorithm from the Random123 library.
 *
 * Initialize by calling set_threefry_array with a user-supplied
 * seed (2^64 - 1 possible seeds). 
 *
 * For multi-processing situations, the user-supplied seed should
 * be different for each process if one seeks to avoid having
 * the same stream of random numbers for each thread/process.
 *
 */
#ifndef RNG_H_
#define RNG_H_

/* Sets RNG seed */
void set_threefry_array(uint64_t user_key);

/* Return double in range [0., 1.) */
double threefryrand();

/* Return random unsigned integer of 64 bits */
uint64_t threefryrand_int();


#endif /* RNG_H_ */
