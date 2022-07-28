/*
Copyright 2022 Mattia Orlandi

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "../include/random.h"

const unsigned int MAX_INT = 4294967295;

// Initialize seed
static unsigned int xs = 42;

/*
 * Set seed of PRNG
 */
void set_prng_seed(unsigned int seed) {
    xs = seed;
}

/*
 * Uniform PRNG for integers (Marsaglia's Xorshift)
 */
unsigned int uniform_randint() {
    xs ^= xs << 13;
    xs ^= xs >> 17;
    xs ^= xs << 5;
    return xs;
}

/*
 * Uniform PRNG for integers in given range
 */
int uniform_randint_range(int min, int max) {
    int r = (int) uniform_randint();
    // Take absolute value after cast
    if (r < 0)
        r = -r;
    return r % (max - min + 1) + min;
}

/*
 * Uniform PRNG for floating points between 0 and 1
 */
fp uniform_rand() {
    return (fp) uniform_randint() / (fp) MAX_INT;
}

/*
 * Uniform PRNG for floating points in given range
 */
fp uniform_rand_range(fp min, fp max) {
    fp r = uniform_rand();
    return r * (max - min) + min;
}

/*
 * Standard normal PRNG (Marsaglia's polar method)
 */
fp standard_rand() {
    static fp cached = 0;
    fp r;

    if (cached != 0) {  // result cached
        r = cached;
        cached = 0;
    } else {
        fp u, v, s;
        // Generate two uniform values between -1 and 1 s.t. the sum of their squares is less than 1
        do {
            u = (fp) uniform_rand() * 2 - 1;
            v = (fp) uniform_rand() * 2 - 1;
            s = u * u + v * v;
        } while (s >= 1);
        // Obtain two normal_mat variables
        r = u * SQRT((-2 * LOG(s)) / s);
        cached = v * SQRT((-2 * LOG(s)) / s);
    }

    return r;
}
