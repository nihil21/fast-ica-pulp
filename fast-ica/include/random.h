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

#ifndef FAST_ICA_RANDOM_H
#define FAST_ICA_RANDOM_H

#include "fp.h"

void set_prng_seed(unsigned int seed);
fp uniform_rand();
fp uniform_rand_range(fp min, fp max);
int uniform_randint_range(int min, int max);
fp standard_rand();

#endif //FAST_ICA_RANDOM_H
