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

#ifndef FAST_ICA_FP_H
#define FAST_ICA_FP_H

#include <math.h>

/*
 * Data type declaration and macros
 */
typedef float fp;
typedef unsigned int uint;

extern const fp PI;

#define SGN(x) (((x) >= 0) ? 1 : -1)
#define MAX(x, y) (((x) >= (y)) ? x : y)
#define POW powf
#define SQRT sqrtf
#define LOG logf
#define EXP expf
#define SIN sinf
#define TANH tanhf
#define ABS fabsf
#define MOD fmodf

#endif //FAST_ICA_FP_H
