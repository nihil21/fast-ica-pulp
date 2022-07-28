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

#ifndef FAST_ICA_FAST_ICA_H
#define FAST_ICA_FAST_ICA_H

#include "../include/matrix.h"

typedef enum FastICAStrategy {
    Parallel,
    Deflation
} FastICAStrategy;

typedef enum GFunc {
    LogCosh,
    Exp,
    Cube
} GFunc;

typedef struct FastICAArgs {
    Matrix *res;
    Matrix *x;
    int n_components;
    bool whiten;
    FastICAStrategy strategy;
    GFunc g_func;
    fp threshold;
    int max_iter;
} FastICAArgs;

void fast_ica(FastICAArgs *arg);

#endif //FAST_ICA_FAST_ICA_H
