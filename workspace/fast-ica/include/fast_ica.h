//
// Created by nihil on 08/10/21.
//

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
    Matrix **res;
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
