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

#ifndef FAST_ICA_SIGNAL_H
#define FAST_ICA_SIGNAL_H

#include "../include/matrix.h"

/*
 * Math operations (in-place)
 */
void sine_mat_(Matrix *m);
void sgn_mat_(Matrix *m);
void mod_mat_(Matrix *m, fp k);

/*
 * Signals
 */
void generate_sine_wave(Matrix *s, fp amp, fp freq, fp phase, int n_samples, fp range);
void generate_square_wave(Matrix *s, fp amp, fp freq, fp phase, int n_samples, fp range);
void generate_sawtooth_wave(Matrix *s, fp amp, fp freq, fp phase, int n_samples, fp range);
void generate_signals(Matrix *s, int n_samples, fp range, bool noise);

#endif //FAST_ICA_SIGNAL_H
