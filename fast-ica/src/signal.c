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

#include "../include/signal.h"

/*
 * Apply sine function to a matrix (in-place)
 */
void sine_mat_(Matrix *m) {
    for (int i = 0; i < m->height; i++) {
        for (int j = 0; j < m->width; j++) {
            MAT_CELL(m, i, j) = SIN(MAT_CELL(m, i, j));
        }
    }
}

/*
 * Apply SGN function to a matrix (in-place)
 */
void sgn_mat_(Matrix *m) {
    for (int i = 0; i < m->height; i++) {
        for (int j = 0; j < m->width; j++) {
            MAT_CELL(m, i, j) = (fp) SGN(MAT_CELL(m, i, j));
        }
    }
}

/*
 * Apply modulus function to a matrix (in-place)
 */
void mod_mat_(Matrix *m, fp k) {
    for (int i = 0; i < m->height; i++) {
        for (int j = 0; j < m->width; j++) {
            MAT_CELL(m, i, j) = MOD(MAT_CELL(m, i, j), k);
        }
    }
}

void generate_sine_wave(Matrix *s, fp amp, fp freq, fp phase, int n_samples, fp range) {
    // Generate a vector with n_samples values from 0 to range, representing time slices
    linspace(s, 0, range, n_samples);
    // Multiply it by angular velocity and translate it by phase (in-place)
    fp omega = 2 * PI * freq;
    scale_(s, omega);
    add_scalar_(s, phase);

    // Generate sine wave (in-place)
    sine_mat_(s);
    // Multiply by amplitude (in-place)
    scale_(s, amp);
}

void generate_square_wave(Matrix *s, fp amp, fp freq, fp phase, int n_samples, fp range) {
    // Generate a vector with n_samples values from 0 to range, representing time slices
    linspace(s, 0, range, n_samples);
    // Multiply it by angular velocity and translate it by phase (in-place)
    fp omega = 2 * PI * freq;
    scale_(s, omega);
    add_scalar_(s, phase);

    // Generate sine wave and take sign (in-place)
    sine_mat_(s);
    sgn_mat_(s);
    // Multiply by amplitude (in-place)
    scale_(s, amp);
}

void generate_sawtooth_wave(Matrix *s, fp amp, fp freq, fp phase, int n_samples, fp range) {
    // Generate a vector with n_samples values from 0 to range, representing time slices
    linspace(s, 0, range, n_samples);
    // Multiply it by angular velocity and translate it by phase (in-place)
    fp omega = 2 * PI * freq;
    scale_(s, omega);
    add_scalar_(s, phase);

    // Generate sawtooth wave (in-place)
    mod_mat_(s, 2 * PI);
    scale_(s, 1 / PI);
    add_scalar_(s, -1.f);
    // Multiply by amplitude (in-place)
    scale_(s, amp);
}

void generate_signals(Matrix *s, int n_samples, fp range, bool noise) {
    // First component: sine wave
    Matrix s1 = extract_row(s, 0);
    generate_sine_wave(&s1, 1.5f, 0.3f, -PI, n_samples, range);
    // Second component: square wave
    Matrix s2 = extract_row(s, 1);
    generate_square_wave(&s2, 1, 0.5f, 0, n_samples, range);
    // Third component: sawtooth wave
    Matrix s3 = extract_row(s, 2);
    generate_sawtooth_wave(&s3, 0.5f, 0.7f, PI, n_samples, range);

    // Apply gaussian noise
    if (noise) {
        Matrix *ns = new_mat(s->height, s->width, false);

        mat_randn(ns);
        scale_(ns, 0.2f);
        add_mat_(s, ns);
        
        free_mat(ns);
    }
}
