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

#include "../include/fast_ica.h"
#include "../include/stat_alloc.h"
#include "../include/linalg.h"
#include "../include/preprocessing.h"

/*
 * LogCosh function
 */
void logcosh_f(const Matrix *gx, const Matrix *gx_prime, const Matrix *x) {
    fp alpha = 1.f;

    for (int i = 0; i < x->height; i++) {
        for (int j = 0; j < x->width; j++) {
            MAT_CELL(gx, i, j) = TANH(alpha * MAT_CELL(x, i, j));
            MAT_CELL(gx_prime, i, j) = alpha * (1 - MAT_CELL(gx, i, j) * MAT_CELL(gx, i, j));
        }
    }
}

/*
 * Exponential function
 */
void exp_f(const Matrix *gx, const Matrix *gx_prime, const Matrix *x) {
    for (int i = 0; i < x->height; i++) {
        for (int j = 0; j < x->width; j++) {
            fp tmp = EXP(-MAT_CELL(x, i, j) * MAT_CELL(x, i, j) / 2);
            MAT_CELL(gx, i, j) = MAT_CELL(x, i, j) * tmp;
            MAT_CELL(gx_prime, i, j) = (1 - MAT_CELL(x, i, j) * MAT_CELL(x, i, j)) * tmp;
        }
    }
}

/*
 * Cube function
 */
void cube_f(const Matrix *gx, const Matrix *gx_prime, const Matrix *x) {
    for (int i = 0; i < x->height; i++) {
        for (int j = 0; j < x->width; j++) {
            fp tmp = MAT_CELL(x, i, j) * MAT_CELL(x, i, j);
            MAT_CELL(gx, i, j) = MAT_CELL(x, i, j) * tmp;
            MAT_CELL(gx_prime, i, j) = 3 * tmp;
        }
    }
}

/*
 * Implement Gram-Schmidt decorrelation
 */
void gram_schmidt_decorrelation(const Matrix *w_i_new, const Matrix *w, const int i) {
    if (i > 0) {
        Matrix tmp_vec = {
            .data = prep_tmp1_data_l1,
            .height = 1,  // row vector
            .width = w->width,
            .offset = w->width,
            .on_cluster = true
        };
        Matrix tmp_mat = {
            .data = prep_tmp3_data_l1,
            .height = w->width,
            .width = w->width,
            .offset = w->width,
            .on_cluster = true
        };
        Matrix w_slice = slice(w, 0, i - 1, 0, w->width - 1);
        mat_mul_t1(&tmp_mat, &w_slice, &w_slice);
        mat_mul(&tmp_vec, w_i_new, &tmp_mat);
        sub_mat_(w_i_new, &tmp_vec);
    }
}

/*
 * Implement symmetric decorrelation
 *
void symmetric_decorrelation(Matrix **w) {
    Matrix *w_wt = mat_mul_t2(*w, *w);
    // Compute eigenvalues and eigenvectors
    Tuple *eigen = solve_eig(w_wt);
    free_mat(w_wt);
    Matrix *eig_vals = eigen->m1;  // column vector
    Matrix *eig_vecs = eigen->m2;
    int n = eig_vals->height;
    Matrix *d = new_mat(n, n, (*w)->on_cluster);
    for (int i = 0; i < n; i++)
        MAT_CELL(d, i, i) = 1 / SQRT(MAT_CELL(eig_vals, i, 0));
    // Compute new weight matrix
    Matrix *tmp1 = mat_mul_t1(eig_vecs, *w);
    Matrix *tmp2 = mat_mul(d, tmp1);
    free_mat(d);
    free_mat(tmp1);

    free_mat(*w);
    *w = mat_mul(eig_vecs, tmp2);
    free_tuple(eigen, true);
    free_mat(tmp2);
    
}

/*
 * Implement FastICA deflationary strategy
 */
Matrix *ica_def(const Matrix *w, const Matrix *x_w, void (*g_func)(const Matrix *, const Matrix *, const Matrix *), fp threshold, int max_iter) {
    int n_units = x_w->height;
    int n_samples = x_w->width;

    // Initialize weights randomly
    mat_randn(w);
    // Iterate over units
    for (int k = 0; k < n_units; k++) {
        // Initialize i-th neuron
        Matrix w_k = extract_row(w, k);  // row vector
        scale_(&w_k, 1 / norm(&w_k));

        for (int i = 0; i < max_iter; i++) {
            // (1, n_units) @ (n_units, n_samples) -> (1, n_samples)
            Matrix ws = {
                .data = ica_tmp_data_l1,
                .height = 1,
                .width = n_samples,
                .offset = n_samples,
                .on_cluster = true
            };
            mat_mul(&ws, &w_k, x_w);
            // Compute G_Ws and G_Ws'
            Matrix g_ws = {
                .data = ica_g_data_l1,
                .height = 1,  // row vector
                .width = n_samples,
                .offset = n_samples,
                .on_cluster = true
            };
            Matrix g_ws_prime = {
                .data = ica_g_prime_data_l1,
                .height = 1,  // row vector
                .width = n_samples,
                .offset = n_samples,
                .on_cluster = true
            };
            g_func(&g_ws, &g_ws_prime, &ws);

            // (1, n_samples) @ (n_units, n_samples).T -> (1, n_samples) @ (n_samples, n_units) -> (1, n_units)
            Matrix a = {
                .data = prep_tmp1_data_l1,
                .height = 1,
                .width = n_units,
                .offset = n_units,
                .on_cluster = true
            };
            mat_mul_t2(&a, &g_ws, x_w);
            scale_(&a, 1 / (fp) n_samples);
            // (1, n_units) * E[(1, n_samples)] -> (1, n_units)
            Matrix b = {
                .data = prep_tmp3_data_l1,
                .height = 1,
                .width = n_units,
                .offset = n_units,
                .on_cluster = true
            };
            scale(&b, &w_k, mean(&g_ws_prime));

            // Compute new weight
            Matrix w_k_new = {
                .data = prep_tmp4_data_l1,
                .height = 1,
                .width = n_units,
                .offset = n_units,
                .on_cluster = true
            };
            sub_mat(&w_k_new, &a, &b);  // (1, n_units)
            // Decorrelate
            gram_schmidt_decorrelation(&w_k_new, w, k);
            // Normalize
            scale_(&w_k_new, 1 / norm(&w_k_new));

            // Compute distance
            fp distance = ABS(dot(&w_k_new, &w_k) - 1.f);

            // Update weight
            copy_mat(&w_k, &w_k_new);

            if (distance < threshold) {
                break;
            }
        }
    }
}

/*
 * Implement FastICA parallel strategy
 *
Matrix *ica_par(Matrix *x_w, void (*g_func)(const Matrix *, const Matrix *, const Matrix *), fp threshold, int max_iter) {
    bool on_cluster = x_w->on_cluster;
    int n_units = x_w->height;
    int n_samples = x_w->width;

    // Initialize weights randomly and decorrelate
    Matrix *w = mat_randn(n_units, n_units, on_cluster);
    symmetric_decorrelation(&w);

    for (int i = 0; i < max_iter; i++) {
        // (n_units, n_units) @ (n_units, n_samples) -> (n_units, n_samples)
        Matrix *ws = mat_mul(w, x_w);
        // Compute G_Ws and G_Ws'
        Tuple *res = g_func(ws);
        free_mat(ws);
        Matrix *g_ws = res->m1;  // (n_units, n_samples)
        Matrix *g_ws_prime = res->m2;  // (n_units, n_samples)

        Matrix *a = new_mat(n_units, n_units, on_cluster);
        Matrix *b = new_mat(n_units, n_units, on_cluster);
        // Iterate over units
        for (int k = 0; k < n_units; k++) {
            // Extract k-th row from G_Ws
            Matrix *g_ws_k = extract_row(g_ws, k);  // row vector
            // (1, n_samples) @ (n_units, n_samples).T -> (1, n_samples) @ (n_samples, n_units) -> (1, n_units)
            Matrix *a_k = mat_mul_t2(g_ws_k, x_w);
            scale_(a_k, 1 / (fp) n_samples);
            free_mat(g_ws_k);

            // Extract k-th row from G_Ws'
            Matrix *g_ws_prime_k = extract_row(g_ws_prime, k);  // row vector
            // Extract k-th row from W
            Matrix *w_k = extract_row(w, k);  // row vector
            // (1, n_units) * E[(1, n_samples)] -> (1, n_units)
            Matrix *b_k = scale(w_k, mean(g_ws_prime_k));
            free_mat(g_ws_prime_k);
            free_mat(w_k);

            // Paste rows
            paste_row(a, a_k, k);
            paste_row(b, b_k, k);
            free_mat(a_k);
            free_mat(b_k);
        }
        free_tuple(res, true);

        // Compute new weight
        Matrix *w_new = sub_mat(a, b);
        free_mat(a);
        free_mat(b);
        // Decorrelate
        symmetric_decorrelation(&w_new);

        // Compute distance
        fp distance = 0;
        Matrix *tmp1 = mat_mul_t2(w_new, w);
        Matrix *tmp2 = diagonal(tmp1);
        free_mat(tmp1);
        for (int ii = 0; ii < tmp2->height; ii++) {
            fp cur_dis = ABS(ABS(MAT_CELL(tmp2, ii, 0)) - 1);
            if (cur_dis > distance)
                distance = cur_dis;
        }
        free_mat(tmp2);

        // Update weights
        free_mat(w);
        w = w_new;

        if (distance < threshold) {
            printf("N. iterations: %d/%d\n", i, max_iter);
            break;
        }
    }

    return w;
}

/*
 * Perform FastICA on a (n_features, n_samples) matrix of observations
 */
void fast_ica(FastICAArgs *arg) {
    // Read input arguments
    Matrix *res = arg->res;
    Matrix *x = arg->x;
    int n_components = arg->n_components;
    bool whiten = arg->whiten;

    int n_features = x->height;
    int n_samples = x->width;

    // Prepare matrices
    Matrix x_w = {
        .data = in_w_data_l1,
        .height = n_components,
        .width = n_samples,
        .offset = n_samples,
        .on_cluster = true
    };
    Matrix white_mat = {
        .data = prep_mat_data_l1,
        .height = n_components,
        .width = n_features,
        .offset = n_features,
        .on_cluster = true
    };
    Matrix mean_vec = {
        .data = prep_vec_data_l1,
        .height = n_features,
        .width = 1,
        .offset = 1,
        .on_cluster = true
    };

    // Whiten data, if specified
    if (whiten) {
        // Avoid under-determined systems
        int max_comp = n_features > n_samples ? n_features : n_samples;
        if (n_components > max_comp)
            n_components = max_comp;

        whitening(&x_w, &white_mat, &mean_vec, x, n_components);
    } else {
        copy_mat(&x_w, x);
    }

    // Select non-quadratic function G
    void (*g)(const Matrix *, const Matrix *, const Matrix *);
    switch (arg->g_func) {
        case LogCosh:
            g = &logcosh_f;
            break;
        case Exp:
            g = &exp_f;
            break;
        case Cube:
            g = &cube_f;
            break;
    }

    Matrix w = {
        .data = prep_tmp2_data_l1,
        .height = COMP,
        .width = COMP,
        .offset = COMP,
        .on_cluster = true
    };
    // Select strategy to estimate weight matrix
    switch (arg->strategy) {
        case Parallel:
            // Execute parallel algorithm
            ica_def(&w, &x_w, g, arg->threshold, arg->max_iter);
            break;
        case Deflation:
            // Execute deflation algorithm
            ica_def(&w, &x_w, g, arg->threshold, arg->max_iter);
            break;
    }

    // Reconstruct signal
    mat_mul(res, &w, &x_w);
}
