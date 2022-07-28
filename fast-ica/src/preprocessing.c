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

#include "../include/preprocessing.h"
#include "../include/stat_alloc.h"
#include "../include/linalg.h"
#include "../include/sorting.h"

/*
 * Array of eigenvalues indexes for sorting
 */
PI_L1 int sort_idx[OBS];

/*
 * Sort eigenvalues and eigenvectors
 */
void eigen_sort(const Matrix *eig_vals, const Matrix *eig_vecs) {
    int n_eig = eig_vals->height;
    int i, j, k;
    
    // Sort eigenvalues in descending order
    quick_sort(eig_vals->data, sort_idx, n_eig, true);

    // Make a copy of the original matrix
    Matrix eig_vecs_copy = {
        .data = prep_tmp2_data_l1,
        .height = OBS,
        .width = OBS,
        .offset = OBS,
        .on_cluster = true
    };
    copy_mat(&eig_vecs_copy, eig_vecs);

    // Reorder eig_vecs according to sort_idx
    for (i = 0; i < n_eig; i++) {
        int idx = sort_idx[i];
        Matrix eig_vec_copy = extract_col(&eig_vecs_copy, idx);
        Matrix eig_vec = extract_col(eig_vecs, i);
        copy_mat(&eig_vec, &eig_vec_copy);
    }
}

/*
 * Perform whitening of input data
 */
void whitening(const Matrix *x_w, const Matrix *white_mat, const Matrix *mean_vec, const Matrix *x, int n_components) {
    bool over_det = n_components < x->height;

    // 1. Compute covariance matrix
    Matrix cov_mat = {
        .data = prep_mat_data_l1,
        .height = OBS,
        .width = OBS,
        .offset = OBS,
        .on_cluster = true
    };
    covariance(&cov_mat, mean_vec, x);

    // 2. Compute eigenvalues and eigenvectors
    Matrix eig_vals = {
        .data = eig_vals_data_l1,
        .height = OBS,
        .width = 1,  // column vector
        .offset = 1,
        .on_cluster = true
    };
    Matrix eig_vecs = {
        .data = eig_vecs_data_l1,
        .height = OBS,
        .width = OBS,
        .offset = OBS,
        .on_cluster = true
    };
    solve_eig(&eig_vals, &eig_vecs, &cov_mat);
    
    // Sort eigenvalues and eigenvectors, if required
    if (over_det)
        eigen_sort(&eig_vals, &eig_vecs);

    // Create diagonal matrix of eigenvalues
    Matrix d = {
        .data = prep_tmp2_data_l1,
        .height = n_components,
        .width = n_components,
        .offset = n_components,
        .on_cluster = true
    };
    eye(&d);
    for (int i = 0; i < n_components; i++)
        MAT_CELL(&d, i, i) = 1 / SQRT(MAT_CELL(&eig_vals, i, 0));

    // 3. Compute whitening matrix
    if (over_det)
        eig_vecs = slice(&eig_vecs, 0, eig_vecs.height - 1, 0, n_components - 1);  // keep only the first n_components columns
    // print_mat(&eig_vecs);
    mat_mul_t2(white_mat, &d, &eig_vecs);

    // 4. Whiten data
    mat_mul(x_w, white_mat, x);
}
