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

#ifndef FAST_ICA_STAT_ALLOC_H
#define FAST_ICA_STAT_ALLOC_H

#include "matrix.h"

#define SAMPLES S_LEN * S_RATE

/*
 * Static array for input data (on L2 and L1)
 */
PI_L2 extern fp in_data_l2[OBS * SAMPLES];
PI_L1 extern fp in_data_l1[OBS * SAMPLES];
PI_L1 extern fp in_w_data_l1[COMP * SAMPLES];  // whitening

/*
 * Static array for output data (on L2 and L1)
 */
PI_L2 extern fp out_data_l2[COMP * SAMPLES];
PI_L1 extern fp out_data_l1[COMP * SAMPLES];

/*
 * Static array for preprocessing step (on L1)
 */
PI_L1 extern fp prep_vec_data_l1[MAX(COMP, OBS)];
PI_L1 extern fp prep_mat_data_l1[MAX(COMP, OBS) * MAX(COMP, OBS)];
PI_L1 extern fp prep_tmp1_data_l1[MAX(COMP, OBS)];
PI_L1 extern fp prep_tmp2_data_l1[MAX(COMP, OBS) * MAX(COMP, OBS)];
PI_L1 extern fp prep_tmp3_data_l1[MAX(COMP, OBS) * MAX(COMP, OBS)];
PI_L1 extern fp prep_tmp4_data_l1[MAX(COMP, OBS) * MAX(COMP, OBS)];
PI_L1 extern fp prep_tmp5_data_l1[MAX(COMP, OBS) * MAX(COMP, OBS)];
PI_L1 extern fp prep_tmp6_data_l1[MAX(COMP, OBS) * MAX(COMP, OBS)];
PI_L1 extern fp prep_tmp7_data_l1[MAX(COMP, OBS) * MAX(COMP, OBS)];

/*
 * Static array for eigenvalues and eigenvectors (on L1)
 */
PI_L1 extern fp eig_vals_data_l1[MAX(COMP, OBS)];
PI_L1 extern fp eig_vecs_data_l1[MAX(COMP, OBS) * MAX(COMP, OBS)];

/*
 * Static array for FastICA (on L1)
 */
PI_L1 extern fp ica_g_data_l1[SAMPLES];
PI_L1 extern fp ica_g_prime_data_l1[SAMPLES];
PI_L1 extern fp ica_tmp_data_l1[SAMPLES];

#endif //FAST_ICA_STAT_ALLOC_H
