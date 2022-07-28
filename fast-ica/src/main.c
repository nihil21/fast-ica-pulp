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

#include <stdio.h>
#include <stdlib.h>
#include "pmsis.h"
#include "../include/random.h"
#include "../include/matrix.h"
#include "../include/stat_alloc.h"
#include "../include/signal.h"
#include "../include/fast_ica.h"

struct ClusterArgs {
    Matrix *x_l2;
    Matrix *res_l2;
};

// Cluster device
struct pi_device cluster_dev;

/*
 * Function executed by each core in the cluster
 */
static void cluster_fn(void *arg) {
    pi_perf_conf(
        1 << PI_PERF_INSTR |
        1 << PI_PERF_CYCLES
    );
    pi_perf_stop();
    pi_perf_reset();
    pi_perf_start();

    // Perform FastICA and measure performance
    fast_ica((FastICAArgs *) arg);

    pi_perf_stop();

    if (pi_core_id() == 0) {
        uint32_t instr_cnt = pi_perf_read(PI_PERF_INSTR);
        uint32_t cycles_cnt = pi_perf_read(PI_PERF_CYCLES);
        printf("Number of Instructions: %d\nClock Cycles: %d\nCPI: %f\n", instr_cnt, cycles_cnt, (float) cycles_cnt / instr_cnt);
    }
}

/*
 * Cluster entry-point
 */
static void cluster_entry(void *arg) {
    // Unpack arguments
    Matrix *x_l2 = ((struct ClusterArgs *) arg)->x_l2;
    Matrix *res_l2 = ((struct ClusterArgs *) arg)->res_l2;

    // Transfer observations to L1 using DMA
    Matrix x_l1 = {
        .data = in_data_l1,
        .height = OBS,
        .width = SAMPLES,
        .offset = SAMPLES,
        .on_cluster = true
    };
    mat_fc2cl(x_l2, &x_l1);

    // Prepare FastICA arguments
    FastICAStrategy strategy = STRATEGY == 0 ? Parallel : Deflation;
    GFunc g_function = G_FUNC == 0 ? LogCosh : (G_FUNC == 1 ? Exp : Cube);
    Matrix res_l1 = {
        .data = out_data_l1,
        .height = COMP,
        .width = SAMPLES,
        .offset = SAMPLES,
        .on_cluster = true
    };
    FastICAArgs fast_ica_arg = {&res_l1, &x_l1, COMP, true, strategy, g_function, 1e-7, MAX_ITER};

    // Spawn team of parallel processes
    pi_cl_team_fork(CORES, cluster_fn, (void *) &fast_ica_arg);

    // Transfer result back to L2 using DMA
    mat_cl2fc(&res_l1, res_l2);
}

/*
 * Run FastICA using the cluster
 */
static void run_fast_ica() {
    // PMSIS data structures
    struct pi_cluster_conf conf;
    struct pi_cluster_task cluster_task;

    // Set seed
    set_prng_seed(42);

    // Create matrix S of original signals (n_components, n_samples)
    Matrix *s = new_mat(COMP, SAMPLES, false);
    generate_signals(s, SAMPLES, S_LEN, NOISE);
    // Standardize signal
    Matrix *s_std = new_mat(s->height, 1, false);
    col_std(s_std, s);
    div_col_(s, s_std);

#if defined(VERB)
    print_mat(s);
#endif

    // Create mixing matrix A (n_observables, n_components)
    fp a_data[] = {1, 1, 1, 0.5f, 2, 1, 1.5f, 1, 2};  // , 0.3f, 0.9f, 1.7f};
    Matrix a = {.data = a_data, .height = OBS, .width = COMP, .offset = COMP, .on_cluster = false};

#if defined(VERB)
    print_mat(&a);
#endif

    // Create observations X by mixing signal S with matrix A (n_observables, n_samples)
    Matrix x_l2 = {
        .data = in_data_l2,
        .height = OBS,
        .width = SAMPLES,
        .offset = SAMPLES,
        .on_cluster = false
    };
    mat_mul(&x_l2, &a, s);

#if defined(VERB)
    print_mat(&x_l2);
#endif

    // Free memory
    free_mat(s);
    free_mat(s_std);

    // Prepare cluster task arguments
    Matrix res_l2 = {
        .data = out_data_l2,
        .height = COMP,
        .width = SAMPLES,
        .offset = SAMPLES,
        .on_cluster = false
    };
    struct ClusterArgs args = {&x_l2, &res_l2};

    // Prepare task
    pi_cluster_task(&cluster_task, cluster_entry, (void *) &args);

    // Initialize cluster device
    pi_cluster_conf_init(&conf);
    pi_open_from_conf(&cluster_dev, &conf);

    // Open cluster device
    if (pi_cluster_open(&cluster_dev))
        pmsis_exit(-1);

    // Execute task on cluster
    pi_cluster_send_task_to_cl(&cluster_dev, &cluster_task);
    
#if defined(VERB)
    print_mat(&res_l2);
#endif

    // Release cluster device
    pi_cluster_close(&cluster_dev);

    pmsis_exit(0);
}

int main() {
    return pmsis_kickoff((void *)run_fast_ica);
}
