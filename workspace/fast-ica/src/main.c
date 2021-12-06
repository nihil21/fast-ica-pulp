#include <stdio.h>
#include <stdlib.h>
#include "pmsis.h"
#include "../include/random.h"
#include "../include/matrix.h"
#include "../include/signal.h"
#include "../include/cluster.h"
#include "../include/fast_ica.h"
#include "../include/preprocessing.h"

// Input/output matrices (on L2 and L1)
PI_L2 static fp x_data[OBS * SAMPLES];
PI_L2 static Matrix x = {.data = x_data, .height = OBS, .width = SAMPLES};
PI_L2 static fp res_data[COMP * SAMPLES];
PI_L2 static Matrix res = {.data = res_data, .height = COMP, .width = SAMPLES};

// Cluster device
struct pi_device cluster_dev;

/*
 * Function to measure the performance of a given function
 */
void perf(void (*func)(FastICAArgs *), FastICAArgs *arg) {
    pi_perf_conf(
        1 << PI_PERF_INSTR |
        1 << PI_PERF_CYCLES
    );
    pi_perf_stop();
    pi_perf_reset();
    pi_perf_start();

    func(arg);

    pi_perf_stop();

    uint32_t instr_cnt = pi_perf_read(PI_PERF_INSTR);
    uint32_t cycles_cnt = pi_perf_read(PI_PERF_CYCLES);
    printf("Number of Instructions: %d\nClock Cycles: %d\nCPI: %f\n", instr_cnt, cycles_cnt, (float) cycles_cnt / instr_cnt);
}

/*
 * Cluster entry-point
 */
static void cluster_entry(void *arg) {
    // Transfer observations to L1 using DMA
    Matrix *x_cl = new_mat(OBS, SAMPLES, true);
    mat_fc2cl(&x, x_cl);

    // Prepare FastICA arguments
    FastICAStrategy strategy = STRATEGY == 0 ? Parallel : Deflation;
    GFunc g_function = G_FUNC == 0 ? LogCosh : (G_FUNC == 1 ? Exp : Cube);
    Matrix *res_cl;
    FastICAArgs fast_ica_arg = {&res_cl, x_cl, COMP, true, strategy, g_function, 1e-4, MAX_ITER};

    // Perform FastICA and measure performance
    perf(fast_ica, &fast_ica_arg);

    // Transfer result back to L2
    mat_cl2fc(res_cl, &res);

    // Free memory
    free_mat(x_cl);
    free_mat(res_cl);
}

static void test_fast_ica() {
    // PMSIS data structures
    struct pi_cluster_conf conf;
    struct pi_cluster_task cluster_task = {0};

    // Set seed
    set_prng_seed(42);

    // Create matrix S of original signals (n_components, n_samples)
    Matrix *s = generate_signals(SAMPLES, WINDOW_SIZE, ADD_NOISE);
    // Standardize signal
    Matrix *s_std = col_std(s);
    div_col_(s, s_std);
    if (VERB)
        print_mat(s);

    // Create mixing matrix A (n_observables, n_components)
    fp a_data[] = {1, 1, 1, 0.5f, 2, 1, 1.5f, 1, 2, 0.3f, 0.9f, 1.7f};
    Matrix *a = from_array(a_data, OBS, COMP, false);
    if (VERB)
        print_mat(a);

    // Create observation X by mixing signal S with matrix A (n_observables, n_samples)
    for (int i = 0; i < a->height; i++) {
        for (int k = 0; k < a->width; k++) {
            for (int j = 0; j < s->width; j++) {
                x.data[i * x.width + j] += MAT_CELL(a, i, k) * MAT_CELL(s, k, j);
            }
        }
    }
    if (VERB)
        print_mat(&x);

    // Prepare task
    pi_cluster_task(&cluster_task, cluster_entry, NULL);

    // Initialize cluster device
    pi_cluster_conf_init(&conf);
    pi_open_from_conf(&cluster_dev, &conf);

    // Open cluster device
    if (pi_cluster_open(&cluster_dev))
        pmsis_exit(-1);

    // Execute task on cluster
    pi_cluster_send_task_to_cl(&cluster_dev, &cluster_task);
    
    if (VERB)
        print_mat(&res);
    
    // Release cluster device
    pi_cluster_close(&cluster_dev);

    // Free memory
    free_mat(s);
    free_mat(s_std);
    free_mat(a);

    pmsis_exit(0);
}

int main() {
    return pmsis_kickoff((void *)test_fast_ica);
}