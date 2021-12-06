APP = fast-ica

APP_SRCS += src/main.c src/utils.c src/random.c src/sorting.c src/tuple.c src/signal.c src/matrix.c src/linalg.c src/preprocessing.c src/fast_ica.c src/cluster.c

# Number of components and observables
COMP ?= 3
OBS ?= 4
APP_CFLAGS += -D COMP=$(COMP) -D OBS=$(OBS)

# Number of cores of the cluster to use
CORES ?= 1
APP_CFLAGS += -D CORES=$(CORES)

# Strategy, G function and maximum number of iterations
STRATEGY ?= 0
G_FUNC ?= 0
MAX_ITER ?= 3000
APP_CFLAGS += -D STRATEGY=$(STRATEGY) -D G_FUNC=$(G_FUNC) -D MAX_ITER=$(MAX_ITER)

# Number of samples, sampling window size and noise
SAMPLES ?= 1000
WINDOW_SIZE ?= 10
ADD_NOISE ?= 1
APP_CFLAGS += -D SAMPLES=$(SAMPLES) -D WINDOW_SIZE=$(WINDOW_SIZE) -D ADD_NOISE=$(ADD_NOISE)

# Verbosity
VERB ?= 0
APP_CFLAGS += -D VERB=$(VERB)

# Optimization flags
APP_CFLAGS += -g -Os

# Math library linking
APP_LDFLAGS += -lm

include $(RULES_DIR)/pmsis_rules.mk
