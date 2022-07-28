# Copyright 2022 Mattia Orlandi
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


APP = fast-ica

APP_SRCS += src/main.c src/fp.c src/utils.c src/random.c src/sorting.c src/matrix.c src/stat_alloc.c src/signal.c src/linalg.c src/preprocessing.c src/fast_ica.c src/cluster.c

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

# Lenght of the recording, sampling rate and noise
S_LEN ?= 10
S_RATE ?= 100
NOISE ?= 1
APP_CFLAGS += -D S_LEN=$(S_LEN) -D S_RATE=$(S_RATE) -D NOISE=$(NOISE)

# Debug and verbosity
ifdef DEBUG
APP_CFLAGS += -D DEBUG
endif
ifdef VERB
APP_CFLAGS += -D VERB
endif

# Optimization flags
APP_CFLAGS += -g -Os

# Math library linking
APP_LDFLAGS += -lm

include $(RULES_DIR)/pmsis_rules.mk
