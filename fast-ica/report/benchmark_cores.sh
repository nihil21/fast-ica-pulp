#!/bin/bash

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

# Set working directory
cd ../

# Iterate over strategies
for c in {1..8}; do
    echo "Cores $c" | tee -a report/report_cores.txt;
    make clean all CORES=$c OBS=3 DEBUG=1
    make run | tee -a report/report_cores.txt
done
