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


import sys
import numpy as np
from matplotlib import pyplot as plt


def plot_signals(s, x, s_) -> None:
    plt.figure(figsize=(15, 10))

    models = [x, s, s_]
    names = [
        "Observations (mixed signal)",
        "True Sources",
        "ICA recovered signals",
    ]
    plots = len(models)

    for ii, (model, name) in enumerate(zip(models, names)):
        plt.subplot(plots, 1, ii + 1)
        plt.title(name)
        for sig in model:
            plt.plot(sig)
        plt.grid()

    plt.tight_layout()
    plt.show()


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: plot_signals.py REPORT.LOG N_COMP N_REC")
    path = sys.argv[1]
    n_comp = int(sys.argv[2])
    n_rec = int(sys.argv[3])

    s = []
    a = []
    x = []
    s_ = []
    with open(path, "r") as f:
        lines = f.readlines()    
        # Read input signal S
        for line in lines[:n_comp]:
            s.append(eval("["+ line.strip().replace(" ", ", ") + "]"))
        s = np.array(s)
        # Read mixing matrix A
        for line in lines[n_comp:n_comp + n_rec]:
            a.append(eval("["+ line.strip().replace(" ", ", ") + "]"))
        a = np.array(a)
        # Read observations X
        for line in lines[n_comp + n_rec:n_comp + 2 * n_rec]:
            x.append(eval("["+ line.strip().replace(" ", ", ") + "]"))
        x = np.array(x)
        # Read restored signal S_
        for line in lines[n_comp + 2 * n_rec:2 * (n_comp + n_rec)]:
            s_.append(eval("["+ line.strip().replace(" ", ", ") + "]"))
        s_ = np.array(s_)
    print(a)
    plot_signals(s, x, s_)


if __name__ == "__main__":
    main()
