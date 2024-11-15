"""
Makes the folder where the single-rotating model grid will be run.

```
python make_single_rot_grid.py single_grid template/single_rot
```

makes a folder `single_grid` with a substructure
`single_grid/mass/period/mass-ratio` containing the work directories.
"""

import sys
import makefolders as mf
import numpy as np
import os
import shutil
import getopt as go

n_m1_low = int((20 - 8) / 1)
n_m1_high = int((50 - 20) / 2)
n_m1_very_high = int((70 - 50) / 2.5) + 1
n_q = int((1 - 0.6) / 0.025)
n_p = int((np.log10(8) - np.log10(0.5)) / 0.04) + 1
m1s_low = np.round(np.linspace(8, 20, n_m1_low, endpoint=False), 4)
m1s_high = np.round(np.linspace(20, 50, n_m1_high, endpoint=False), 4)
m1s_very_high = np.round(np.linspace(50, 70, n_m1_very_high, endpoint=True), 4)
m1s = np.concatenate((m1s_low, m1s_high, m1s_very_high))
m1s = m1s[m1s <= 20]
# m1s = m1s[m1s <= 50]
qs = np.round(np.linspace(0.6, 1.0, n_q, endpoint=False), 3)
ps = np.round(np.geomspace(0.5, 8.0, n_p, endpoint=True), 4)
# ps = ps[ps > 1.00]
# ps = ps[ps < 3.0]
print('would make:')
print(len(m1s), 'm1s', m1s)
print(len(ps), 'ps', ps)
print(len(qs), 'qs', qs)
print('total:', len(m1s) * len(qs) * len(ps), 'runs')

if not input("continue? (Y/n) ") == "Y":
    exit()
opts, args = go.getopt(sys.argv[1:], 'o')
try:
    dest = args[0]
    template = args[1]
except IndexError:
    print('give destination + source template')
    exit()

if '-o' not in args:  # -o to overwrite without confirmation
    if os.path.exists(dest):
        if input('already exists! Deleting here what\'s present. Continue? (Y/n) ') == "Y":
            shutil.rmtree(dest)
        else:
            exit(0)

# print('making', len(m1s) * len(qs) * len(ps), 'cases')
count = 0
for m1 in m1s:
    for p in ps:
        for q in qs:
            mm1 = str(np.round(float(m1), 2))
            qq = str(np.round(float(q), 3))
            pp = str(np.round(float(p), 2))
            mf.make_single(str(m1), str(p), str(q), template, dest)
            count += 1
print("made", count, "cases")
