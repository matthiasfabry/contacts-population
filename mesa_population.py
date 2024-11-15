"""
Defines an object containing a set of mesa history files on which to do calculations. Essentially a glorified array,
but useful nonetheless.
"""
from __future__ import annotations
import glob
import os
import sys
import pathlib as pl
import numpy as np
import mesa_data as md


def _get_index(val, values):
    if isinstance(val, float) or isinstance(val, int):
        ix = np.flatnonzero(values == float(val))
    elif isinstance(val, slice):
        if val.start is None and val.stop is None:
            ix = np.arange(0, len(values), 1, dtype=int)
        elif val.start is None:
            ix = np.flatnonzero(values < float(val.stop))
        elif val.stop is None:
            ix = np.flatnonzero(float(val.start) <= values)
        else:
            ix = np.flatnonzero(np.logical_and(float(val.start) <= values,
                                               values <= float(val.stop)))
    else:
        raise IndexError('give float or slice for indexing m1, q or p')
    return ix


class MesaPopulation:
    # top_dir/part/template/m/p/q/ET/ is the default structure where mesa work directories live
    def __init__(self, ddir, structure='*/*/*/*/*/*/', primonly=False, columns=None):

        if ddir[-1] != '/':
            ddir += '/'
        folders = glob.glob(ddir + structure)  # look in all folders downstream
        if len(folders) == 0:
            print('no mesa data found in this directory')
            return
        # determine m1s, qs, ps
        m1s = []
        qs = []
        ps = []
        for folder in folders:
            m, p, q = folder.rsplit('/')[-5:-2]
            m = float(m)
            q = float(q)
            p = float(p)
            if m not in m1s:
                m1s.append(m)
            if q not in qs:
                qs.append(q)
            if p not in ps:
                ps.append(p)
        self.m1s = np.sort(np.array(m1s))
        self.qs = np.sort(np.array(qs))
        self.ps = np.sort(np.array(ps))
        self.hists = np.empty((len(m1s), len(ps), len(qs), 3, 2), dtype=object)
        for i, folder in enumerate(folders):
            m, p, q, et = folder.rsplit('/')[-5:-1]
            print('loading folder {} of {}'.format(i, len(folders)), end='\r')
            if 'no_ET' in folder and os.path.exists(folder + 'same_as_ET.txt'):
                folder = str(pl.Path(folder).parent.joinpath('ET')) + '/'
                print('same as ET, looking there', folder)
            if primonly:
                try:
                    primhist = md.MesaData(glob.glob(folder + 'LOGS1/history.data')[0],
                                           read_data_cols=columns)
                except IndexError:
                    print('nothing found for', folder)
                    primhist = None
                except KeyError:
                    print(folder, 'does not contain a certain key')
                    primhist = None
                sechist = None
            else:
                try:
                    primhist, sechist = md.get_binary_hist(folder, read_data_cols=columns)
                except IndexError:
                    print('nothing found for', folder)
                    primhist, sechist = None, None
                except KeyError:
                    print(folder, 'does not contain a certain key')
                    primhist, sechist = None, None
            im = np.nonzero(self.m1s == float(m))[0]
            iq = np.nonzero(self.qs == float(q))[0]
            ip = np.nonzero(self.ps == float(p))[0]
            if et == 'ET':
                ie = 1
            elif et == 'single_rot':
                ie = 2
            else:
                ie = 0
            # print(im, iq, ip)
            self.hists[im, ip, iq, ie, 0] = primhist
            self.hists[im, ip, iq, ie, 1] = sechist

    def __getitem__(self, item):
        m, p, q, e, c = item
        im = _get_index(m, self.m1s)
        iq = _get_index(q, self.qs)
        ip = _get_index(p, self.ps)
        if not (len(im) > 0 and len(iq) > 0 and len(ip) > 0 and
                (e in ('ET', 'no_ET', 'single_rot', 'all')) and (c in ('prim', 'sec', 'both'))):
            print("nothing found for ", item, im, ip, iq)
            print(len(im), len(iq), len(ip),
                  (e in ('ET', 'no_ET', 'single_rot', 'all')), (c in ('prim', 'sec', 'both')))
            return None
        if e == 'ET':
            ie = [1]
        elif e == 'no_ET':
            ie = [0]
        elif e == 'single_rot':
            ie = [2]
        else:
            ie = np.array([0, 1, 2])
        if c == 'prim':
            ic = [0]
        elif c == 'sec':
            ic = [1]
        else:
            ic = np.array([0, 1])

        ret = self.hists[np.ix_(im, ip, iq, ie, ic)]
        if ret.size == 1:
            return ret[0, 0, 0, 0, 0]
        else:
            return ret


def merge_pops(this, other):
    for im in range(len(this.hists)):
        for ip in range(len(this.hists[im])):
            for iq in range(len(this.hists[im, ip])):
                for ie in range(len(this.hists[im, ip, iq])):
                    for ic in range(2):
                        if this.hists[im, ip, iq, ie, ic] is None and \
                                other.hists[im, ip, iq, ie, ic] is not None:
                            this.hists[im, ip, iq, ie, ic] = other.hists[im, ip, iq, ie, ic]


if __name__ == '__main__':
    pop = MesaPopulation(sys.argv[1])
    print(pop.ps)
