"""
Loads a single mesa history or profile file, along with wrappers to easily plot certain properties.
Copied and adapted from the mkipp package. See https://github.com/orlox/mkipp.
"""
import glob
import re
from collections import OrderedDict

import h5py
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma


def _read_file(file, is_hdf5, columns, column_names):
    if not is_hdf5:
        # read data
        data = np.loadtxt(file, skiprows=6,
                          usecols=tuple([columns[k] for k in column_names]),
                          unpack=True)
    else:
        file = h5py.File(file, "r")
        data = file['data_vals'][:, sorted([columns[k] for k in column_names])]
        data = data.transpose()
        file.close()
    return data


# class to extract a mesa data file.
class MesaData:
    def __init__(self, history_file, only_read_header=False, read_data=True,
                 read_data_cols=None,
                 clean_data=True, is_hdf5=False):
        if read_data_cols is None:
            read_data_cols = []
        self.filename = history_file
        self.is_hdf5 = is_hdf5
        # header is a dictionary with the general info from the second and
        # third line of file
        self.header = {}
        self.header_num = {}
        self.data = None
        self.num_points = 0
        self.read_columns = None
        # columns is a dictionary which gives the column number (minus 1)
        # corresponding to the key
        self.columns = {}
        columns = []
        if not is_hdf5:
            file = open(self.filename, "r")
            # first line is not used
            file.readline()
            # following two lines have header data
            header_names = file.readline().split()
            header_vals = file.readline().split()
            for i, header_name in enumerate(header_names):
                # need to properly account for these new columns
                if header_name in ["version_number", "compiler", "build",
                                   "MESA_SDK_version", "date", "math_backend"]:
                    continue
                self.header[header_name] = float(header_vals[i])
                self.header_num[header_name] = i
                i += 1
            if only_read_header:
                file.close()
                return
            # next line is empty
            file.readline()
            # following two lines have column data
            nums = file.readline().split()
            names = file.readline().split()
            for i, name in enumerate(names):
                self.columns[name] = int(nums[i]) - 1
                columns.append(name)
            file.close()
        else:
            file = h5py.File(self.filename, "r")
            header_names = file['header_names'][:]
            header_vals = file['header_vals'][:]
            for i in range(len(header_names)):
                key = header_names[i].decode('utf-8')
                self.header[key] = header_vals[i]
                self.header_num[key] = i
            columns = file['data_names'][:].tolist()
            for i, col in enumerate(columns):
                self.columns[col.decode('utf-8')] = i
                columns[i] = col.decode('utf-8')
            file.close()

        if not read_data:
            return

        if len(read_data_cols) == 0:
            read_data_cols = columns
        self.read_data(read_data_cols, clean_data=clean_data)
        # self.revise_stop_condition()

    def revise_stop_condition(self):
        if 'stopping_condition' in self.data.keys():
            with open(self.filename.rsplit('/', maxsplit=2)[0] + '/term.txt', 'r') as out:
                lns = out.readlines()
                for line in lns:
                    if line.find('stopping because of problems') != -1:
                        self.data['stopping_condition'][-1] = -99
                        return

    def __repr__(self):
        return self.filename

    def read_data(self, column_names, clean_data=True, sample_every_n=1):
        # always include model_number if its part of the data
        if "model_number" not in column_names and "model_number" in \
                self.columns:
            column_names.append("model_number")

        # be sure there are no repeated column names
        # (could use set but that breaks the order of the columns, which is
        # needed if I want to save the file)
        column_names = list(OrderedDict.fromkeys(column_names))
        self.read_columns = column_names

        data = _read_file(self.filename, self.is_hdf5, self.columns, column_names)

        self.data = {}
        # Be careful in case only one column is required
        if len(column_names) > 1:
            for i, column in enumerate(column_names):
                self.data[column] = np.atleast_1d(data[i])
        else:
            self.data[column_names[0]] = data

        # clean redos
        if clean_data and "model_number" in self.columns and len(
                self.data["model_number"]) > 1:
            # create a mask
            model_number = self.data["model_number"]
            mask = np.zeros(len(model_number))
            max_model_number = model_number[-1]
            # last entry is valid, start from there and remove repeats
            for i in range(len(model_number) - 2, -1, -1):
                if model_number[i] >= max_model_number:
                    # exclude this point
                    mask[i] = 1
                else:
                    max_model_number = model_number[i]

            if sum(mask) > 0:
                for column in column_names:
                    self.data[column] = ma.masked_array(self.data[column],
                                                        mask=mask).compressed()

        # subsample points
        if sample_every_n > 1 and "model_number" in self.columns and len(
                self.data["model_number"]) > 2:
            # keep first and last entry
            # create a mask
            model_number = self.data["model_number"]
            mask = np.zeros(len(model_number))
            for i in range(1, len(model_number) - 1):
                if (i + 1) % sample_every_n != 0:
                    # exclude this point
                    mask[i] = 1

            if sum(mask) > 0:
                for column in column_names:
                    self.data[column] = ma.masked_array(self.data[column],
                                                        mask=mask).compressed()

        # count number of points using first entry in dict
        self.num_points = len(self.data[self.read_columns[0]])

    def append_other_hist(self, other_history_file, only_read_header=False,
                          read_data=True, read_data_cols=None, clean_data=True,
                          is_hdf5=False):
        # read other data
        other = MesaData(other_history_file, only_read_header, read_data,
                         read_data_cols, clean_data, is_hdf5)
        # merge data
        for i in range(len(self.read_columns)):
            if self.read_columns[i] in other.read_columns:
                # found corresponding cols, merge them
                self.data[self.read_columns[i]] = np.concatenate((
                    self.data[self.read_columns[i]],
                    other.data[self.read_columns[i]]))
            else:
                # not found, add nans for that column
                self.data[self.read_columns[i]] = np.concatenate(
                    (self.data[self.read_columns[i]],
                     np.full(len(other.data[other.read_columns[0]]), np.nan)))

        # go over columns in other not in self
        for j in range(len(other.read_columns)):
            if other.read_columns[j] not in self.read_columns:
                self.data[other.read_columns[j]] = np.concatenate((
                    np.full(self.num_points, np.nan),
                    other.data[other.read_columns[j]]))
                self.read_columns.append(other.read_columns[j])

        self.num_points = len(self.data[self.read_columns[0]])

    def get(self, key):
        return self.data[key]

    def save_as_hdf5(self, filename, header_str_dtype="S28",
                     data_str_dtype="S40", compression_opts=4):
        f = h5py.File(filename, "w")
        dset_header_names = f.create_dataset("header_names",
                                             (len(self.header),),
                                             dtype=header_str_dtype)
        dset_header_vals = f.create_dataset("header_vals", (len(self.header),),
                                            dtype="d")
        for key in self.header:
            dset_header_names[self.header_num[key]] = np.string_(key)
            dset_header_vals[self.header_num[key]] = self.header[key]
        dset_column_names = f.create_dataset("data_names",
                                             (len(self.read_columns),),
                                             dtype=data_str_dtype)
        dset_column_vals = f.create_dataset("data_vals", (
            self.num_points, len(self.read_columns)), dtype="d",
                                            compression='gzip',
                                            compression_opts=compression_opts)
        for k, key in enumerate(self.read_columns):
            dset_column_names[k] = np.string_(key)
            dset_column_vals[:, k] = self.data[key]
        f.close()

    # creates a mesa look-alike output file
    # prints all integers as doubles
    # not the most efficient code but I don't care
    def save_as_ascii(self, filename, header_str_format="{0:>28}",
                      header_double_format="{0:>28.16e}",
                      data_str_format="{0:>40}",
                      data_double_format="{0:>40.16e}"):
        f = open(filename, "w")
        for i in range(len(list(self.header))):
            f.write(header_str_format.format(i + 1))
        f.write("\n")
        # create an ordered list of keys
        header_keys = []
        for i in range(len(list(self.header))):
            for key in self.header:
                if self.header_num[key] == i:
                    header_keys.append(key)
                    break
        for i, key in enumerate(header_keys):
            f.write(header_str_format.format(key))
        f.write("\n")
        for i, key in enumerate(header_keys):
            f.write(header_double_format.format(self.header[key]))
        f.write("\n")
        f.write("\n")

        for i in range(len(list(self.read_columns))):
            f.write(data_str_format.format(i + 1))
        f.write("\n")

        for i, key in enumerate(self.read_columns):
            f.write(data_str_format.format(key))
        for k in range(self.num_points):
            f.write("\n")
            for i, key in enumerate(self.read_columns):
                f.write(data_double_format.format(self.data[key][k]))

        f.close()

    def trim(self, param, minpar, maxpar):
        try:
            vals = self.get(param)
        except KeyError:
            print('cannot find {} in mesadata'.format(param))
            return
        mask = np.empty(len(vals), dtype=bool)
        for i in range(len(vals)):
            mask[i] = minpar <= vals[i] <= maxpar
        for key in self.data:
            self.data[key] = self.data[key][mask]

    def trim_PMS(self, step_req: int = 1, age_in: float = None):
        print('trimming', repr(self))
        if age_in is None:
            try:
                log_l = self.get('log_L')
                log_lnuc = self.get('log_Lnuc')
                age = self.get('star_age')
            except KeyError:
                print(repr(self), 'cannot obtain PMS mask')
                return
            zamsix = 0
            count = 0
            while count < step_req:
                zamsix += 1
                if abs(log_l[zamsix] - log_lnuc[zamsix]) > 0.005:
                    count = 0
                else:
                    count += 1
            # print('zams found at model no: ', self.get('model_number')[zamsix])
            # print(logL[zamsix], logLnuc[zamsix], age[zamsix])
        else:
            try:
                age = self.get('star_age')
            except KeyError:
                print('cannot obtain PMS mask')
                return
            zamsix = np.argwhere(
                age_in - age == np.min(abs(age - age_in)))[0, 0]
            # print('trimming from model no', self.get('model_number')[zamsix])

        for key in self.data:
            self.data[key] = self.data[key][zamsix:]
        return age[zamsix], zamsix


# reads the profiles.index files in the folders specified by the logs_dirs
# array and returns
# an array containing paths to the individual profile files, after cleaning
# up redos and backups
def get_profile_paths(logs_dirs=None):
    if logs_dirs is None:
        logs_dirs = ["LOGS/"]
    elif isinstance(logs_dirs, str):
        if logs_dirs[-1] != '/':
            logs_dirs += '/'
        logs_dirs = [logs_dirs]
    profile_paths = []
    for log_dir in logs_dirs:
        model_number, paths = np.loadtxt(log_dir + "profiles.index",
                                         skiprows=1, usecols=(0, 2),
                                         unpack=True)
        mask = np.zeros(len(paths))
        max_model_number = model_number[-1]
        # last entry is valid, start from there and remove repeats
        for i in range(len(model_number) - 2, -1, -1):
            if model_number[i] >= max_model_number:
                mask[i] = 1
            else:
                max_model_number = model_number[i]

        if sum(mask) > 0:
            paths = ma.masked_array(paths, mask=mask).compressed()
        profile_paths.extend(
            [log_dir + "profile" + str(int(i)) + ".data" for i in paths])
    return profile_paths


def mesa_prof_data(lst):
    ret = list()
    for i in lst:
        ret.append(MesaData(i))
    return ret


def plot(mdata, x, y, *args, mask=None, xscale=1., yscale=1., labels=True, **kwargs):
    xdata = mdata.get(x) / xscale
    if mask is None:
        ydata = mdata.get(y)
    else:
        ydata = np.ma.masked_where(mask, mdata.get(y))

    plt.plot(xdata, ydata / yscale, *args, **kwargs)
    if labels:
        plt.xlabel(r'{}'.format(x))
        plt.ylabel(r'{}'.format(y))


def plot_profiles(profs, x, y, **kwargs):
    for i in range(len(profs)):
        plt.plot(profs[i].get(x), profs[i].get(y),
                 c=plt.get_cmap('inferno')(i), **kwargs)
    plt.xlabel('${{}}$'.format(x))
    plt.ylabel('${{}}$'.format(y))


def plot_HRD(mhist: MesaData, *args, **kwargs):
    plt.plot(mhist.get('log_Teff'), mhist.get('log_L'), *args, **kwargs)
    if not plt.gca().xaxis_inverted():
        plt.gca().invert_xaxis()
    plt.xlabel(r'$\log T_{\rm eff}[K]$')
    plt.ylabel(r'$\log L/L_\odot$')


def plot_kiel(mhist: MesaData, *args, **kwargs):
    plt.plot(mhist.get('log_Teff'), mhist.get('log_g'), *args, **kwargs)
    if not plt.gca().yaxis_inverted():
        plt.gca().invert_yaxis()
    if not plt.gca().xaxis_inverted():
        plt.gca().invert_xaxis()
    plt.xlabel(r'$\log T_{\rm eff}[K]$')
    plt.ylabel(r'$\log g [cgs]$')


def plot_T_rho(mprof: MesaData, *args, **kwargs):
    plt.plot(mprof.get('logRho'), mprof.get('logT'), *args, **kwargs)
    plt.ylabel(r'$\log T/$K')
    plt.xlabel(r'$\log \frac{\rho}{\rm g/cm^3}$')


def plot_T_rho_lines():
    plt.plot([-8, 2],
             [(np.log10(3.2e7) + (-8 - np.log10(0.7)) / 3.0),
              (np.log10(3.2e7) + (2 - np.log10(0.7)) / 3.0)],
             color='gray', ls='-.', alpha=0.7)
    hburn = np.loadtxt(
        '/Users/matthiasf/software/mesa-r15140/data/star_data/plot_info'
        '/hydrogen_burn.data')
    plt.plot(hburn[hburn[:, 0] < 3, 0], hburn[hburn[:, 0] < 3, 1],
             color='gray', ls='-.', alpha=0.7)


def sort_human(ll):
    def convert(text): return float(text) if text.isdigit() else text

    def alphanum(key): return [convert(c) for c in re.split('([0-9]*)', key)]

    return sorted(ll, key=alphanum)


def get_binary_profs(ddir):
    if ddir[-1] != '/':
        ddir += '/'
    return mesa_prof_data(get_profile_paths(ddir + "LOGS1")), \
        mesa_prof_data(get_profile_paths(ddir + "LOGS2"))


def get_binary_hist(ddir, read_data_cols=None):
    if ddir[-1] != '/':
        ddir += '/'
    primhist = MesaData(glob.glob(ddir + 'LOGS1/history.data')[0], read_data_cols=read_data_cols)
    sechist = MesaData(glob.glob(ddir + 'LOGS2/history.data')[0], read_data_cols=read_data_cols)
    return primhist, sechist


def determine_contact(primdata):
    if primdata is None:
        return
    pcon = primdata.get('rl_relative_overflow_1') >= 0
    scon = primdata.get('rl_relative_overflow_2') >= 0
    return np.logical_and(pcon, scon)


def age_intervals(data, dt):
    ages = data.get('star_age')
    k = 0
    m = np.zeros(len(ages), dtype=bool)
    for j in range(len(ages)):
        if ages[j] >= k * dt - ages[0]:
            m[j] = True
            k += 1
    return m


def make_contact_HRD(primdata, secdata, primcontact, primmarks, seccontact,
                     secmarks, label, c1, c2):
    plot_HRD(primdata, c1 + '--', label='prim ' + label)
    plot_HRD(secdata, c2 + '--', label='sec ' + label)
    plt.plot(np.ma.masked_where(primcontact, primdata.get('log_Teff')),
             np.ma.masked_where(primcontact, primdata.get('log_L')),
             c1, linewidth=2)
    plt.plot(np.ma.masked_where(seccontact, secdata.get('log_Teff')),
             np.ma.masked_where(seccontact, secdata.get('log_L')),
             c2, linewidth=2)
    plt.plot(primdata.get('log_Teff')[primmarks],
             primdata.get('log_L')[primmarks], 'ko', markerfacecolor='none')
    plt.plot(secdata.get('log_Teff')[secmarks], secdata.get('log_L')[secmarks],
             'ko', markerfacecolor='none')


def plot_contact_HRD(primdata, secdata, label='', c1='r', c2='b'):
    make_contact_HRD(primdata, secdata, determine_contact(primdata),
                     age_intervals(primdata, 1e5),
                     determine_contact(secdata),
                     age_intervals(secdata, 1e5),
                     label, c1, c2)


def plot_mtrans(hist, *args, **kwargs):
    plot(hist, 'm', 'lg_mtransfer_rate', *args, **kwargs)
    plt.ylim((-10, -1))


def plot_period(hist, *args, **kwargs):
    plot(hist, 'star_age', 'period_days', *args, **kwargs)


def find_RL_index(mprof, mhist, star=1):
    try:
        rl = mprof.header['rl'] / 6.957e10
    except KeyError:
        rl = mhist.get('rl_{}'.format(star))[
            mhist.get('model_number') == mprof.header['model_number']]
    ix = np.argmin(abs(rl - mprof.get('radius')))
    return ix


def plot_q(data, *args, **kwargs):
    plot(data, 'star_age', 'q', *args, **kwargs)


def routrl(qq):
    sig = 62.9237 / (15.9839 + qq ** 0.2240)
    return 1 + 3.3752 / (1 + ((np.log(qq) + 1.0105) / sig) ** 2) / (
            9.0087 + qq ** -0.4022)


def plot_radius_evolution(data, *args, prim=True, xscale=1, **kwargs):
    if prim:
        rl = data.get('rl_1')
    else:
        rl = data.get('rl_2')
    plt.plot(data.get('star_age') / xscale, data.get('radius') / rl, *args,
             **kwargs)
