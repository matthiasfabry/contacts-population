"""
Uses a mesa work directory template to build two _actual_ mesa work directories. For a given initial primary mass,
period and mass ration, `make_pair` builds one with energy transfer enabled, one with energy transfer disabled.
`make_single` builds the one with energy transfer disabled and using the single-rotating star deformation corrections.
"""
import numpy as np
import shutil as sh
import os


def no_ds(ddir, names):
    ignore = []
    if '.DS_Store' in names:
        ignore.append('.DS_Store')
    return ignore


def make_pair(mm1, pp, qq, src, ddest=''):
    mm1 = str(np.round(float(mm1), 2))
    qq = str(np.round(float(qq), 3))
    pp = str(np.round(float(pp), 2))

    if len(ddest) > 0 and ddest[-1] != '/':
        ddest += '/'

    template = src.rsplit('/', maxsplit=1)[1]
    if template == '':
        template = src
    for option in ['ET', 'no_ET']:
        path = ddest + template + '/' + mm1 + '/' + pp + "/" + qq + '/' + option
        sh.copytree(src, path, dirs_exist_ok=True, ignore=no_ds)
        # create inlist with variable parameters
        new_inlist = open(path + '/inlist_extra', 'w')
        new_inlist.write("&binary_controls\n")
        new_inlist.write("m1 = " + mm1 + "d0\n")
        new_inlist.write("m2 = " + str(np.round(float(qq)*float(mm1), 4)) + "d0\n")
        new_inlist.write("initial_period_in_days = " + pp + "d0\n")
        new_inlist.write("/\n")
        new_inlist.close()
        # ET inlist
        new_inlist = open(path + '/inlist_et', 'w')
        new_inlist.write('&controls\n')
        if option == "ET":
            new_inlist.write('use_other_energy = .true.\n')
        else:
            new_inlist.write('use_other_energy = .false.\n')
        new_inlist.write("/\n")
        new_inlist.close()

    path = ddest + template + '/' + mm1 + '/' + pp + '/' + qq
    sh.copy('run_both', path)
    os.chmod(path + '/run_both', 0o777)
    # path = ddest + template + '/' + mm1 + '/' + pp
    # sh.copy('run_all', path)
    # os.chmod(path + '/run_all', 0o777)


def make_single(mm1, pp, qq, src, ddest=''):
    mm1 = str(np.round(float(mm1), 2))
    qq = str(np.round(float(qq), 3))
    pp = str(np.round(float(pp), 2))

    if len(ddest) > 0 and ddest[-1] != '/':
        ddest += '/'

    template = src.rsplit('/', maxsplit=1)[1]
    if template == '':
        template = src

    path = ddest + template + '/' + mm1 + '/' + pp + "/" + qq + "/single_rot"
    sh.copytree(src, path, dirs_exist_ok=True, ignore=no_ds)
    # create inlist with variable parameters
    new_inlist = open(path + '/inlist_extra', 'w')
    new_inlist.write("&binary_controls\n")
    new_inlist.write("m1 = " + mm1 + "d0\n")
    new_inlist.write("m2 = " + str(np.round(float(qq)*float(mm1), 4)) + "d0\n")
    new_inlist.write("initial_period_in_days = " + pp + "d0\n")
    new_inlist.write("/\n")
    new_inlist.close()


if __name__ == '__main__':
    import getopt as go
    import sys
    
    opts, args = go.getopt(sys.argv[1:], '')
    src, dest, m1, p, q = args
    
    make_pair(m1, p, q, src, dest)
