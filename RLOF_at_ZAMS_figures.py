"""
Makes figure A.1 showing the evolution of systems that have RLOF at ZAMS
"""


import glob
import matplotlib.pyplot as plt
import numpy as np
import mesa_data as md
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines

plt.style.use('mesa')


def make_rlof_zams_plot():
    cases = glob.glob('full_grid/full_grid_high_mass/grid/20.0/*/0.8/')
    cases = md.sort_human(cases)
    primhists_et = []
    primhists_no_et = []
    sechists_et = []
    sechists_no_et = []
    for i, case in enumerate(cases[:]):  # make copy so we don't delete while looping
        print(case)
        primhist, sechist = md.get_binary_hist(case + 'ET')
        if primhist.get('had_RLOF_at_ZAMS')[-1] == 1:
            primhists_et.append(primhist)
            sechists_et.append(sechist)
        else:
            cases.remove(case)
        if len(glob.glob(case + 'ET/had_no_ET.txt')) == 0:
            primhist, sechist = md.get_binary_hist(case + 'no_ET')
        if primhist.get('had_RLOF_at_ZAMS')[-1]:
            primhists_no_et.append(primhist)
            sechists_no_et.append(sechist)

    plt.figure()
    gs = gridspec.GridSpec(2, 1, hspace=0.02)
    axs = np.empty(2, dtype=object)

    # ET models
    # mass ratio
    axs[0] = plt.subplot(gs[0])
    axs[0].set_title('ET', y=0.75)
    tax = plt.twinx()
    for primhist_ET in primhists_et:
        plt.sca(tax)
        md.plot_q(primhist_ET, xscale=1e6)
        # period
        plt.sca(axs[0])
        md.plot(primhist_ET, 'star_age', 'period_days', '--', xscale=1e6)
    tax.hlines(1, min(primhists_et[-1].get('star_age') / 1e6),
               1.1 * max(primhists_et[-1].get('star_age') / 1e6), color='gray', linewidth=1)

    axs[0].set_ylabel(r'$p/{\rm d}$')
    tax.set_ylabel(r'$q$')
    lgd = plt.legend()
    ax = lgd.axes
    handles, labels = ax.get_legend_handles_labels()
    handles.append(mlines.Line2D([0], [0], linestyle='-', color='gray'))
    labels.append("$q$")
    handles.append(mlines.Line2D([0], [0], linestyle='--', color='gray'))
    labels.append(r"$p/{\rm d}$")
    lgd._legend_box = None
    lgd._init_legend_box(handles, labels)
    lgd._set_loc(lgd._loc)
    lgd.set_title(lgd.get_title().get_text())
    # NO ET
    axs[1] = plt.subplot(gs[1], sharex=axs[0])
    axs[1].set_title('no ET', y=0.75)
    tax2 = plt.twinx()
    for primhist_NO_ET in primhists_no_et:
        plt.sca(tax2)
        md.plot_q(primhist_NO_ET, xscale=1e6)

        plt.sca(axs[1])
        md.plot(primhist_NO_ET, 'star_age', 'period_days', '--', xscale=1e6)
    tax2.hlines(1, min(primhists_et[-1].get('star_age') / 1e6),
                1.1 * max(primhists_et[-1].get('star_age') / 1e6), color='gray', linewidth=1)
    axs[1].set_ylabel(r'$p/{\rm d}$')
    axs[1].set_xlabel(r'age (Myr)')
    tax2.set_ylabel(r'$q$')
    axs[1].set_xlim((1.1e-2, 5))
    axs[1].set_xscale('log')
    plt.tight_layout()



