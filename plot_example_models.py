"""
Makes figures 2-3-4 by plotting the mass transfer rate, radius and mass-ratio evolution of selected systems.
"""
import glob
import matplotlib.pyplot as plt
import numpy as np
import mesa_data as md
import matplotlib.lines as lines
import matplotlib.gridspec as gridspec
import os

plt.style.use('mesa')


# %%
def get_hists(case):
    primhist_ET, sechist_ET = md.get_binary_hist(case + 'ET')
    print(primhist_ET.get('stopping_condition')[-1])
    if os.path.exists(case + 'no_ET/same_as_ET.txt'):
        primhist_NO_ET, sechist_NO_ET = primhist_ET, sechist_ET
    else:
        primhist_NO_ET, sechist_NO_ET = md.get_binary_hist(case + 'no_ET')
    print(primhist_NO_ET.get('stopping_condition')[-1])
    primhist_single, sechist_single = md.get_binary_hist(case + 'single_rot')
    print(primhist_single.get('stopping_condition')[-1])

    age, _ = primhist_ET.trim_PMS(step_req=3)
    primhist_ET.data['star_age'] -= age
    sechist_ET.trim_PMS(age_in=age)
    sechist_ET.data['star_age'] -= age
    age, _ = primhist_NO_ET.trim_PMS(step_req=3)
    primhist_NO_ET.data['star_age'] -= age
    sechist_NO_ET.trim_PMS(age_in=age)
    sechist_NO_ET.data['star_age'] -= age
    age, _ = primhist_single.trim_PMS(step_req=3)
    primhist_single.data['star_age'] -= age
    sechist_single.trim_PMS(age_in=age)
    sechist_single.data['star_age'] -= age
    return primhist_ET, sechist_ET, primhist_NO_ET, sechist_NO_ET, primhist_single, sechist_single


def plot_case(prim_ET, sec_ET,
              prim_NO_ET, sec_NO_ET,
              prim_single, sec_single,
              name='example1', show_wind=False):
    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1, hspace=0.0, height_ratios=(3, 2, 2))
    axs = np.empty((3, 1), dtype=object)
    myr = 1e6
    # mdot
    axs[0, 0] = plt.subplot(gs[0, 0])
    md.plot(prim_ET, 'star_age', 'lg_mtransfer_rate', 'r', label='ET', xscale=myr)
    if show_wind:
        md.plot(prim_ET, 'star_age', 'lg_wind_mdot_1', 'r--', xscale=myr)
    plt.ylim(bottom=-8, top=-2)
    plt.setp(axs[0, 0].get_xticklabels(), visible=False)
    # radius evolution
    axs[1, 0] = plt.subplot(gs[1, 0], sharex=axs[0, 0])
    md.plot_radius_evolution(prim_ET, 'r', xscale=myr)
    md.plot_radius_evolution(sec_ET, 'r--', prim=False, xscale=myr)
    plt.hlines(1, min(prim_ET.get('star_age') / 1e6),
               1.1 * max(prim_ET.get('star_age') / 1e6), color='gray', linewidth=1)
    plt.setp(axs[1, 0].get_xticklabels(), visible=False)

    axs[2, 0] = plt.subplot(gs[2, 0], sharex=axs[1, 0])
    md.plot_q(prim_ET, 'r', xscale=myr)
    plt.hlines(1, min(prim_ET.get('star_age') / 1e6),
               1.2 * max(prim_ET.get('star_age') / 1e6), color='gray', linewidth=1)
    # NO ET
    # mdot
    plt.sca(axs[0, 0])
    md.plot(prim_NO_ET, 'star_age', 'lg_mtransfer_rate', 'b', label='no ET', xscale=myr, zorder=-1)
    if show_wind:
        md.plot(prim_NO_ET, 'star_age', 'lg_wind_mdot_1', 'b--', xscale=myr, zorder=-1)

    # radius evolution, mass ratio
    plt.sca(axs[1, 0])
    md.plot_radius_evolution(prim_NO_ET, 'b', xscale=myr, zorder=-1)
    md.plot_radius_evolution(sec_NO_ET, 'b--', prim=False, xscale=myr, zorder=-1)

    plt.sca(axs[2, 0])
    md.plot_q(prim_NO_ET, 'b', xscale=myr, zorder=-1)

    # single
    plt.sca(axs[0, 0])
    md.plot(prim_single, 'star_age', 'lg_mtransfer_rate', 'm', label='single rot.', xscale=myr, zorder=-2)
    legend_elements = [lines.Line2D([0], [0], color='red', label='ET'),
                       lines.Line2D([0], [0], color='blue', label='no ET'),
                       lines.Line2D([0], [0], color='magenta', label='single rot.')]
    if show_wind:
        md.plot(prim_single, 'star_age', 'lg_wind_mdot_1', 'm--', xscale=myr, zorder=-2)
        legend_elements.extend(
            [lines.Line2D([0], [0], color='gray', label=r'$\dot{M}_{\rm trans}$'),
             lines.Line2D([0], [0], color='gray', linestyle='--', label=r'$\dot{M}_{\rm wind}$')])
    plt.legend(handles=legend_elements, loc='upper right', fontsize=5)
    plt.ylabel(r'$\log \frac{\dot{M}}{M_\odot / {\rm yr}}$')

    # radius evolution, mass ratio
    plt.sca(axs[1, 0])
    md.plot_radius_evolution(prim_single, 'm', xscale=myr, zorder=-2)
    md.plot_radius_evolution(sec_single, 'm--', prim=False, xscale=myr, zorder=-2)
    plt.ylabel(r'$R/R_{\rm RL}$')
    legend_elements = [lines.Line2D([0], [0], color='gray', label='primary'),
                       lines.Line2D([0], [0], color='gray', linestyle='--', label='secondary')]
    plt.legend(handles=legend_elements, fontsize=5)

    plt.sca(axs[2, 0])
    md.plot_q(prim_single, 'm', xscale=myr, zorder=-2)

    plt.xlabel(r'age (Myr)')
    plt.ylabel(r'$q$')
    # plt.savefig('pngs/'+name)
    # plt.close()


# %%  class I
case1 = glob.glob('full_grid/26.0/1.82/0.625/')[0]
plot_case(*get_hists(case1), name='case1')
# %%  class II
case2 = glob.glob('full_grid/30.0/1.52/0.9/')[0]
plot_case(*get_hists(case2), name='case2')
# %%  class IV
case4 = glob.glob('full_grid/62.5/1.52/0.825/')[0]
hists = get_hists(case4)
plot_case(*hists, name='case4', show_wind=True)
