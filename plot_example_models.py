"""
Makes figures 2-3-4-B.1-C.1 by plotting the mass transfer rate, radius and mass-ratio evolution of selected systems.
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
case1 = glob.glob('cases/24.0/1.15/0.7/')[0]
plot_case(*get_hists(case1), name='case1')
# %%  class II
case2 = glob.glob('cases/30.0/1.52/0.9/')[0]
plot_case(*get_hists(case2), name='case2')
# %%  class IV
case4 = glob.glob('cases/62.5/1.52/0.825/')[0]
hists = get_hists(case4)
plot_case(*hists, name='case4', show_wind=True)
# %% appendix models
def plot_whole(prim_ET, sec_ET, prim_NO_ET, sec_NO_ET, prim_whole, sec_whole):
    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1, hspace=0.0, height_ratios=(3, 2, 2))
    axs = np.empty((3, 1), dtype=object)
    myr = 1e6

    # mdot
    axs[0, 0] = plt.subplot(gs[0, 0])
    md.plot(prim_ET, 'star_age', 'lg_mtransfer_rate', 'r', label='ET', xscale=myr)
    plt.ylim(bottom=-8, top=-2)
    plt.setp(axs[0, 0].get_xticklabels(), visible=False)
    # radius evolution
    axs[1, 0] = plt.subplot(gs[1, 0], sharex=axs[0, 0])
    md.plot_radius_evolution(prim_ET, 'r', xscale=myr)
    md.plot_radius_evolution(sec_ET, 'r--', prim=False, xscale=myr)
    plt.hlines(1, min(prim_ET.get('star_age') / 1e6), 1.1 * max(prim_ET.get('star_age') / 1e6), color='gray',
               linewidth=1)
    plt.setp(axs[1, 0].get_xticklabels(), visible=False)

    axs[2, 0] = plt.subplot(gs[2, 0], sharex=axs[1, 0])
    md.plot_q(prim_ET, 'r', xscale=myr)
    plt.hlines(1, min(prim_ET.get('star_age') / 1e6), 1.1 * max(prim_ET.get('star_age') / 1e6), color='gray',
               linewidth=1)

    # NO ET
    # mdot
    plt.sca(axs[0, 0])
    md.plot(prim_NO_ET, 'star_age', 'lg_mtransfer_rate', 'b', label='no ET', xscale=myr, zorder=-2)

    # radius evolution, mass ratio
    plt.sca(axs[1, 0])
    md.plot_radius_evolution(prim_NO_ET, 'b', xscale=myr, zorder=-2)
    md.plot_radius_evolution(sec_NO_ET, 'b--', prim=False, xscale=myr, zorder=-2)
    plt.ylabel(r'$R/R_{\rm RL}$')

    plt.sca(axs[2, 0])
    md.plot_q(prim_NO_ET, 'b', xscale=myr, zorder=-2)

    # whole
    plt.sca(axs[0, 0])
    md.plot(prim_whole, 'star_age', 'lg_mtransfer_rate', 'orange', label='single rot.', xscale=myr)
    legend_elements = [lines.Line2D([0], [0], color='orange', label='whole envelope'),
                       lines.Line2D([0], [0], color='red', label='ET'),
                       lines.Line2D([0], [0], color='blue', label='no ET')]
    plt.ylabel(r'$\log \frac{\dot{M}}{M_\odot / {\rm yr}}$')
    plt.legend(handles=legend_elements, loc='upper right', fontsize=5)

    # radius evolution, mass ratio
    plt.sca(axs[1, 0])
    md.plot_radius_evolution(prim_whole, 'orange', xscale=myr)
    md.plot_radius_evolution(sec_whole, 'orange', ls='--', prim=False, xscale=myr)
    legend_elements = [lines.Line2D([0], [0], color='gray', label='primary'),
                       lines.Line2D([0], [0], color='gray', linestyle='--', label='secondary')]
    plt.legend(handles=legend_elements, fontsize=5)

    plt.sca(axs[2, 0])
    md.plot_q(prim_whole, 'orange', xscale=myr)
    plt.xlabel(r'age (Myr)')
    plt.ylabel(r'$q$')
    # plt.savefig('pngs/'+name)
    # plt.close()


case5 = glob.glob('cases/30.0/1.52/0.9/')[0]
prim_ET, sec_ET, prim_no_ET, sec_no_ET, _, _ = get_hists(case5)
prim_whole, sec_whole = md.get_binary_hist(glob.glob('appendices/whole_envelope/30.0/1.52/0.9/ET')[0])
age, _ = prim_whole.trim_PMS(step_req=3)
prim_whole.data['star_age'] -= age
sec_whole.trim_PMS(age_in=age)
sec_whole.data['star_age'] -= age
plot_whole(prim_ET, sec_ET, prim_no_ET, sec_no_ET, prim_whole, sec_whole)

#%%  low mass case
prim, sec = md.get_binary_hist(glob.glob('appendices/low_mass/1.0/0.9/0.6/ET/')[0])
prim_n, sec_n = md.get_binary_hist(glob.glob('appendices/low_mass/1.0/0.9/0.6/no_ET/')[0])

age, _ = prim.trim_PMS(step_req=3)
prim.data['star_age'] -= age
sec.trim_PMS(age_in=age)
sec.data['star_age'] -= age

age, _ = prim_n.trim_PMS(step_req=3)
prim_n.data['star_age'] -= age
sec_n.trim_PMS(age_in=age)
sec_n.data['star_age'] -= age

fig = plt.figure()
gs = gridspec.GridSpec(4, 2, hspace=0.0, height_ratios=(2, 2, 2, 2), width_ratios=(2, 5), wspace=0.03)
axs = np.empty((4, 2), dtype=object)
gyr = 1e9

# mdot
axs[0, 0] = plt.subplot(gs[0, 0])
md.plot(prim_n, 'star_age', 'lg_mtransfer_rate', 'b', label='no ET', xscale=gyr)
md.plot(prim, 'star_age', 'lg_mtransfer_rate', 'r', label='ET', xscale=gyr)
plt.ylim(bottom=-12, top=-5)
plt.setp(axs[0, 0].get_xticklabels(), visible=False)
plt.ylabel(r'$\log \frac{\dot{M}}{M_\odot / {\rm yr}}$')
plt.xlabel('')

axs[0, 1] = plt.subplot(gs[0, 1])
plt.gca().yaxis.set_tick_params(labelright=True, labelleft=False)
md.plot(prim_n, 'star_age', 'lg_mtransfer_rate', 'b', label='no ET', xscale=gyr)
md.plot(prim, 'star_age', 'lg_mtransfer_rate', 'r', label='ET', xscale=gyr)
plt.ylabel('')
plt.ylim(bottom=-12, top=-5)
plt.xlim(left=2.62, right=3.02)
plt.setp(axs[0, 1].get_xticklabels(), visible=False)
plt.legend(loc='upper left')

# radius evolution
axs[1, 0] = plt.subplot(gs[1, 0], sharex=axs[0, 0])
plt.hlines(1, min(prim.get('star_age') / gyr), 1.1 * max(prim.get('star_age') / gyr), color='gray', linewidth=1)
md.plot_radius_evolution(prim_n, 'b', xscale=gyr)
md.plot_radius_evolution(sec_n, 'b--', prim=False, xscale=gyr)
md.plot_radius_evolution(prim, 'r', xscale=gyr, label='primary')
md.plot_radius_evolution(sec, 'r--', prim=False, xscale=gyr, label='secondary')
plt.setp(axs[1, 0].get_xticklabels(), visible=False)
plt.ylabel(r'$R/R_{\rm RL}$')

axs[1, 1] = plt.subplot(gs[1, 1], sharex=axs[0, 1])
plt.gca().yaxis.set_tick_params(labelright=True, labelleft=False)
plt.hlines(1, min(prim.get('star_age') / gyr), 1.1 * max(prim.get('star_age') / gyr), color='gray', linewidth=1)
md.plot_radius_evolution(prim_n, 'b', xscale=gyr)
md.plot_radius_evolution(sec_n, 'b--', prim=False, xscale=gyr)
md.plot_radius_evolution(prim, 'r', xscale=gyr)
md.plot_radius_evolution(sec, 'r--', prim=False, xscale=gyr)
plt.setp(axs[1, 0].get_xticklabels(), visible=False)
plt.ylim(bottom=0.80, top=1.27)
legend_elements = [lines.Line2D([0], [0], color='gray', label='primary'),
                   lines.Line2D([0], [0], color='gray', linestyle='dotted', label='secondary')]
plt.legend(handles=legend_elements, fontsize=5, loc='upper left')

axs[2, 0] = plt.subplot(gs[2, 0], sharex=axs[1, 0])
plt.hlines(1, min(prim.get('star_age') / gyr), 1.1 * max(prim.get('star_age') / gyr), color='gray', linewidth=1)
md.plot_q(prim_n, 'b', xscale=gyr)
md.plot_q(prim, 'r', xscale=gyr)
plt.setp(axs[2, 0].get_xticklabels(), visible=False)
plt.ylabel(r'$q$')

axs[2, 1] = plt.subplot(gs[2, 1], sharex=axs[1, 1])
plt.gca().yaxis.set_tick_params(labelright=True, labelleft=False)
plt.hlines(1, min(prim.get('star_age') / gyr), 1.1 * max(prim.get('star_age') / gyr), color='gray', linewidth=1)
md.plot_q(prim_n, 'b', xscale=gyr)
md.plot_q(prim, 'r', xscale=gyr)
# thermal timescale average of q
q = prim.get('q')
dt = 10 ** prim.get('log_dt')
ix = np.argmax(q)
t_kh = prim.get('kh_timescale')[ix]
q = q[ix:]
dt = dt[ix:]
running_q = np.empty(q.size)
for j in range(len(running_q)):
    t = 0
    k = j
    while t < t_kh / 2:
        t += dt[k]
        if k == 0:
            break
        k -= 1
    m = j
    t = 0
    while t < t_kh / 2:
        t += dt[m]
        m += 1
        if m == len(running_q):
            break
    running_q[j] = np.sum(q[k:m] * dt[k:m]) / np.sum(dt[k:m])
plt.plot(prim.get('star_age')[ix:] / gyr, running_q, 'k', label=r'$\tilde{q}_{\rm KH}$')
plt.ylabel('')
plt.xlabel(r'age (Gyr)')
plt.legend(loc='upper left', fontsize=5)


axs[3, 0] = plt.subplot(gs[3, 0], sharex=axs[0, 0])
plt.sca(axs[3, 0])
md.plot(prim, 'star_age', 'period_days', 'r', xscale=gyr)
md.plot(prim_n, 'star_age', 'period_days', 'b', xscale=gyr)
plt.ylabel(r'$p/{\rm d}$')
plt.xlabel(r'age (Gyr)')
axs[3, 1] = plt.subplot(gs[3, 1], sharex=axs[0, 1])
md.plot(prim, 'star_age', 'period_days', 'r', xscale=gyr)
md.plot(prim_n, 'star_age', 'period_days', 'b', xscale=gyr)
plt.ylabel('')
plt.gca().yaxis.set_tick_params(labelright=True, labelleft=False)
plt.ylim(bottom=0.16, top=0.39)
plt.xlabel(r'age (Gyr)')

plt.subplots_adjust(bottom=0.15, left=0.17, top=0.95)
