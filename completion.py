import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pop_synth as ps
import mesa_population as mp
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

plt.style.use('mesa')

no_result_color = mcolors.CSS4_COLORS['grey']
error_color = mcolors.CSS4_COLORS['red']
survive_color = mcolors.CSS4_COLORS['gold']
initial_L2OF_color = mcolors.CSS4_COLORS['steelblue']
ms_merger_color = mcolors.CSS4_COLORS['chocolate']
no_caseA_color = mcolors.CSS4_COLORS['beige']
max_mdot_color = mcolors.CSS4_COLORS['navy']

msc_hatch = '///'
zams_rlof_hatch = '\\\\\\'

error_patch = mpatches.Patch(color=error_color, label='error')

survive_patch = mpatches.Patch(color=survive_color, label='survives MS')
initial_l2of_patch = mpatches.Patch(color=initial_L2OF_color, label='initial L2OF')
ms_merger = mpatches.Patch(color=ms_merger_color, label='merger')
no_result_patch = mpatches.Patch(color=no_result_color, label='uncompleted')
max_mdot_patch = mpatches.Patch(color=max_mdot_color,
                                label=r'$\dot{M}_{\rm trans} \geq 10^{-1} M_\odot/{\rm yr}$')
ms_contact_patch = mpatches.Patch(facecolor='w', hatch=msc_hatch, label='contact')
zams_RLOF_patch = mpatches.Patch(facecolor='w', hatch=zams_rlof_hatch, label='RLOF at ZAMS')


def my_colors(x):  # should be int between [0, clrs.N)
    if x == 0:  # some stop condition I have not covered, this is an error at this point.
        return error_color
    elif x == 1:  # survives to TAMS
        return survive_color
    elif x == 2:  # L2OF at ZAMS
        return initial_L2OF_color
    elif x == 3:  # L2OF during MS
        return ms_merger_color
    elif x == 4:  # no case A interaction
        return no_caseA_color
    elif x == 5:
        return max_mdot_color
    elif x == 6:  # convergence errors
        return error_color


vcolors = np.vectorize(my_colors)

completionclrs = vcolors(np.arange(8))
indexmap = {-1: 0, 0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, -99: 6}
completionmap = mcolors.ListedColormap(completionclrs)


def plot_completion_p_m1(p: mp.MesaPopulation, q, et_value, ax, do_legend=False,
                         do_y=True, do_x=True):
    x = p.m1s
    y = p.ps
    X, Y = np.meshgrid(x, y)
    Z = np.empty(X.shape, dtype=int)
    contact = np.empty(X.shape, dtype=bool)
    zamsrlof = np.empty(X.shape, dtype=bool)
    for i in range(len(x)):
        for j in range(len(y)):
            this = p[x[i], y[j], q, et_value, 'prim']
            if this is None:
                Z[j, i] = 4
                contact[j, i] = False
                zamsrlof[j, i] = False
            else:
                Z[j, i] = int(indexmap[this.get('stopping_condition')[-1]])
                contact[j, i] = bool(this.get('had_contact')[-1])
                zamsrlof[j, i] = bool(this.get('had_RLOF_at_ZAMS')[-1])
    plt.sca(ax)
    plt.title(r'${}$'.format(q))
    masks = [contact & zamsrlof, contact & (~zamsrlof), (~contact) & zamsrlof,
             (~contact) & (~zamsrlof)]
    hatches = ['', zams_rlof_hatch, msc_hatch, msc_hatch + zams_rlof_hatch]
    for mask, hatch in zip(masks, hatches):
        plt.pcolor(X, Y, ma.array(Z, mask=mask), hatch=hatch, cmap=completionmap,
                   edgecolors='black', vmin=0, vmax=completionmap.N)
    plt.gca().set_yscale('log')
    if do_legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                   handles=[error_patch, survive_patch, initial_l2of_patch,
                            ms_merger, max_mdot_patch,
                            ms_contact_patch, zams_RLOF_patch])
    if do_x:
        plt.xlabel(r'$M_{\rm 1,ini}/M_\odot$')
    if do_y:
        plt.ylabel(r'$p$/d')


def plot_completion_p_q(p: mp.MesaPopulation, m1_value, et_value, ax, do_legend=False,
                        do_y=True, do_x=True, multi_panel=True, legend_fs=10):
    x = p.qs
    y = p.ps
    X, Y = np.meshgrid(x, y)
    Z = np.empty(X.shape, dtype=int)
    contact = np.empty(X.shape, dtype=bool)
    zamsrlof = np.empty(X.shape, dtype=bool)
    for i in range(len(x)):
        for j in range(len(y)):
            this = p[m1_value, y[j], x[i], et_value, 'prim']
            if this is None:
                Z[j, i] = 4
                contact[j, i] = False
                zamsrlof[j, i] = False
            else:
                try:
                    Z[j, i] = int(indexmap[this.get('stopping_condition')[-1]])
                except IndexError:
                    print(this.get('stopping_condition'), m1_value, y[j], x[i])
                contact[j, i] = bool(this.get('had_contact')[-1])
                zamsrlof[j, i] = bool(this.get('had_RLOF_at_ZAMS')[-1])
    plt.sca(ax)
    plt.title(r'${}$'.format(m1_value), rotation=90)
    masks = [contact & zamsrlof, contact & (~zamsrlof), (~contact) & zamsrlof,
             (~contact) & (~zamsrlof)]
    hatches = [msc_hatch + zams_rlof_hatch, msc_hatch, zams_rlof_hatch, '']
    for mask, hatch in zip(masks, hatches):
        plt.pcolor(X, Y, ma.array(Z, mask=~mask), hatch=hatch, cmap=completionmap,
                   edgecolors='black', vmin=0, vmax=completionmap.N, linewidths=0.1)
    plt.gca().set_yscale('log')
    if do_legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                   handles=[error_patch, survive_patch, initial_l2of_patch,
                            ms_merger, max_mdot_patch,
                            ms_contact_patch, zams_RLOF_patch],
                   fontsize=legend_fs)
    if multi_panel:
        plt.xticks([])
    if do_x:
        plt.xlabel(r'$q_{\rm init} \in [0.6, 0.975]$')
    if do_y:
        plt.ylabel(r'$p_{\rm init}$/d')

    plt.gca().yaxis.set_major_formatter(ps.log_p_fmt)
    plt.gca().yaxis.set_minor_formatter(ps.log_p_fmt)


def compute_completion_rate(pop):
    print("relative number of failed models")
    count = 0.0
    for model in pop[:, :, :, "ET", "prim"].flatten():
        if model.get('stopping_condition')[-1] == -1 or model.get('stopping_condition')[-1] == -99:
            count += 1
    print(count / pop[:, :, :, "ET", "prim"].size)
    count = 0.0
    for model in pop[:, :, :, "no_ET", "prim"].flatten():
        if model.get('stopping_condition')[-1] == -1 or model.get('stopping_condition')[-1] == -99:
            count += 1
    print(count / pop[:, :, :, "no_ET", "prim"].size)
