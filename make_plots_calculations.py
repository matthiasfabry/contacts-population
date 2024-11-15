"""
Makes most of the plots of the paper. It first builds the population from the model grids (or loads the npz files if
you don't want to wait a bunch). Then it performs the population synthesis by computing the appropriate weighted sums
to build the 2D probability densities.
If you use the `npz` files, these are already weighted using the priors in the paper.
"""
import completion as comp
import mesa_population as mp
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pop_synth as ps
import RLOF_at_ZAMS_figures as zamsrlof
import observational_data as obs
import scipy.stats as sst
import matplotlib.colors as mcol

plt.style.use('mesa')  # comment out if you don't have this style defined

# %%
inferno = plt.colormaps['inferno']
red = inferno(0.6)
orange = inferno(0.8)
blue = plt.colormaps['viridis'](0.4)
mymap = mcolors.ListedColormap(inferno.colors[int(0.2 * inferno.N):])


def common_pq_labels(aaxs):
    aaxs[1, 3].set_ylabel(r'$\mathcal{P}_{\rm contact}(q_{\rm obs}, p_{\rm obs})$')
    aaxs[0, 0].set_ylabel(r'$\mathcal{P}_{\rm contact}(q_{\rm obs})$')
    aaxs[1, 0].set_xlabel(r'$q_{\rm obs}$')
    aaxs[1, 1].set_xlabel(r'$q_{\rm obs}$')
    aaxs[1, 2].set_xlabel(r'$q_{\rm obs}$')
    aaxs[1, 0].set_ylabel(r'$p_{\rm obs}$ (d)')
    aaxs[1, 0].set_yscale('log')
    aaxs[0, 0].set_yscale('log')
    aaxs[1, 0].yaxis.set_major_formatter(ps.log_p_fmt)
    aaxs[1, 0].yaxis.set_minor_formatter(ps.log_p_fmt)
    aaxs[1, 0].xaxis.set_label_coords(0.5, -0.1)
    aaxs[1, 1].xaxis.set_label_coords(0.5, -0.1)
    aaxs[1, 2].xaxis.set_label_coords(0.5, -0.1)


# %%  get all history files, might take a while
columns = ['model_number', 'star_age', 'stopping_condition', 'had_RLOF_at_ZAMS', 'had_contact',
           'q', 'log_dt', 'log_L', 'star_1_mass', 'star_2_mass', 'period_days',
           'rl_relative_overflow_1', 'rl_relative_overflow_2', 'lg_mtransfer_rate', 'kh_timescale'
           ]

pop = mp.MesaPopulation('full_grid/', structure='*/*/*/*/', columns=columns)
have_population = True

#%%  completion for a couple M1s (Fig. 1). Classes are overplotted using Inkscape

f, axs = plt.subplots(1, 2, figsize=(6.76, 2.535), sharey='all',
                      gridspec_kw={'wspace': 0.05, 'left': 0.08, 'bottom': 0.2, 'right': 0.75})

comp.plot_completion_p_q(pop, 26.0, 'ET', axs[0], do_legend=False, do_x=False, multi_panel=False)
plt.title(r'$M_{{1, \rm init}} = 26 M_\odot$', rotation=0)
plt.xlabel(r'$q_{{\rm init}}$')
plt.xticks(pop.qs[::4])
plt.xticks(pop.qs[2::4], minor=True)
comp.plot_completion_p_q(pop, 62.5, 'ET', axs[1], do_legend=True, do_x=False, do_y=False,
                         legend_fs=6, multi_panel=False)
plt.title(r'$M_{{1, \rm init}} = 62.5 M_\odot$', rotation=0)
plt.xlabel(r'$q_{{\rm init}}$')
plt.xticks(pop.qs[::4])
plt.xticks(pop.qs[2::4], minor=True)
plt.savefig('pngs/completion.png', dpi=300)

# %%  make completion plot for each m1_init (Figs B.1 B.2 B.3)
plt.switch_backend('Agg')  # do the following plots in non-interactive backend,
# otherwise we're limited by screen size

modes = ['ET', 'no_ET', 'single_rot']
for et in modes:
    figs = 2
    panels = len(pop.m1s) // figs
    for j in range(figs):
        fig, axs = plt.subplots(1, panels, figsize=(10, 6), sharey='all', gridspec_kw={'wspace': 0})
        axs = np.atleast_1d(axs)
        plt.suptitle(r'$M_{{1, \rm init}} / M_\odot$, {}'.format(et.replace('_', ' ')))
        for i in range(panels):
            m1 = sorted(pop.m1s)[j * panels + i]
            comp.plot_completion_p_q(pop, m1, et, ax=axs[i],
                                     do_legend=i == panels - 1, do_y=i == 0, do_x=i == panels // 2)
        for ax in axs:
            for axis in ['top', 'bottom', 'left', 'right']:
                ax.spines[axis].set_linewidth(0.5)

        plt.savefig('pngs/comp_{}_{}.png'.format(et, j), dpi=300)
        plt.close(fig)

plt.switch_backend('MacOSX')  # or another backend for your system

# %% zams rlof evolution plot (Fig. A.1)
zamsrlof.make_rlof_zams_plot()
plt.savefig('pngs/p_q_zamsrlof_evolution.png', dpi=300)

# %% Here the real population synthesis starts

# 1) zams rlof mass-ratio distribution
# define bins we are going to use for pop synth
qbinsize = 0.025
n_q = int((1.0 - 0.4) / qbinsize) + 1 + 1  # there's one extra "edge"
qbin_edges = np.linspace(0.4, 1.0, n_q)
qbins = qbin_edges[:-1] + qbinsize / 2

logpbinsize = 0.04
n_logp = int((np.log10(3) - np.log10(0.5)) / logpbinsize) + 1 + 1
pbin_edges = np.geomspace(0.5, 3, n_logp)
pbins = 10 ** (np.log10(pbin_edges[:-1]) + logpbinsize / 2)

# %% compute distributions from models and bins
if have_population:
    pq_zams_et, pq_zams_no_et, pq_zams_single_rot = ps.get_population_summed_weighted_probabilities(
        pop, qbin_edges, pbin_edges, ps.get_q_contact,
        ps.get_p_contact, ps.not_rlof_or_L2OF_at_ZAMS)

    # save for quicker access later
    np.savez('npzs/pq_zams_et', pq_zams_et)
    np.savez('npzs/pq_zams_no_et', pq_zams_no_et)
    np.savez('npzs/pq_zams_single_rot', pq_zams_single_rot)
else:
    # or load the npz files that are pre-computed
    pq_zams_et = np.load('zams_pq_et.npz')['arr_0']
    pq_zams_no_et = np.load('zams_pq_no_et.npz')['arr_0']
    pq_zams_single_rot = np.load('zams_pq_single_rot.npz')['arr_0']
# %%
# plot pq diagram
axs = ps.plot_distributions([pq_zams_et, pq_zams_no_et, pq_zams_single_rot], qbins, pbins,
                            pbin_edges,
                            need_to_regularize_y=True, cmap=mymap,
                            levels=np.array([0.99, 0.97, 0.9]))

common_pq_labels(axs)
axs[1, 0].set_xlim((0.56, 1.024))
axs[1, 0].set_ylim((0.46, 2.3))
axs[0, 0].set_ylim((7e-4, 2.1))
plt.savefig('pngs/pq_dist_zams_rlof.png', dpi=300)
# %% cumulative
plt.figure()
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_zams_et, 'r', label='ET')
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_zams_no_et, 'b', label='no ET',
                          zorder=-1)
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_zams_single_rot, 'm',
                          label='single rotation', zorder=-2)
plt.gca().set_xlim((0.54, 1.02))
plt.gca().set_ylim((1e-4, 1.5))
plt.gca().set_yscale('log')
plt.legend()
plt.xlabel(r'$q_{\rm obs}$')
plt.ylabel(r'CDF$(q_{\rm obs})$')
plt.savefig('pngs/q_dist_zams_rlof.png', dpi=300)

# %%  wider range of the bins for the rest of the analysis
qbinsize = 0.025
qmin = 0.2
qmax = 1
n_q = int((qmax - qmin) / qbinsize) + 2
qbin_edges = np.linspace(qmin, qmax, n_q)
qbins = qbin_edges[:-1] + qbinsize / 2

logpbinsize = 0.04
pmin = 0.5
pmax = 8
n_logp = int((np.log10(pmax) - np.log10(pmin)) / logpbinsize) + 2
pbin_edges = np.geomspace(pmin, pmax, n_logp)
pbins = 10 ** (np.log10(pbin_edges[:-1]) + logpbinsize / 2)


# %%
# distributions for models that do not RLOF at zams
if have_population:
    pq_et, pq_no_et, pq_single_rot = ps.get_population_summed_weighted_probabilities(
        pop, qbin_edges, pbin_edges, ps.get_q_contact,
        ps.get_p_contact, ps.rlof_or_L2OF_at_ZAMS)

    np.savez('npzs/pq_et', pq_et)
    np.savez('npzs/pq_no_et', pq_no_et)
    np.savez('npzs/pq_single_rot', pq_single_rot)
else:
    pq_et = np.load('npzs/pq_et.npz')['arr_0']
    pq_no_et = np.load('npzs/pq_no_et.npz')['arr_0']
    pq_single_rot = np.load('npzs/pq_single_rot.npz')['arr_0']
# %% PDF plot

plt.figure()
counts, bins = np.histogram(obs.get_q(obs.alldata), bins=qbin_edges)
plt.bar(bins[:-1], counts, width=bins[-1] - bins[-2], align='edge', color='r')
plt.xlabel(r'$q$')
plt.ylabel(r'Count')
plt.gca().spines['left'].set_color('r')
plt.gca().spines['right'].set_color('r')

plt.gca().tick_params(which="both", axis='y', colors='red')
plt.gca().yaxis.label.set_color('red')
# %%
plt.twinx()
new_pq_no_et = ps.collapse_xbins(pq_no_et, 2)
new_pq_et = ps.collapse_xbins(pq_et, 2)
ps.make_marg_y_plot(plt.gca(), qbins, pbin_edges, new_pq_no_et,
                    color=mcol.colorConverter.to_rgba((0, 0, 0), 0), edgecolor='b',
                    label='no energy transfer')
ps.make_marg_y_plot(plt.gca(), qbins, pbin_edges, new_pq_et,
                    color=mcol.colorConverter.to_rgba((0, 0, 0), 0), edgecolor='goldenrod',
                    label='with energy transfer')
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['right'].set_color('goldenrod')
plt.gca().tick_params(which="both", axis='y', colors='goldenrod')
plt.legend()
plt.ylabel(r'Probability')
plt.gca().yaxis.label.set_color('goldenrod')
# %%
axs = ps.plot_distributions([pq_et, pq_no_et, pq_single_rot], qbins, pbins, pbin_edges,
                            need_to_regularize_y=True, cmap=mymap)

common_pq_labels(axs)
axs[1, 0].set_xlim((0.21, 1.024))
axs[1, 0].set_ylim((0.55, 8.1))
axs[0, 0].set_ylim((7e-4, 3))
for ax in axs[1, 1:-1]:
    plt.setp(ax.get_yminorticklabels(), visible=False)
cs = []
for data, color in [(obs.mwdata, 'r'), (obs.mcdata, 'b')]:
    qs = obs.get_q(data)
    for ax in axs[1, :-1]:
        cs.append(ax.errorbar(qs, obs.get_p(data), xerr=obs.get_q_error(data),
                              ms=3, fmt='o', elinewidth=0.5, mfc='w', mew=0.5, mec=color,
                              ecolor=color))
axs[1, 0].legend([cs[0], cs[3]], ['MW data', 'MC data'], fontsize=5, loc='lower left')
plt.savefig('pngs/pq_dist.png', dpi=300)
#%%
print(sst.pmean(qbins, 1, weights=ps.marginalize_over_y(qbins, pbin_edges, pq_et)))
print(sst.pmean(qbins, 1, weights=ps.marginalize_over_y(qbins, pbin_edges, pq_no_et)))
print(sst.pmean(qbins, 1, weights=ps.marginalize_over_y(qbins, pbin_edges, pq_single_rot)))
# %%
lbinsize = 0.025
lmin = 0.05
lmax = 2.0
n_l = int((lmax - lmin) / lbinsize) + 1
lbin_edges = np.linspace(lmin, lmax, n_l)
lbins = lbin_edges[:-1] + lbinsize / 2
#%%
lq_et, lq_no_et, lq_single_rot = ps.get_population_summed_weighted_probabilities(
    pop, qbin_edges, lbin_edges, ps.get_q_contact,
    ps.get_lum_ratio_contact, ps.rlof_or_L2OF_at_ZAMS
)

np.savez('npzs/lq_et', lq_et)
np.savez('npzs/lq_no_et', lq_no_et)
np.savez('npzs/lq_single_rot', lq_single_rot)

#%%
lq_et = np.load('npzs/lq_et.npz')['arr_0']
lq_no_et = np.load('npzs/lq_no_et.npz')['arr_0']
lq_single_rot = np.load('npzs/lq_single_rot.npz')['arr_0']
# %%
axs = ps.plot_distributions([lq_et, lq_no_et, lq_single_rot], qbins, lbins, lbin_edges,
                            cmap=mymap, levels=None, marg=False, figsize=(6.76, 1.8))

for i in range(3):
    l1, = axs[0, i].plot(qbins, qbins, 'r--')
    l2, = axs[0, i].plot(qbins, qbins ** 2.5, 'b--')

axs[0, -1].set_ylabel(r'$\mathcal{P}_{\rm contact}(q_{\rm obs}, L_2/L_1)$')
axs[0, 0].set_yscale('linear')
axs[0, 0].set_ylabel(r'$L_2/L_1$')
for ax in axs[0, :-1]:
    ax.set_xlabel(r'$q_{\rm obs}$')

# overplot observational data
cs = []
for data, color in [(obs.mwdata, 'r'), (obs.mcdata, 'b')]:
    # data = data[np.logical_and(data['rrl1'] > 1, data['rrl2'] > 1)]
    qs = obs.get_q(data)
    for ax in axs[0, :-1]:
        cs.append(ax.errorbar(qs, obs.get_lum_ratio(data),
                              xerr=obs.get_q_error(data),
                              yerr=obs.get_lr_error(data),
                              ms=3, fmt='o', elinewidth=0.5,
                              mfc='w', mew=0.5, mec=color, ecolor=color))
axs[0, 0].legend([cs[0], cs[len(axs[0])], l1, l2], ['MW data', 'MC data', r'$q^1$', r'$q^{2.5}$'])
axs[0, 0].set_xlim((0.26, 1.05))
axs[0, 0].set_ylim((-0.05, 1.72))
axs[0, 0].xaxis.set_label_coords(0.5, -0.1)
axs[0, 1].xaxis.set_label_coords(0.5, -0.1)
axs[0, 2].xaxis.set_label_coords(0.5, -0.1)

plt.savefig('pngs/lq_dist.png', dpi=300)


# %% "high mass only"

def highmass(hhist):
    return hhist.get('star_1_mass')[0] < 19.9 or ps.rlof_or_L2OF_at_ZAMS(hhist)


pq_et_high_m, pq_no_et_high_m, pq_single_rot_high_m = \
    ps.get_population_summed_weighted_probabilities(
        pop, qbin_edges, pbin_edges, ps.get_q_contact, ps.get_p_contact, highmass)

# np.savez('npzs/pq_et_high_m', pq_et_high_m)
# np.savez('npzs/pq_no_et_high_m', pq_no_et_high_m)
# np.savez('npzs/pq_single_rot_high_m', pq_single_rot_high_m)
#%%
pq_et_high_m = np.load('npzs/pq_et_high_m.npz')['arr_0']
pq_no_et_high_m = np.load('npzs/pq_no_et_high_m.npz')['arr_0']
pq_single_rot_high_m = np.load('npzs/pq_single_rot_high_m.npz')['arr_0']
# %%
axs = ps.plot_distributions([pq_et_high_m, pq_no_et_high_m, pq_single_rot_high_m],
                            qbins, pbins, pbin_edges,
                            need_to_regularize_y=True, cmap=mymap,
                            levels=np.array([0.99, 0.95, 0.9, 0.7, 0.5]))
common_pq_labels(axs)
axs[1, 0].set_xlim((0.26, 1.024))
axs[1, 0].set_ylim((0.7, 7.9))
axs[0, 0].set_ylim((7e-4, 2.2))
for ax in axs[1, 1:-1]:
    plt.setp(ax.get_yminorticklabels(), visible=False)

# %% overplot observational data

cs = []
for data, color in [(obs.mwdata[obs.mwdata['m1'] > 19.9], 'r'),
                    (obs.mcdata[obs.mcdata['m1'] > 19.9], 'b')]:
    qs = obs.get_q(data)
    for ax in axs[1, :-1]:
        cs.append(ax.errorbar(qs, obs.get_p(data),
                              xerr=obs.get_q_error(data),
                              yerr=obs.get_lr_error(data),
                              ms=3, fmt='o', elinewidth=0.5,
                              mfc='w', mew=0.5, mec=color, ecolor=color))
axs[1, 0].legend([cs[0], cs[len(axs[0])]], ['MW data', 'MC data'], fontsize=5)
plt.savefig('pngs/pq_dist_high_mass.png', dpi=300)
# %% cumulative
plt.figure()
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_et, 'r', label='ET')
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_no_et, 'b', label='no ET',
                          zorder=-1)
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_single_rot, 'm',
                          label='single rotation', zorder=-2)
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_et_high_m, 'r--',
                          label=r'ET, $M_{1, \rm init} \geq 20M_\odot$')
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_no_et_high_m, 'b--',
                          label=r'no ET, $M_{1, \rm init} \geq 20M_\odot$')
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_single_rot_high_m, 'm--',
                          label=r'single rot., $M_{1, \rm init} \geq 20M_\odot$')
plt.gca().set_xlim((0.23, 1.02))
plt.gca().set_ylim((1e-2, 1.5))
plt.gca().set_yscale('log')
obsqs_all = np.sort(obs.get_q(obs.alldata))
obsqs = np.sort(obs.get_q(obs.alldata[obs.alldata['m1']>19.9]))
obsls = np.sort(obs.get_lum_ratio(obs.alldata[obs.alldata['m1']>19.9]))
plt.gca().ecdf(obsqs_all, color='gray', label='observed systems')
plt.gca().ecdf(obsqs, color='gray', ls='dashed', label=r'observed systems, $M_1 > 20 M_\odot$')

# qq = np.linspace(0.2, 1, 100)
# plt.plot(qq, qq ** 3, color=orange, label=r'$q^3$', zorder=-1)
plt.legend(fontsize=5)
plt.xlabel(r'$q_{\rm obs}$')
plt.ylabel(r'CDF$(q_{\rm obs})$')
plt.savefig('pngs/q_dist_high_mass.png', dpi=300)

# %% KS tests


# create random variable that distributes per the marginalized probabilities
class QDist(sst.rv_continuous):
    def __init__(self, qbase, qdata):
        super().__init__(a=0, b=1)
        self.qdata = qdata
        self.qbase = qbase
        self.ecdf = np.cumsum(qdata)

    def _cdf(self, x, *args):
        if isinstance(x, float):
            q_ix = len(self.qbase[self.qbase < x]) - 1  # index of biggest q that is smaller than x
            return self.ecdf[q_ix]
        cdfs = np.empty(x.size)
        for i, xval in enumerate(x):
            q_ix = len(self.qbase[self.qbase < xval]) - 1
            cdfs[i] = self.ecdf[q_ix]
        return cdfs


# %%

qdist_et_high_m = QDist(qbins, ps.marginalize_over_y(qbins, pbin_edges, pq_et_high_m))
qdist_no_et_high_m = QDist(qbins, ps.marginalize_over_y(qbins, pbin_edges, pq_no_et_high_m))
qdist_single_rot_high_m = QDist(qbins, ps.marginalize_over_y(qbins, pbin_edges, pq_single_rot_high_m))
qdist_et = QDist(qbins, ps.marginalize_over_y(qbins, pbin_edges, pq_et))
qdist_no_et = QDist(qbins, ps.marginalize_over_y(qbins, pbin_edges, pq_no_et))
qdist_single_rot = QDist(qbins, ps.marginalize_over_y(qbins, pbin_edges, pq_single_rot))
ldist_et = QDist(lbins, ps.marginalize_over_x(qbin_edges, lbins, lq_et))
ldist_no_et = QDist(lbins, ps.marginalize_over_x(qbin_edges, lbins, lq_no_et))
ks_et_high_m = sst.ks_1samp(obsqs, qdist_et_high_m.cdf, method='exact')
ks_no_et_high_m = sst.ks_1samp(obsqs, qdist_no_et_high_m.cdf, method='exact')
ks_single_rot_high_m = sst.ks_1samp(obsqs, qdist_single_rot_high_m.cdf, method='exact')
ks_no_et = sst.ks_1samp(obsqs, qdist_no_et.cdf, method='exact')
ks_et = sst.ks_1samp(obsqs, qdist_et.cdf, method='exact')
ks_single_rot = sst.ks_1samp(obsqs, qdist_single_rot.cdf, method='exact')
ks_l_et = sst.ks_1samp(obsls, ldist_et.cdf, method='exact')
ks_l_no_et = sst.ks_1samp(obsls, ldist_no_et.cdf, method='exact')
print('ET high mass only,', ks_et_high_m.statistic, ks_et_high_m.pvalue)
print('no ET high mass only,', ks_no_et_high_m.statistic, ks_no_et_high_m.pvalue)
print('single rot high mass only,', ks_single_rot_high_m.statistic, ks_single_rot_high_m.pvalue)
print('ET all mass,', ks_et.statistic, ks_et.pvalue)
print('no ET all mass,', ks_no_et.statistic, ks_no_et.pvalue)
print('single rot all mass,', ks_single_rot.statistic, ks_single_rot.pvalue)
print('lum ratio ET', ks_l_et.statistic, ks_l_et.pvalue)
print('lum ratio no ET', ks_l_no_et.statistic, ks_l_no_et.pvalue)


# %% high p; do pcut search
pcuts = np.array([0.95, 1.05, 1.15, 1.26, 1.38, 1.52, 1.66])
pvals_et = np.empty(pcuts.size)
pvals_no_et = np.empty(pcuts.size)
pvals_single_rot = np.empty(pcuts.size)
for i, pcut in enumerate(pcuts):
    def highperiod(hhist):
        return hhist.get('period_days')[0] <= pcut or ps.rlof_or_L2OF_at_ZAMS(hhist)


    pq_et_high_p, pq_no_et_high_p, pq_single_rot_high_p = \
        ps.get_population_summed_weighted_probabilities(
            pop, qbin_edges, pbin_edges, ps.get_q_contact, ps.get_p_contact, highperiod)

    qdist_et_high_p = QDist(qbins, ps.marginalize_over_y(qbins, pbin_edges, pq_et_high_p))
    qdist_no_et_high_p = QDist(qbins, ps.marginalize_over_y(qbins, pbin_edges, pq_no_et_high_p))
    qdist_single_rot_high_p = QDist(qbins, ps.marginalize_over_y(qbins, pbin_edges, pq_single_rot_high_p))
    ks_et_high_p = sst.ks_1samp(obsqs, qdist_et_high_p.cdf, method='exact')
    ks_no_et_high_p = sst.ks_1samp(obsqs, qdist_no_et_high_p.cdf, method='exact')
    ks_single_rot_high_p = sst.ks_1samp(obsqs, qdist_single_rot_high_p.cdf, method='exact')
    print(pcut)
    print('ET high period only,', ks_et_high_p.statistic, ks_et_high_p.pvalue)
    pvals_et[i] = ks_et_high_p.pvalue
    print('no ET high period only,', ks_no_et_high_p.statistic, ks_no_et_high_p.pvalue)
    pvals_no_et[i] = ks_no_et_high_p.pvalue
    print('single rot high period only,', ks_single_rot_high_p.statistic, ks_single_rot_high_p.pvalue)
    pvals_single_rot[i] = ks_single_rot_high_p.pvalue

# %%
plt.figure()
plt.plot(pcuts, pvals_et, 'r', label='ET')
plt.plot(pcuts, pvals_no_et, 'b', label='no ET')
plt.plot(pcuts, pvals_single_rot, 'm', label='single rotation')
plt.legend()
plt.xlabel(r'$p_{\rm cut}({\rm d})$')
plt.ylabel(r'KS $p$-value')
plt.savefig('pngs/p-values.png', dpi=300)
# %%  use best pcut
pcut = 1.26


def highperiod(hhist):
    return hhist.get('period_days')[0] <= pcut or ps.rlof_or_L2OF_at_ZAMS(hhist)


pq_et_high_p, pq_no_et_high_p, pq_single_rot_high_p = \
    ps.get_population_summed_weighted_probabilities(
        pop, qbin_edges, pbin_edges, ps.get_q_contact, ps.get_p_contact, highperiod)

# %%
np.savez('npzs/pq_et_high_p', pq_et_high_p)
np.savez('npzs/pq_no_et_high_p', pq_no_et_high_p)
np.savez('npzs/pq_single_rot_high_p', pq_single_rot_high_p)
# %%
axs = ps.plot_distributions([pq_et_high_p, pq_no_et_high_p, pq_single_rot_high_p],
                            qbins, pbins, pbin_edges,
                            need_to_regularize_y=True, cmap=mymap,
                            levels=np.array([0.99, 0.95, 0.9, 0.7, 0.5]))
common_pq_labels(axs)
axs[1, 0].set_xlim((0.24, 1.024))
axs[1, 0].set_ylim((0.7, 7.9))
axs[0, 0].set_ylim((3e-3, 1.2))

cs = []
for data, color in [(obs.mwdata, 'r'), (obs.mcdata, 'b')]:
    qs = obs.get_q(data)
    for ax in axs[1, :-1]:
        cs.append(ax.errorbar(qs, obs.get_p(data), xerr=obs.get_q_error(data),
                              ms=3, fmt='o', elinewidth=0.5, mfc='w', mew=0.5, mec=color,
                              ecolor=color))
axs[1, 0].legend([cs[0], cs[len(axs[0])]], ['MW data', 'MC data'], loc='lower left', fontsize=5)
for ax in axs[1, 1:-1]:
    plt.setp(ax.get_yminorticklabels(), visible=False)
plt.savefig('pngs/pq_dist_pcut.png', dpi=300)
# %% cumulative
plt.figure()
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_et_high_p, 'r--',
                          label=r'ET, $p_{{\rm init}} > {}{{\rm d}}$'.format(pcut))
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_no_et_high_p, 'b--',
                          label=r'no ET, $p_{{\rm init}} > {}{{\rm d}}$'.format(pcut))
ps.make_marg_cumul_y_plot(plt.gca(), qbin_edges[1:], pbin_edges, pq_single_rot_high_p, 'm--',
                          label=r'single rot., $p_{{\rm init}} > {}{{\rm d}}$'.format(pcut))
plt.gca().set_xlim((0.23, 1.02))
plt.gca().set_ylim((3e-2, 1.3))
plt.gca().set_yscale('log')
obsqs = np.sort(obs.get_q(obs.alldata))
plt.gca().ecdf(obsqs, color='gray', label='observed systems')
qq = np.linspace(0.2, 1, 100)
plt.legend()
plt.xlabel(r'$q_{\rm obs}$')
plt.ylabel(r'CDF$(q_{\rm obs})$')
plt.savefig('pngs/q_dist_pcut.png', dpi=300)
# %%
plt.figure()
opik = 1 / pbins / (np.log10(8. / 0.5))
plt.plot(pbins, opik, label=r'\"Opik\'s law')
plt.ylabel('Probability')
plt.xlabel(r'$P_{\rm init} ({\rm d})$')
# %%
plt.fill_between(pbins, 0, opik, where=pbins <= 1.4, alpha=0.6, edgecolor='b')
plt.fill_between(pbins, 0, opik, where=np.logical_and(pbins >= 1.3, pbins < 2.0), alpha=0.6,
                 facecolor='red', edgecolor='r')
# %%
maxwell = 2 * pbins ** 2 / 8 * np.exp(-pbins ** 2 / 8)
plt.plot(pbins, maxwell, zorder=10)
plt.fill_between(pbins, 0, maxwell, where=pbins <= 2.1, alpha=0.6, edgecolor='b', zorder=9)
plt.fill_between(pbins, 0, maxwell, where=np.logical_and(pbins >= 2., pbins < 6), alpha=0.6,
                 facecolor='red', edgecolor='r', zorder=9)
# %%

steps = 0

for model in pop.hists.flatten():
    steps += model.get('model_number')[-1]

print(steps)
print(steps / len(pop.hists.flatten()))
print(steps / (20 * 365.25))
