"""
all kinds of functions useful to do population synthesis on mesa models.
you are not supposed to execute code here, import this file instead.
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as spinp
import mesa_data as md
import matplotlib.colors as mcolors

inferno = plt.colormaps['inferno']


# red = inferno(0.6)
# orange = inferno(0.8)
# blue = plt.colormaps['viridis'](0.4)


def log_p_fmt(x, pos):
    if np.round(x, 1) in [0.6, 0.8, 1, 2, 4, 6]:
        return "{}".format(np.round(x, 1), "1.1f")
    else:
        return ""


def get_t_contact(hhist):
    """
    get timesteps of all points of the model that it spent in contact,
    along with total time in contact
    :param hhist: MesaData object containing the history of the primary star's evolution
    :return: array containing dt(t) and a number Sum(dt)
    """
    if hhist is None:
        raise ValueError("history file is None")
    contact = md.determine_contact(hhist)
    dt = 10 ** hhist.get('log_dt')[contact]  # yr
    return np.array(dt), np.sum(dt)


def get_q_contact(hhist1, _=None):  # must have same calling signature as get_lum_ratio_contact
    """
    get mass ratio and mass ratio change of all points of the model that it spent in contact
    :param hhist1: MesaData object containing the history of the primary star's evolution
    :param _: dummy paramater
    :return: array containing (q(t), dq(t))
    """
    if hhist1 is None:
        raise ValueError("history file is None")
    contact = md.determine_contact(hhist1)
    qq = hhist1.get('q')[contact]
    qq = np.minimum(qq, 1. / qq)  # observed q
    dq = qq - np.minimum(np.roll(hhist1.get('q'), 1)[contact],
                         1. / np.roll(hhist1.get('q'), 1)[contact])
    return np.array((qq, dq))


def get_p_contact(hhist1, _=None):  # must have same calling signature as get_lum_ratio_contact
    """
    get period and period change of all points of the model that it spent in contact
    :param hhist1: MesaData object containing the history of the primary star's evolution
    :param _: dummy parameter
    :return: array containing (q(t), dq(t))
    """
    if hhist1 is None:
        raise ValueError("history file is None")
    contact = md.determine_contact(hhist1)
    pp = hhist1.get('period_days')[contact]  # days
    dp = pp - np.roll(hhist1.get('period_days'), 1)[contact]
    return np.array((pp, dp))


def get_lum_ratio_contact(hhist1, hhist2):
    """
    get luminosity ratio and its change of all points of the model that it spent in contact
    :param hhist1: MesaData object containing the history of the primary star's evolution
    :param hhist2: MesaData object containing the history of the secondary
    :return: array containing (L2/L1(t), d(L2/L1)(t))
    """
    if hhist1 is None:
        raise ValueError("primary history file is None")
    if hhist2 is None:
        raise ValueError("secondary history file is None")
    contact = md.determine_contact(hhist1)
    l2l1 = 10 ** (hhist2.get('log_L')[contact] - hhist1.get('log_L')[contact])
    dl = l2l1 - 10 ** (np.roll(hhist2.get('log_L'), 1)[contact] -
                       np.roll(hhist1.get('log_L'), 1)[contact])
    qq = hhist1.get('star_2_mass')[contact] / hhist1.get('star_1_mass')[contact]
    for i in range(len(qq)):
        if qq[i] > 1:  # secondary is more massive, would flip q, so also flip l2l1
            l2l1[i] = 1 / l2l1[i]
            dl[i] = l2l1[i] - 10 ** (np.roll(hhist1.get('log_L'), 1)[contact][i] -
                                     np.roll(hhist2.get('log_L'), 1)[contact][i])
    return np.array((l2l1, dl))


def get_d2t_dx_dy(hhist1, hhist2, xbin_edges, ybin_edges, x_lookup, y_lookup):
    """
    For given set of period and mass ratio bins, analyze the history files to compute how much time
    the model spent in contact in each bin.
    :param xbin_edges: array of x bin edges
    :param ybin_edges: array of y bin edges
    :param hhist1: MesaData object containing the history of the primary star's evolution
    :param hhist2: MesaData object containing the history of the seconary star's evolution
    :param x_lookup: function that returns (x, dx) of the property measured on x-axis
    :param y_lookup: function that return (y, dy) of the property measured on y-axis
    :return: 2d array of d^2t_dxdy for each of the x and y bins
    """
    x_dx = x_lookup(hhist1, hhist2)
    y_dy = y_lookup(hhist1, hhist2)
    dt, total_time = get_t_contact(hhist1)
    result = np.zeros((len(ybin_edges) - 1, len(xbin_edges) - 1))
    for i in range(len(x_dx[0])):  # timestep index
        j = len(xbin_edges) - 2  # x index
        k = len(ybin_edges) - 2  # y index

        xhere = x_dx[0][i]
        previous_x = xhere - x_dx[1][i]
        max_xhere = max(xhere, previous_x)
        min_xhere = min(xhere, previous_x)

        yhere = y_dy[0][i]
        previous_y = yhere - y_dy[1][i]
        max_yhere = max(yhere, previous_y)
        min_yhere = min(yhere, previous_y)

        while j > 0 and xbin_edges[j] > max_xhere:  # look for rightmost x bin
            j -= 1
        jmax = j
        while j > 0 and min_xhere < xbin_edges[j]:  # leftmost x bin
            j -= 1
        jmin = j

        while k > 0 and ybin_edges[k] > max_yhere:  # look for rightmost y bin
            k -= 1
        kmax = k
        while k > 0 and min_yhere < ybin_edges[k]:
            k -= 1
        kmin = k

        # calculate fractions on x axis
        xfracs = np.zeros((len(ybin_edges) - 1, len(xbin_edges) - 1))
        remaining = 1
        for j in range(jmin, jmax + 1):
            if max_xhere < xbin_edges[j + 1] or j == len(xbin_edges) - 2:
                xfracs[:, j] = remaining
            else:
                frachere = (xbin_edges[j + 1] - max(xbin_edges[j], min_xhere)) / \
                           (max_xhere - min_xhere)  # linear fraction that goes in bin j
                xfracs[:, j] = frachere
                remaining -= frachere

        # calculate fractions on y axis
        yfracs = np.zeros((len(ybin_edges) - 1, len(xbin_edges) - 1))
        remaining = 1
        for k in range(kmin, kmax + 1):
            if max_yhere < ybin_edges[k + 1] or k == len(ybin_edges) - 2:
                yfracs[k, :] = remaining
            else:
                frachere = (ybin_edges[k + 1] - max(ybin_edges[k], min_yhere)) / \
                           (max_yhere - min_yhere)  # linear fraction that goes in bin k
                yfracs[k, :] = frachere
                remaining -= frachere  # remove what we have accounted for here

        # if (jmin != jmax or kmin != kmax):  # debug
        #     print(np.sum(pfracs * qfracs))  # should add to 1
        #     print(pfracs, qfracs)
        result += dt[i] * yfracs * xfracs  # multiply fracs by dt_contact

    assert np.isclose(np.sum(result), total_time)  # make sure we counted correctly!

    # divide by bin sizes
    xsize, ysize = np.meshgrid(xbin_edges[1:] - xbin_edges[:-1], ybin_edges[1:] - ybin_edges[:-1])

    return result / (xsize * ysize)


def get_x_y_probabilities(population, xbin_edges, ybin_edges,
                          x_lookup, y_lookup, ignore_condition):
    """
    for a given MesaPopulation, and bin edges of both p and q, compute
    the probability of finding each models in the bins defined by the edges
    :param population: MesaPopulation
    :param ybin_edges: array of bin edges on y-axis
    :param xbin_edges: array of bin edges on x-axis
    :param x_lookup: function that returns (x, dx) for a model
    :param y_lookup: function that return (y, dy) for a model
    :param ignore_condition: function evaluatied on the history of the primary star that tells
        whether to ignore it in the probability density
    :return: `keys`: array of shape (len(m1s), len(ps), len(qs), 2, 4), holding (m1, p, q, et) for
        each of the simulations in the population
        and `data`: array of shape (len(m1s), len(ps), len(qs), 2, len(ppbin_edges) - 1,
        len(qqbin_edges) - 1) holding the array of probabilities to find each model
        (indexed by the first four indices), in each bin of x and y (indexed by the last 2 indices).
    """
    supershape = population.hists.shape[:-1]  # don't need prim/sec dimension
    keysshape = np.append(supershape, 4)
    datashape = np.append(supershape, [len(ybin_edges) - 1, len(xbin_edges) - 1])
    params = np.zeros(keysshape)
    data = np.zeros(datashape)
    print(data.shape)
    for i in range(len(population.m1s)):  # loop over init primary masses
        for k in range(len(population.ps)):  # loop over p_init
            for j in range(len(population.qs)):  # loop over q_init
                for m in range(3):  # loop over ET/no ET/single rot
                    hhist1 = population.hists[i, k, j, m, 0]
                    hhist2 = population.hists[i, k, j, m, 1]
                    m1here = population.m1s[i]
                    phere = population.ps[k]
                    qhere = population.qs[j]
                    et_here = m  # 1 for ET, 0 for no ET, 2 for single rot
                    params[i, k, j, m] = [m1here, phere, qhere,
                                          et_here]  # save the initial parameters
                    if hhist1 is None or hhist2 is None:
                        print('none detected', i, k, j, m)
                        data[i, k, j, m] = 0

                        continue

                    if ignore_condition(hhist1):
                        continue  # if here, times can remain zero and will not contribute
                    data[i, k, j, m] = \
                        get_d2t_dx_dy(hhist1, hhist2, xbin_edges, ybin_edges, x_lookup, y_lookup)
    return params, data


def sum_population(data):
    summed_data = np.zeros(data.shape[-2:])
    for i in range(data.shape[0]):  # loop over init primary masses
        for k in range(data.shape[1]):  # loop over p_init
            for j in range(data.shape[2]):  # loop over q_init
                for m in range(data.shape[3]):  # loop over ET
                    summed_data += data[i, k, j, m]
    return summed_data


def normalize(data):
    return data / np.sum(data)


def get_population_summed_weighted_probabilities(population, xbin_edges, ybin_edges,
                                                 x_lookup, y_lookup, ignore_condition):
    """
    does the fetching, weighting, filtering, summing and normalizing of the population data
    see `get_x_y_probabilities` for explanation of arguments
    :return: arrays of the normalized probabilities for the ET, no ET and single rot grids
    """
    params, raw = get_x_y_probabilities(population, xbin_edges, ybin_edges,
                                        x_lookup, y_lookup, ignore_condition)
    weigted, _ = weight_data(raw, params)
    summed_et = sum_population(weigted[:, :, :, [1]])
    summed_no_et = sum_population(weigted[:, :, :, [0]])
    summed_single = sum_population(weigted[:, :, :, [2]])
    return normalize(summed_et), normalize(summed_no_et), normalize(summed_single)


# prior-weighting functions
def q_weight(qq):  # mass ratio distribution
    return qq ** 0  # flat


def m1_weight(mm1, ms):  # IMF
    ix = np.nonzero(ms == mm1)
    assert len(ix) == 1
    m_low, m_high = _get_bounds(ix[0][0], mm1, ms)
    assert m_low < m_high
    return m_low ** (-1.35) - m_high ** (-1.35)  # integrated salpeter


def p_weight(pp):  # period_distribution
    return np.log10(pp) ** 0  # log-flat


def _get_bounds(ii, param, population):
    if ii != 0:
        low = 0.5 * (param + population[ii - 1])  # in middle between i and i-1
    else:
        low = param - 0.5 * (population[ii + 1] - param)  # extrapolate lower edge
    if ii != len(population) - 1:
        high = 0.5 * (param + population[ii + 1])
    else:
        high = param + 0.5 * (param - population[ii - 1])  # extrapolate upper edge
    return low, high


def total_weight(mm1, pp, qq, ms):
    wm = m1_weight(mm1, ms)
    wp = p_weight(pp)
    wq = q_weight(qq)
    return wm * wp * wq  # assumes constant star-formation rate


def weight_data(data_in, population_params):
    ms = np.unique(population_params[:, :, :, :, 0])  # masses live here
    ps = np.unique(population_params[:, :, :, :, 1])  # period live here
    qs = np.unique(population_params[:, :, :, :, 2])  # qs live here
    weightshere = np.zeros(population_params.shape[:-1])
    data_out = np.zeros(data_in.shape)
    for i in range(len(ms)):  # loop over init primary masses
        for k in range(len(ps)):  # loop over p_init
            for j in range(len(qs)):  # loop over q_init
                for m in range(3):  # loop over ET
                    m1here = population_params[i, k, j, m, 0]
                    phere = population_params[i, k, j, m, 1]
                    qhere = population_params[i, k, j, m, 2]
                    w = total_weight(m1here, phere, qhere, ms)
                    weightshere[i, k, j, m] = w
                    data_out[i, k, j, m] = w * data_in[i, k, j, m]
    return data_out, weightshere


# some ignore functions
def not_rlof_or_L2OF_at_ZAMS(hhist):
    return hhist.get('had_RLOF_at_ZAMS')[-1] != 1 or hhist.get('stopping_condition')[-1] == 2.0


def rlof_or_L2OF_at_ZAMS(hhist):
    return hhist.get('had_RLOF_at_ZAMS')[-1] == 1 or hhist.get('stopping_condition')[-1] == 2.0


def L2OF_at_ZAMS(hhist):
    return hhist.get('stopping_condition')[-1] == 2.0


def marginalize_over_x(xbin_edges, ybins, data):
    hist = np.zeros(ybins.shape)
    for j in range(hist.size):
        for i in range(xbin_edges.size - 1):
            hist[j] += data[j, i] * (xbin_edges[i + 1] - xbin_edges[i])
    hist /= np.sum(hist)
    return hist


def marginalize_over_y(xbins, ybin_edges, data):
    hist = np.zeros(xbins.shape)
    for i in range(hist.size):
        for j in range(ybin_edges.size - 1):
            hist[i] += data[j, i] * (ybin_edges[j + 1] - ybin_edges[j])
    hist /= np.sum(hist)
    return hist


def collapse_xbins(data, n):
    new_sh = (data.shape[0], data.shape[1] // n)
    new_data = np.zeros(new_sh)
    for i in range(new_data.shape[1]):
        for j in range(n):
            new_data[:, i] += data[:, n * i + j]
    return new_data


# plotting functions
def make_marg_y_plot(aax, xbins, ybin_edges, data, **kwargs):
    hist = marginalize_over_y(xbins, ybin_edges, data)
    aax.bar(xbins, hist, width=xbins[-1] - xbins[-2], **kwargs)


def make_marg_cumul_y_plot(aax, xbins, ybin_edges, toplot, *args, **kwargs):
    hist = marginalize_over_y(xbins, ybin_edges, toplot)
    cumul = np.cumsum(hist)
    aax.plot(xbins, cumul, *args, **kwargs)


def make_contours(aax, x, y, toplot, regularized_data, hdis, *args, **kwargs):
    # Effectively double integral over dxdy
    density = np.sort(regularized_data.flatten())
    cdf = (density / density.sum()).cumsum()

    # The highest-density interval is equivalent to an inverse empirical (non-normal) CDF
    # evaluated at the lower-bound quantile.
    include_quantile = hdis
    exclude_quantile = 1 - include_quantile
    idensity = cdf.searchsorted(exclude_quantile)
    qdensity = density[idensity]

    print(qdensity, idensity)
    cons = aax.contour(x, y, toplot, levels=np.unique(qdensity), *args, **kwargs)

    fmt = {}
    for ll, s in zip(cons.levels, hdis):
        fmt[ll] = s  # translate density back to hdi for the label

    aax.clabel(cons, cons.levels, fmt=fmt, inline=True, fontsize=5, manual=True)
    return cons


def regularize_y(xbins, ybins, data):
    regular = np.linspace(ybins.min(), ybins.max(), len(ybins))
    qr, pr = np.meshgrid(xbins, regular)
    xis = np.dstack((qr, pr))
    regularized = spinp.interpn(
        points=(xbins, ybins),
        values=data.T,
        xi=xis
    )
    return regularized


def plot_distributions(data: list, xbins, ybins, ybin_edges,
                       need_to_regularize_y=False,
                       cmap=plt.colormaps['inferno'],
                       levels=np.array([0.99, 0.95, 0.90, 0.8]),
                       minprob=1e-5,
                       marg=True,
                       figsize=(6.76, 2.25)):
    x, y = np.meshgrid(xbins, ybins)
    panels = len(data)
    regularized_data = [None] * len(data)
    for i in range(len(data)):
        if need_to_regularize_y:
            regularized_data[i] = regularize_y(xbins, ybins, data[i])
        else:
            regularized_data[i] = data[i]

    fig = plt.figure(figsize=figsize)
    widths = tuple([12] * panels + [1])
    rows = 2 if marg else 1
    heights = (2, 7) if marg else None
    gs = fig.add_gridspec(rows, panels + 1, left=0.11, right=0.89, bottom=0.15, top=0.90,
                          wspace=0.02, hspace=0.05, width_ratios=widths,
                          height_ratios=heights)
    axs = np.empty((rows, panels + 1), dtype=object)
    pq_row = 1 if marg else 0
    for i in range(panels):
        if i != 0:
            sharey = axs[pq_row, 0]
        else:
            sharey = None
        axs[pq_row, i] = fig.add_subplot(gs[pq_row, i], sharex=sharey, sharey=sharey)
        axs[pq_row, i].set_ylim(0.6, 3)
        p1 = axs[pq_row, i].pcolormesh(x, y, np.ma.array(data[i], mask=data[i] <= minprob),
                                       norm=mcolors.LogNorm(vmin=minprob), cmap=cmap)
        # axs[pq_row, i].set_ylim(top=4)
        if levels is not None:
            c1 = make_contours(axs[pq_row, i], x, y, data[i], regularized_data[i], levels,
                               colors='k', linewidths=0.5)
        if i != 0:
            plt.setp(axs[pq_row, i].get_yticklabels(), visible=False)

    axs[pq_row, -1] = fig.add_subplot(gs[pq_row, -1])
    cbar = plt.colorbar(p1, cax=axs[pq_row, -1], use_gridspec=True)

    colors = ['r', 'b', 'm']
    titles = ['ET', 'no ET', 'single rotation']
    for i in range(panels):
        if marg:
            axs[0, i] = fig.add_subplot(gs[0, i], sharey=axs[0, 0], sharex=axs[1, i], zorder=-i)
            make_marg_y_plot(axs[0, i], xbins, ybin_edges, data[i], color=colors[i], edgecolor='k')
            plt.setp(axs[0, i].get_xticklabels(), visible=False)
            if i != 0:
                plt.setp(axs[0, i].get_yticklabels(), visible=False)
        axs[0, i].set_title(titles[i])

    return axs
