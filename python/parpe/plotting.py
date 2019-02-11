import matplotlib.pyplot as plt
import numpy as np

def plotCostTrajectory(costTrajectory,
                       color=None,
                       scaleToIteration=1,
                       legend=True,
                       legend_loc='upper right',
                       ax=None, log=True):
    """Plot the provided cost trajectory

    Arguments:
    costTrajectory: ndarray(numIterations x numStarts): cost over iterations for all optimizations

    Returns:

    """

    if ax is None:
        ax = plt.subplots()[1]

    for start in range(costTrajectory.shape[1]):
        y = costTrajectory[:,start]
        y = y[~np.isnan(y)]
        c = color[start] if color and len(color) > 1 else color
        if log:
            ax.semilogy(y, c=c)
        else:
            ax.plot(y, c=c)

    # set viewport

    minCost = np.nanmin(costTrajectory)
    maxCost = np.nanmax(costTrajectory[scaleToIteration:,:])
    #ax.set_ylim(bottom=(1 - np.sign(minCost) * 0.02) * minCost,
    #            top=(1 + np.sign(maxCost) * 0.02) * maxCost)
    ax.set_ylim(bottom=minCost - 0.01 * np.abs(maxCost - minCost),
                top=maxCost + 0.01 * np.abs(maxCost - minCost))

    if legend:
        ax.legend(['start %d'%i for i in range(costTrajectory.shape[1]) ], loc=legend_loc)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Cost')

    return ax


def plotDoseResponseLogDose(conc, mes, sim, title, ax):
    ax.semilogx(conc, mes, '.', label='mes')
    ax.semilogx(conc, sim, '.', label='sim')
    ax.set_ylabel('Proliferation')
    ax.set_xlabel('[Drug]')
    ax.set_title(title)
    #ax.legend();


def plotDoseResponseCategorical(conc, mes, sim, title, ax):
    """
    e.g. #
    fig, ax = plt.subplots()
    plotDoseResponse(np.array([10, 3, 5]), np.array([1, 2, 3]), np.array([1, 3, 2]), 'test', ax)
    """
    order = np.argsort(conc)
    concStr = [str(conc[i]) for i in order]

    ax.scatter(concStr, np.array(mes)[order], label='mes')
    ax.scatter(concStr, np.array(sim)[order], label='sim')
    ax.vlines(concStr, np.array(mes)[order], np.array(sim)[order])

    ax.set_ylabel('Proliferation')
    ax.set_xlabel('[Drug]')
    ax.set_title(title)
    #ax.legend();


def plotCorrelations(ymes, ysim):
    """
    Plot correlation of measured and simulated data

    Arguments:
    ----------
    ymes: @type numpy.ndarray
        measured values n_condition x nt x ny
    ysim: @type numpy.ndarray
        simulated values n_condition x nt x ny
    """
    for iy in range(ysim.shape[2]):
        plotCorrelation(ymes[:, :, iy], ysim[:, :, iy],
                        title='Observable %d' % iy)


def plotCorrelation(ymes, ysim, title=None, alpha=1.0, legend=False, square=True):
    """
    Plot correlation of measured and simulated data

    Arguments:
    ----------
    ymes: @type numpy.ndarray
        measured values n_condition x nt
    ysim: @type numpy.ndarray
        simulated values n_condition x nt
    """
    fig, ax = plt.subplots()
    for icondition in range(ysim.shape[0]):
        x = ymes[icondition, :]
        y = ysim[icondition, :]
        #x, y = flatten_filter_nan(x, y)
        r = correlation_coefficient(x, y)
        ax.scatter(x, y, label='Condition %d, r=%.3f' % (icondition, r), alpha=alpha)

    ax.set_xlabel('measurement (AU)')
    ax.set_ylabel('simulation (AU)')

    if square:
        square_plot_equal_ranges(ax)

    if title:
        plt.title(title)
    if legend:
        plt.legend()

    return ax


def square_plot_equal_ranges(ax, lim=None):
    """Square plot with equal range"""

    #ax.set_aspect('equal', adjustable='box')

    ax.axis('square')

    if lim is None:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        lim = [np.min([xlim[0], ylim[0]]),
               np.max([xlim[1], ylim[1]])]

    ax.set_xlim(lim)
    ax.set_ylim(lim)

    # Same tick mark on x and y
    ax.yaxis.set_major_locator(ax.xaxis.get_major_locator())

    return ax


def plotTrajectoryFits(ymes, ysim, timepoints):
    """For each simulation condition create a plot with time-course for measured and simulated values for all observables.

    Arguments:
    ----------
    ymes: @type numpy.ndarray
        measured values n_condition x nt x ny
    ysim: @type numpy.ndarray
        simulated values n_condition x nt x ny

    """
    for icondition in range(ysim.shape[0]):
        plotTrajectoryFit(ymes[icondition].T,
                          ysim[icondition].T,
                          timepoints,
                          title='Condition %d' % icondition)
        plt.show()


def plotTrajectoryFit(ymes, ysim, timepoints, title=None):
    """Create a plot with time-course for measured and simulated values for all observables.

    Arguments:
    ----------
    ymes: @type numpy.ndarray
        measured values n_observable x nt
    ysim: @type numpy.ndarray
        simulated values n_observable x nt

    """
    fig, ax = plt.subplots()
    for iy in range(ysim.shape[0]):
        ax.plot(timepoints, ysim[iy],
                label='$y_%d$ sim' % (iy), alpha=0.7, c='C%d' % iy)
        ax.plot(timepoints, ymes[iy],
                label='$y_%d$ mes' % (iy), linestyle='dotted', marker='o',
                c='C%d' % iy)
    plt.xlabel('$t$ (s)')
    plt.ylabel('$y_i(t)$ (AU)')
    if title:
        plt.title(title)
    plt.legend()

    return ax


def flatten_filter_nan(a, b):
    """Flatten arrays a and b removing entries for which either a or b is NaN"""

    assert (a.shape == b.shape)

    a = a.flatten()
    b = b.flatten()

    mask = np.isfinite(a + b)
    a = a[mask]
    b = b[mask]

    return a, b


def correlation_coefficient(a, b):
    """Get correlation coefficient for flattened a and b"""

    if not a.size or not b.size:
        return np.nan

    a, b = flatten_filter_nan(a, b)

    return np.corrcoef(a, b)[0, 1]


def plotCorrelationDensity(ymes, ysim, title=None,
                           title_append_corr=True,
                           normalize_percentile=100,
                           nbins=50,
                           contour=True):
    """
    Plot correlation of measured and simulated data

    Arguments:
    ----------
    ymes: @type numpy.ndarray
        measured values n_condition x nt
    ysim: @type numpy.ndarray
        simulated values n_condition x nt
    """
    assert (ymes.shape == ysim.shape)
    from scipy.stats import kde

    x, y = flatten_filter_nan(ymes, ysim)

    fig, ax = plt.subplots()

    # set viewport
    ymin = np.nanmin([x, y])
    ymax = np.nanmax([x, y])

    #ax.hexbin(x, y, gridsize=nbins, cmap=plt.cm.BuGn_r)

    # normalize colors
    vmax = np.percentile([x, y], normalize_percentile)

    # get kernel density estimate
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[ymin:ymax:nbins * 1j, ymin:ymax:nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    ax.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud',
                  cmap=None, vmax=vmax)  # plt.cm.BuGn_r

    if contour:
        ax.contour(xi, yi, zi.reshape(xi.shape))

    square_plot_equal_ranges(ax, lim=[ymin, ymax])

    ax.set_xlabel('measurement (AU)')
    ax.set_ylabel('simulation (AU)')

    r = correlation_coefficient(x, y)
    if title is None:
        title = 'r=%.3f' % (r)
    elif title_append_corr:
        title += ' (r=%.3f)' % (r)
    ax.set_title(title)

    return ax


def plotCorrelationsByObsDensity(ymes, ysim):
    """
    Plot correlation of measured and simulated data

    Arguments:
    ----------
    ymes: @type numpy.ndarray
        measured values n_condition x nt x ny
    ysim: @type numpy.ndarray
        simulated values n_condition x nt x ny
    """
    for iy in range(ysim.shape[2]):
        plotCorrelationDensity(ymes[:, :, iy], ysim[:, :, iy],
                        title='Observable %d' % iy)


def plotWaterfall(finalCost):
    """Plot "waterfall plot"

    Sorted scatter plot of optimization results.

    Arguments:
    finalCost: ndarray(numStarts) of final cost

    Returns:

    """

    # x axis should have integer labels
    from matplotlib.ticker import MaxNLocator
    ax = plt.figure().gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    order = np.argsort(finalCost)
    plt.scatter(range(finalCost.size), finalCost[:, order])

    plt.xlabel('Sorted start index')
    plt.ylabel('Final cost')


def plotCorrelationBox(data):
    """Create boxplot of correlation by observable

    Arguments:
        data: Correlations as obtained from getCorrTable
    """

    mask = ~np.isnan(data)
    filtered_data = [d[m] for d, m in zip(data.T, mask.T)]

    plt.figure(figsize=(15, 4))
    plt.boxplot(filtered_data)


def plotCorrelationBoxMulti(datasets, labels,
                            legend_loc='lower center'):
    """
    Plot correlation boxplots of multi-start results for multiple datasets side by side

    data1: nstarts x ny
    """
    for x in datasets:
        assert (datasets[0].shape[1] == x.shape[1])

    num_groups = len(datasets)
    num_observations = datasets[0].shape[1]

    filtered_data = []
    for x in datasets:
        mask = ~np.isnan(x)
        filtered = [d[m] for d, m in zip(x.T, mask.T)]
        filtered_data.append(filtered)

    def set_box_color(bp, color):
        '''
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)
        '''
        plt.setp(bp['medians'], color='b')
        for patch in bp['boxes']:
            patch.set_facecolor(color)

    fig = plt.figure(figsize=(24, 8))
    ax = fig.gca()

    for j, x in enumerate(filtered_data):
        bp = ax.boxplot(x,
                        positions=[i * num_groups + j - (num_groups / 2 - 0.5)
                                   for i in range(num_observations)],
                        patch_artist=True)
        col = 'C' + str(j)
        set_box_color(bp, col)
        ax.plot([], c=col, label=labels[j])

    ax.set_xlabel('Model output index')
    ax.set_ylabel('Correlation coefficient')

    ax.legend(loc=legend_loc)

    ticks = ['%d' % i for i in range(num_observations)]
    ax.set_xticks(range(0, len(ticks) * num_groups, num_groups))
    ax.set_xticklabels(ticks)
    ax.set_xlim(-2, num_observations * num_groups - 2)

    # from matplotlib.ticker import MaxNLocator
    # ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    return ax
