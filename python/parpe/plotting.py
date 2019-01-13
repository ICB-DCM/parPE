import matplotlib.pyplot as plt
import numpy as np

def plotCostTrajectory(costTrajectory,
                       color=None,
                       scaleToIteration=1,
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
    ax.set_ylim(bottom=(1 - np.sign(minCost) * 0.02) * minCost,
                top=(1 + np.sign(maxCost) * 0.02) * maxCost)

    ax.legend(['start %d'%i for i in range(costTrajectory.shape[1]) ], loc=legend_loc)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Cost')

    return ax


def plotDoseResponseLogDose(conc, mes, sim, title, ax):
    ax.semilogx(conc, mes, '.', label='mes')
    ax.semilogx(conc, sim, '.', label='sim')
    ax.set_ylabel('Proliferation')
    ax.set_xlabel('[Drug]');
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
    ax.set_xlabel('[Drug]');
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


def plotCorrelation(ymes, ysim, title=None):
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
        r = correlation_coefficient(x, y)
        ax.scatter(x, y, label='Condition %d, r=%.3f' % (icondition, r))

    ax.set_xlabel('measurement (AU)')
    ax.set_ylabel('simulation (AU)')

    # Square plot with equal range
    ax.axis('square')
    # ax.set_aspect('equal', adjustable='box')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    lim = [np.min([xlim[0], ylim[0]]),
           np.min([xlim[1], ylim[1]])]
    ax.set_xlim(lim)
    ax.set_ylim(lim)

    if title:
        plt.title(title)

    plt.legend()


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


def correlation_coefficient(a, b):
    """Get correlation coefficient for flattened a and b"""

    assert (a.shape == b.shape)

    a = a.flatten()
    b = b.flatten()

    mask = np.isfinite(a + b)
    a = a[mask]
    b = b[mask]

    if not a.size or not b.size:
        return np.nan

    return np.corrcoef(a, b)[0, 1]
