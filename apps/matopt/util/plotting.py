import matplotlib.pyplot as plt


def makeBondLengthsHistogram(BondLengths, bins=None, visible=False, save=False, filename=None):
    """

    Args:
        BondLengths: param bins:  (Default value = None)
        visible: Default value = False)
        save: Default value = False)
        filename: Default value = None)
        bins:  (Default value = None)

    Returns:

    """
    plt.hist(BondLengths, bins=bins)
    if (visible):
        plt.show()
    if (save):
        plt.savefig(filename)
    plt.close()
