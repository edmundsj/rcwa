from matplotlib import pyplot as plt


class Results:
    def __init__(self, results_dict):
        self.inner_dict = results_dict

    def __getitem__(self, key):
        return self.inner_dict[key]

    def keys(self):
        return self.inner_dict.keys()

    def items(self):
        return self.inner_dict.items()

    def values(self):
        return self.inner_dict.values()

    def plot(self, x='wavelength', y='RTot', c=None, fig=None, ax=None, show=False):
        """
        :param x: Variable to plot along the x-axis
        :param y: Variable to plot along the y-axis
        :param c: Variable to plot vs. x/y as distinct curves
        :param fig: Figure to use for plotting. If None, will create with pyplot interface
        :param ax: Axes to use for  plotting. If None, will create with pyplot interface.
        :param show: Whether to show the plot using the pyplot interface. False by default.

        :returns fig, ax: Figure and Axes objects created with matplotlib pyplot interface
        """
        if fig is None and ax is None:
            fig, ax = plt.subplots()
        elif fig is not None and ax is None:
            ax = fig.add_subplot()

        x_data = self[x]
        if hasattr(y, '__iter__') and not isinstance(y, str):
            y_data = [self[yy] for yy in y]
        else:
            y_data = [self[y]]

        for dat in y_data:
            ax.plot(x_data, dat)
        ax.legend(y)
        ax.set_xlabel(x)
        ax.set_ylabel(y)

        if show:
            plt.show()

        return fig, ax

