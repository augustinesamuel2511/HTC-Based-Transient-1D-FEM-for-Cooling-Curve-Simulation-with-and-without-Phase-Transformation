import matplotlib.pyplot as plt
import numpy as np


class GraphPlotter:
    def __init__(self, title="2D Plot", xlabel="X-axis", ylabel="Y-axis"):
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.fig, self.ax = plt.subplots()
        self.plots = []  # store plotted data

    def plot_single(self, x, y, label="Data 1", **kwargs):
        """Plot single x and y data."""
        self.ax.plot(x, y, label=label, **kwargs)
        self.plots.append((x, y, label))

    def plot_single_multi_y(self, x, y_list, labels=None, **kwargs):
        """Plot single x with multiple y datasets."""
        for i, y in enumerate(y_list):
            label = labels[i] if labels and i < len(labels) else f"Data {i + 1}"
            self.ax.plot(x, y, label=label, **kwargs)
            self.plots.append((x, y, label))

    def plot_multi_xy(self, x_list, y_list, labels=None, **kwargs):
        """Plot multiple x and y datasets."""
        for i, (x, y) in enumerate(zip(x_list, y_list)):
            label = labels[i] if labels and i < len(labels) else f"Data {i + 1}"
            self.ax.plot(x, y, label=label, **kwargs)
            self.plots.append((x, y, label))

    def show(self, legend=True, grid=True):
        """Show the plot with legend and grid."""
        self.ax.set_title(self.title)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        if legend:
            self.ax.legend()
        if grid:
            self.ax.grid(True)
        plt.show()

    def get_data(self):
        """Return plotted data as list of (x, y, label)."""
        return self.plots

    def export_to_txt(self, filename="graph_data.txt"):
        """Export plotted data to a text file."""
        with open(filename, "w") as f:
            for x_vals, y_vals, lbl in self.plots:
                f.write(f"Label: {lbl}\n")
                f.write("X\tY\n")
                for xv, yv in zip(x_vals, y_vals):
                    f.write(f"{xv:.6f}\t{yv:.6f}\n")
                f.write("\n")
        print(f"Data exported to {filename}")

