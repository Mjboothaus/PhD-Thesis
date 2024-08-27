# File: pyoz_plot.py

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from typing import Dict, Any, Optional


class PyOZPlotter:
    def __init__(self, config: Dict[str, Any], system: Dict[str, Any], constants: Any, r: np.ndarray):
        self.config = config
        self.system = system
        self.constants = constants
        self.r = r
        self.fig, self.axes = plt.subplots(2, 2, figsize=(12, 10))
        self.fig.suptitle("PyOZ Solver Results")
        self.initialize_plots()

    def initialize_plots(self):
        self.plots = {
            "U_r": [
                ax.plot([], [], label=f"U_{i}{j}")[0]
                for ax in self.axes.flat[:2]
                for i in range(self.system["ncomponents"])
                for j in range(i, self.system["ncomponents"])
            ],
            "G_r": [
                ax.plot([], [], label=f"G_{i}{j}")[0]
                for ax in self.axes.flat[2:]
                for i in range(self.system["ncomponents"])
                for j in range(i, self.system["ncomponents"])
            ],
        }

        for ax in self.axes.flat:
            ax.set_xlabel("r")
            ax.legend()

        self.axes[0, 0].set_ylabel("U(r)")
        self.axes[0, 1].set_ylabel("U_erf(r)")
        self.axes[1, 0].set_ylabel("G(r)")
        self.axes[1, 1].set_ylabel("g(r)")

    def update_plot(
        self,
        U_r: Optional[np.ndarray] = None,
        U_erf: Optional[np.ndarray] = None,
        G_r: Optional[np.ndarray] = None,
        g_r: Optional[np.ndarray] = None,
    ):
        data = [U_r, U_erf, G_r, g_r]
        for (i, ax), d in zip(enumerate(self.axes.flat), data):
            if d is not None:
                for j, line in enumerate(ax.get_lines()):
                    comp1, comp2 = divmod(j, self.system["ncomponents"])
                    line.set_data(self.r, d[comp1, comp2])
            ax.relim()
            ax.autoscale_view()

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def animate(self, frame_data):
        self.update_plot(*frame_data)

    def show(self):
        plt.show()

    def save(self, filename: str):
        self.fig.savefig(filename)


# Usage:
# plotter = PyOZPlotter(config, system, constants, r)
# plotter.update_plot(U_r=U_ij, U_erf=U_erf_ij['real'], G_r=G_r_ij, g_r=g_r_ij)
# plotter.show()  # or plotter.save('pyoz_results.png')
