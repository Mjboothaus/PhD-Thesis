"""
PyOZ - Ornstein-Zernike Equation Solver

This program numerically solves the Ornstein-Zernike equation.
"""

import sys
import numpy as np
from typing import List, Optional
import toml
from dataclasses import dataclass

from input import InputParser
from potential import PotentialCalculator
from dft import DFT
from closure import ClosureCalculator
from solver import Solver
from property import PropertyCalculator
from plot import PyOZPlotter


@dataclass
class Fluid:
    ncomponents: int
    dens: dict
    closure: str
    potential: List[str]
    sigma: List[float]


@dataclass
class Numerics:
    npoints: int
    deltar: float
    deltak: float
    max_iter: int
    convergence_crit: float
    max_dsqn: float
    do_nr: bool
    mix_param: float
    nr_convergence_factor: float
    do_graphics: bool


class Config:
    def __init__(self, config_dict: dict):
        self.fluid = Fluid(**config_dict["system"])
        self.numerics = Numerics(**config_dict["control"])
        self.parameters = config_dict["parameters"]
        self.constants = config_dict["constants"]
        self.output = config_dict["output"]


class PyOZ:
    def __init__(self, config: Config):
        self.fluid = config.fluid
        self.numerics = config.numerics
        self.parameters = config.parameters
        self.constants = config.constants
        self.output = config.output
        self.initialize_system()

    def initialize_system(self):
        self.r = np.array([(x + 1) * self.numerics.deltar for x in range(self.numerics.npoints)])
        self.k = np.array([(x + 1) * self.numerics.deltak for x in range(self.numerics.npoints)])

        print("Initializing DFT routines")
        self.dft = DFT(self.numerics.npoints, self.numerics.deltar, self.numerics.deltak, self.r, self.k)
        self.dft.print_status()
        print("")

        self.potential_calculator = PotentialCalculator(self)
        self.closure_calculator = ClosureCalculator(self.fluid)
        self.solver = Solver(self.fluid.ncomponents)
        self.property_calculator = PropertyCalculator(self)

        if self.numerics.do_graphics:
            self.plotter = PyOZPlotter(self)

    def solve(self):
        U_ij, U_ij_individual, dU_ij_individual, U_discontinuity, U_erf_ij = (
            self.potential_calculator.calculate_potentials(self.dft, self.r, self.k)
        )
        modMayerFunc = self.potential_calculator.def_modMayerFunc(
            U_ij, U_ij_individual, U_discontinuity, U_erf_ij["real"]
        )

        G_r_ij = np.zeros((self.fluid.ncomponents, self.fluid.ncomponents, self.numerics.npoints))
        C_f_ij = np.zeros_like(G_r_ij)
        Cs_f_ij = np.zeros_like(G_r_ij)

        E_ij = np.eye(self.fluid.ncomponents)

        converged = False
        for iteration in range(self.numerics.max_iter):
            G_o_ij = G_r_ij.copy()

            cs_r_ij, g_r_ij = self.closure_calculator.calculate_gamma_term(
                self.r, modMayerFunc, U_discontinuity, G_r_ij
            )

            for i in range(self.fluid.ncomponents):
                for j in range(self.fluid.ncomponents):
                    Cs_f_ij[i, j], C_f_ij[i, j] = self.dft.dfbt(
                        cs_r_ij[i, j], norm=self.fluid.dens["ij"][i, j], corr=-U_erf_ij["fourier"][i, j]
                    )

            H_f_ij = self.solver.solve(
                E_ij - self.dft.ft_convolution_factor * C_f_ij, C_f_ij, self.numerics.npoints
            )

            S = E_ij + H_f_ij
            G_f_ij = S - E_ij - Cs_f_ij

            for i in range(self.fluid.ncomponents):
                for j in range(self.fluid.ncomponents):
                    G_r_ij[i, j] = self.dft.idfbt(
                        G_f_ij[i, j], norm=self.fluid.dens["ij"][i, j], corr=-U_erf_ij["real"][i, j]
                    )

            if self.check_convergence(G_o_ij, G_r_ij):
                converged = True
                break

            if self.numerics.do_graphics:
                self.plotter.update_plot(U_r=U_ij, U_erf=U_erf_ij["real"], G_r=G_r_ij, g_r=g_r_ij)

        if converged:
            print(f"Converged after {iteration + 1} iterations")
        else:
            print("Failed to converge")

        properties = self.property_calculator.calculate_properties(self.r, g_r_ij, S)
        return properties

    def check_convergence(self, G_old: np.ndarray, G_new: np.ndarray) -> bool:
        norm_dsqn = np.linalg.norm(G_new - G_old)
        return norm_dsqn <= self.numerics.convergence_crit


def load_toml_config(file_path: str) -> Config:
    with open(file_path, "r") as f:
        config_dict = toml.load(f)
    return Config(config_dict)


def main(argv: Optional[List[str]] = None):
    print("\npyOZ - iterative solver of the Ornstein-Zernike equation")
    print("Refactor code based on that by:\n")
    print("Lubos Vrbka, 2008-2009\n")

    if argv is None:
        argv = sys.argv[1:]

    input_parser = InputParser()
    cmdline_args = input_parser.parse_cmdline(argv)

    if cmdline_args.get("config"):
        # Use TOML configuration file
        config = load_toml_config(cmdline_args["config"])
    else:
        # Use command line arguments
        config_dict = input_parser.parse_input(cmdline_args)
        config = Config(config_dict)

    solver = PyOZ(config)
    properties = solver.solve()

    print("Calculation completed. Results:")
    for prop, value in properties.items():
        print(f"{prop}: {value}")

    if solver.numerics.do_graphics:
        solver.plotter.show()


if __name__ == "__main__":
    main()
