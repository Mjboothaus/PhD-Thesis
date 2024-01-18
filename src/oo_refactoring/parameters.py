# src/parameters.py
import toml
from dataclasses import dataclass
import numpy as np


@dataclass
class Fluid:
    name: str
    symbol: str
    component: list[str]
    valence: np.ndarray
    charge: np.ndarray
    temperature: float
    concentration: np.ndarray
    epsilon_r: float
    n_component: int
    n_pair: int
    index: int
    charge_pair: np.ndarray
    rho: np.ndarray
    beta: float
    epsilon: float
    cr_filename: str = ""


class FluidInitialiser:
    def __init__(self, config_path):
        self.config_path = config_path
        self.fluids = self.load_fluids()

    def load_fluids(self):
        with open(self.config_path, "r") as f:
            try:
                config = toml.load(f)
            except Exception:
                raise FileNotFoundError(f"{self.config_path} not found.")

        fluids = {}
        for fluid_name, attributes in config.items():
            # Calculate derived values here
            valence = np.array(attributes.get("valence", []))
            temperature = attributes.get("temperature", 0.0)
            concentration = np.array(attributes.get("concentration", []))
            epsilon_r = attributes.get("epsilon_r", 0.0)
            beta = 1.0 / (Boltzmann * temperature)
            epsilon = 4.0 * np.pi * epsilon_r * epsilon_0
            rho = concentration * Avogadro / 1.0e27
            charge = valence * elementary_charge
            n_component = len(attributes.get("component", []))
            n_pair = int((n_component + 1) * n_component / 2)
            charge_pair = np.zeros(n_pair)  # Placeholder for actual calculation

            fluids[fluid_name] = Fluid(
                name=attributes.get("name", ""),
                symbol=attributes.get("symbol", ""),
                component=attributes.get("component", []),
                valence=np.array(attributes.get("valence", [])),
                temperature=attributes.get("temperature", 0.0),
                concentration=np.array(attributes.get("concentration", [])),
                epsilon_r=attributes.get("epsilon_r", 0.0),
                n_component=len(attributes.get("component", [])),
                n_pair=int(
                    (len(attributes.get("component", [])) + 1)
                    * len(attributes.get("component", []))
                    / 2
                ),
                index=attributes.get("index", 0),
                cr_filename=attributes.get("cr_filename", ""),
                beta=beta,
                epsilon=epsilon,
                rho=rho,
                charge=charge,
                charge_pair=charge_pair,
            )
            # Write derived values back to the TOML file for reference
            self.write_derived_values(
                fluid_name, beta, epsilon, rho, charge, charge_pair
            )

        return fluids

    def write_derived_values(self, fluid_name, beta, epsilon, rho, charge, charge_pair):
        # This method writes the derived values back to the TOML file
        # with comments indicating they are for information only
        derived_data = {
            "beta": beta,
            "epsilon": epsilon,
            "rho": rho.tolist(),
            "charge": charge.tolist(),
            "charge_pair": charge_pair.tolist(),
        }
        with open(self.config_path, "a") as f:
            f.write(f"\n# Derived values for {fluid_name} - for info only\n")
            toml.dump({fluid_name: derived_data}, f)

    def get_fluid(self, symbol):
        # Normalize the symbol to match the TOML configuration
        symbol = symbol.lower()
        return self.fluids.get(symbol)


# Usage example:
# fluid_initialiser = FluidInitialiser(config_path="path/to/fluids.toml")
# lennard_jones_fluid = fluid_initialiser.get_fluid("lj")
