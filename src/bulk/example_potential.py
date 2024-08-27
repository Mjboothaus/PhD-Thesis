import numpy as np
import matplotlib.pyplot as plt
from pyoz_potential import PotentialCalculator

# Example configuration and parameters
config = {
    "npoints": 1000,
    "alpha": 0.5,  # For Ewald summation in Coulomb potential
}

system = {"ncomponents": 2, "potential": ["hs", "lj", "coulomb"]}

parameters = {
    "sigma": [1.0, 1.2],  # Diameter for each component
    "epsilon": [1.0, 1.5],  # LJ well depth for each component
    "charge": [1.0, -1.0],  # Charges for each component
}

constants = type("Constants", (), {"beta": 1.0})()  # Simple object to hold beta value

# Create distance arrays
r = np.linspace(0.1, 10, config["npoints"])
k = np.linspace(0.1, 10, config["npoints"])


# Create a mock DFT object (replace with actual DFT implementation if available)
class MockDFT:
    def __init__(self):
        pass

    def dfbt(self, func, norm, corr):
        return func, func  # Mock implementation


dft = MockDFT()

# Create PotentialCalculator instance
calculator = PotentialCalculator(config, system, parameters, constants)

# Calculate potentials
U_ij, U_ij_individual, dU_ij_individual, U_discontinuity, U_erf_ij = calculator.calculate_potentials(
    dft, r, k
)

# Calculate modified Mayer functions
modMayerFunc = calculator.def_modMayerFunc(U_ij, U_ij_individual, U_discontinuity, U_erf_ij["real"])

# Plotting
plt.figure(figsize=(15, 10))

# Plot total potential
plt.subplot(2, 2, 1)
plt.title("Total Potential")
for i in range(system["ncomponents"]):
    for j in range(i, system["ncomponents"]):
        plt.plot(r, U_ij[i, j], label=f"U_{i}{j}")
plt.xlabel("r")
plt.ylabel("U(r)")
plt.legend()

# Plot individual potentials
plt.subplot(2, 2, 2)
plt.title("Individual Potentials")
for pot_type in U_ij_individual:
    plt.plot(r, U_ij_individual[pot_type][0, 1], label=pot_type)
plt.xlabel("r")
plt.ylabel("U(r)")
plt.legend()

# Plot Erf-corrected potential
plt.subplot(2, 2, 3)
plt.title("Erf-corrected Potential")
plt.plot(r, U_erf_ij["real"][0, 1], label="Real")
plt.plot(k, U_erf_ij["fourier"][0, 1], label="Fourier")
plt.xlabel("r / k")
plt.ylabel("U_erf(r) / U_erf(k)")
plt.legend()

# Plot modified Mayer functions
plt.subplot(2, 2, 4)
plt.title("Modified Mayer Functions")
for key in modMayerFunc:
    plt.plot(r, modMayerFunc[key][0, 1], label=key)
plt.xlabel("r")
plt.ylabel("f(r)")
plt.legend()

plt.tight_layout()
plt.show()
