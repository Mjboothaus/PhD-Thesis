from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import optimize


@dataclass
class DiscretisationScheme:
    n_point: int
    n_component: int
    n_pair: int
    z_cutoff: float
    grid_size: float
    tolerance: float
    max_iteration: int


class NumericsManager:
    self.n_point = n_point
    self.z_cutoff = z_cutoff
    self.n_component = n_component
    self.n_pair = n_pair
    self.tolerance = tolerance
    self.max_iteration = max_iteration
    self.grid_size = self.z_cutoff / (self.n_point - 1)
    self.discretisation = DiscretisationScheme(
        self.n_point,
        self.n_component,
        self.n_pair,
        self.z_cutoff,
        self.grid_size,
        self.tolerance,
        self.max_iteration,
    )

    @staticmethod
    def perform_calculation(model, parameters):
        pass

    @staticmethod
    def calculate_integrals():
        pass


class Solver:
    @staticmethod
    def solve_model(
        opt_func, tw_initial, fluid, model, discrete, beta_phiw, beta_psi_charge
    ):
        charge_pair = fluid.charge_pair
        n_component = fluid.n_component
        rho = fluid.rho
        f1 = model.f1
        f2 = model.f2
        z = model.z
        z_index = model.z_index
        n_point = discrete.n_point
        tolerance = discrete.tolerance
        max_iteration = discrete.max_iteration
        tw_args = (
            beta_phiw,
            beta_psi_charge,
            charge_pair,
            rho,
            f1,
            f2,
            z,
            n_component,
            n_point,
            z_index,
        )

        return optim.root(
            opt_func,
            tw_initial,
            args=tw_args,
            method="krylov",
            jac=None,
            callback=None,
            options={"disp": True, "maxiter": max_iteration, "fatol": tolerance},
        )


# # src/oo_refactoring/numerics.py
# import numpy as np
# from scipy.integrate import solve_ivp
# from scipy.optimize import minimize

# class NumericsManager:
#     @staticmethod
#     def perform_calculation(model, parameters):
#         # Example of solving a differential equation
#         # Define the differential equation as a function
#         def differential_equation(t, y, args):
#             # Placeholder for the actual differential equation
#             # dy/dt = f(t, y, args)
#             return np.sin(t) - y * args

#         # Initial conditions and time span
#         y0 = np.array([parameters.initial_condition])
#         t_span = (0, parameters.time_end)

#         # Solve the differential equation
#         sol = solve_ivp(differential_equation, t_span, y0, args=(parameters.some_parameter,))

#         # Example of optimization
#         # Define the function to minimize
#         def function_to_optimize(x, args):
#             # Placeholder for the actual function to minimize
#             return (x - args) ** 2

#         # Perform the optimization
#         optimization_result = minimize(function_to_optimize, x0=np.array([parameters.initial_guess]), args=(parameters.some_parameter,))

#         # Example of integration
#         # Define the integrand function
#         def integrand(x):
#             # Placeholder for the actual integrand
#             return np.exp(-x ** 2)

#         # Perform the integration
#         integral_result = np.trapz(integrand(model.z), model.z)

#         # Collect all results in a dictionary
#         results = {
#             'differential_solution': sol.y,
#             'optimization_result': optimization_result.x,
#             'integral_result': integral_result,
#         }

#         return results
