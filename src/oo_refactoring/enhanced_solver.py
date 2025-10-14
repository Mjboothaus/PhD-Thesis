"""
Enhanced solver with better optimization tracking and visualization.
"""

import numpy as np
from scipy import optimize
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
import logging

@dataclass
class OptimizationHistory:
    """Tracks optimization progress and metrics."""
    iterations: List[int] = field(default_factory=list)
    function_values: List[float] = field(default_factory=list)
    residuals: List[float] = field(default_factory=list)
    step_sizes: List[float] = field(default_factory=list)
    convergence_metrics: Dict[str, List[float]] = field(default_factory=lambda: {
        'max_residual': [],
        'mean_residual': [],
        'norm_residual': []
    })
    
    def update(self, iteration: int, x: np.ndarray, fun: np.ndarray, step: float):
        """Update history with new iteration data."""
        self.iterations.append(iteration)
        self.function_values.append(float(np.mean(np.abs(fun))))
        residual = float(np.linalg.norm(fun))
        self.residuals.append(residual)
        self.step_sizes.append(step)
        
        # Update convergence metrics
        self.convergence_metrics['max_residual'].append(float(np.max(np.abs(fun))))
        self.convergence_metrics['mean_residual'].append(float(np.mean(np.abs(fun))))
        self.convergence_metrics['norm_residual'].append(residual)

    def get_convergence_summary(self) -> Dict[str, float]:
        """Get summary of final convergence state."""
        if not self.iterations:
            return {}
        
        return {
            'total_iterations': self.iterations[-1],
            'final_residual': self.residuals[-1],
            'final_function_value': self.function_values[-1],
            'final_step_size': self.step_sizes[-1],
            'max_residual': self.convergence_metrics['max_residual'][-1],
            'mean_residual': self.convergence_metrics['mean_residual'][-1],
            'norm_residual': self.convergence_metrics['norm_residual'][-1]
        }

    def plot_convergence(self, display: bool = True) -> Optional[go.Figure]:
        """Create interactive convergence plots."""
        if not self.iterations:
            return None
            
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                "Residuals",
                "Function Values",
                "Step Sizes",
                "All Metrics (Log Scale)"
            )
        )
        
        # Residuals plot
        for metric_name, values in self.convergence_metrics.items():
            fig.add_trace(
                go.Scatter(
                    x=self.iterations,
                    y=values,
                    name=metric_name.replace('_', ' ').title(),
                    mode='lines+markers'
                ),
                row=1, col=1
            )
        
        # Function values plot
        fig.add_trace(
            go.Scatter(
                x=self.iterations,
                y=self.function_values,
                name='Function Value',
                mode='lines+markers'
            ),
            row=1, col=2
        )
        
        # Step sizes plot
        fig.add_trace(
            go.Scatter(
                x=self.iterations,
                y=self.step_sizes,
                name='Step Size',
                mode='lines+markers'
            ),
            row=2, col=1
        )
        
        # Combined log plot
        fig.add_trace(
            go.Scatter(
                x=self.iterations,
                y=self.residuals,
                name='Residual (Log)',
                mode='lines'
            ),
            row=2, col=2
        )
        fig.add_trace(
            go.Scatter(
                x=self.iterations,
                y=self.function_values,
                name='Function Value (Log)',
                mode='lines'
            ),
            row=2, col=2
        )
        fig.add_trace(
            go.Scatter(
                x=self.iterations,
                y=self.step_sizes,
                name='Step Size (Log)',
                mode='lines'
            ),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            height=800,
            showlegend=True,
            title_text="Optimization Convergence"
        )
        
        # Set log scale for combined plot
        fig.update_yaxes(type="log", row=2, col=2)
        
        if display:
            st.plotly_chart(fig, use_container_width=True)
        
        return fig

class EnhancedSolver:
    """Enhanced solver with better optimization tracking and visualization."""
    
    def __init__(self):
        self.history = OptimizationHistory()
        self.iteration = 0
        self.last_x = None
        
    def callback(self, x: np.ndarray, *args, **kwargs) -> None:
        """Callback function to track optimization progress."""
        if self.last_x is None:
            step = 0.0
        else:
            step = float(np.linalg.norm(x - self.last_x))
        
        # Get function value at current point
        fun = self.opt_func(x, *self.opt_args)
        
        # Update history
        self.history.update(self.iteration, x, fun, step)
        
        # Update state
        self.iteration += 1
        self.last_x = x.copy()
        
    def solve_model(
        self,
        opt_func,
        tw_initial,
        fluid,
        model,
        discrete,
        beta_phiw,
        beta_psi_charge,
        *args,
        **kwargs
    ):
        """Solve model with enhanced tracking and visualization."""
        self.opt_func = opt_func
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
        
        self.opt_args = (
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

        # Reset state
        self.iteration = 0
        self.last_x = None
        self.history = OptimizationHistory()
        
        # Run optimization
        result = optimize.root(
            opt_func,
            tw_initial,
            args=self.opt_args,
            method="krylov",
            jac=None,
            callback=self.callback,
            options={
                "disp": True,
                "maxiter": max_iteration,
                "fatol": tolerance
            },
        )
        
        # Add convergence information to result
        result.convergence_history = self.history
        result.convergence_summary = self.history.get_convergence_summary()
        
        return result