## Ornstein-Zernike Equation Solver

This code implements a numerical solver for the Ornstein-Zernike (OZ) equation, which is fundamental in the statistical mechanical theory of liquids. In this case it is for bulk fluids (hence the sub-directory name.)

### Key Equations

#### Total Interaction Potential

   The total pair potential is the sum of individual contributions:

   $$U_{ij}(r) = U_{hs}(r) + U_{lj}(r) + U_{coulomb}(r)$$

#### Modified Mayer Functions

   $$f_{ij}(r) = e^{-\beta U_{ij}(r)} - 1$$

#### Ornstein-Zernike (OZ) Equation

   In Fourier space:

   $$h_{ij}(k) = c_{ij}(k) + \sum_k \rho_k c_{ik}(k) h_{kj}(k)$$

   Where $h_{ij}(r)$ is the total correlation function, $c_{ij}(r)$ is the direct correlation function, and $\rho_k$ is the number density of component $k$.

#### Closure Relations
   
1. The Hypernetted Chain (HNC) closure:

   $$c_{ij}(r) = \exp[-\beta U_{ij}(r) + \gamma_{ij}(r)] - \gamma_{ij}(r) - 1$$

   where $\gamma_{ij}(r) = h_{ij}(r) - c_{ij}(r)$ is the indirect correlation function.

2. Percus-Yevick (PY) closure:
    $$c_{ij}(r)={ij}(r)[1+γ{ij}(r)]{ij}(r)=f{ij}(r)[1+γ_{ij}(r)]$$

    where $f_{ij}(r) = e^{-\beta U_{ij}(r)} - 1$ is the Mayer function.

3. Mean Spherical Approximation (MSA) closure:

    $$-\beta U_{ij}(r)&\text{for }r\geq \sigma_{ij}\\-1&\text{for }r<\sigma_{ij}\end{cases}$$
    
    where $\sigma_{ij}$ is the contact distance between particles of type $i$ and $j$.

#### Ng Correction

   For long-range potentials:

   $$c_{ij}(r) = c^s_{ij}(r) - \beta U^{long}_{ij}(r)$$

   Where $c^s_{ij}(r)$ is the short-range direct correlation function.

#### Structure Factor

   $$S_{ij}(k) = \delta_{ij} + \sqrt{\rho_i \rho_j} h_{ij}(k)$$


### Numerical Methods

#### Newton-Raphson/Conjugate Gradient Method

Solves the linear system $AX = B$ where:

   - $A$ is a linear operator
   - $X$ is $\delta\gamma$
   - $B$ is the difference between input and output $\gamma$


#### Discrete Fourier Transforms

   Used to switch between real and reciprocal space.

#### Picard Iteration

   $$\gamma^{new}_{ij}(r) = (1-\alpha)\gamma^{old}_{ij}(r) + \alpha\gamma^{calculated}_{ij}(r)$$

   Where $\alpha$ is the mixing parameter.

#### Convergence Criterion 

   Based on the relative change in $\gamma_{ij}(r)$ between iterations:

   $$DSQN = \frac{\sum_{i,j} \int [\gamma^{new}_{ij}(r) - \gamma^{old}_{ij}(r)]^2 dr}{\sum_{i,j} \int [\gamma^{old}_{ij}(r)]^2 dr}$$

This solver provides a powerful tool for studying the structure of complex fluid systems, allowing for the calculation of various thermodynamic and structural properties based on the obtained correlation functions.