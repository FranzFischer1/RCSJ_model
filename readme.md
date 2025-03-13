Below is a **sample README** that provides a detailed derivation and explanation of an RCSJ (Resistively and Capacitively Shunted Josephson Junction) simulation in Python. You can adapt its structure and content as needed for your GitHub project.

---

# RCSJ Simulation in Python

This repository contains a Python implementation of the **RCSJ (Resistively and Capacitively Shunted Josephson Junction)** model, including thermal noise and an optional AC drive. The simulation demonstrates how to numerically integrate the stochastic differential equation (SDE) for the junction's phase over time, using a Heun‐type (explicit trapezoidal) stochastic integrator. Optionally, it uses [Numba](https://numba.pydata.org/) to JIT-compile the core integration loop for speed.

## 1. Introduction and Background

A Josephson junction can be idealized by two superconducting electrodes, separated by a thin insulating barrier that permits tunneling of Cooper pairs. In the **RCSJ model**, a junction of critical current \(I_c\) is shunted by a resistor \(R\) and a capacitor \(C\), all in parallel. An external current \(I_{\text{bias}}(t)\) drives the circuit, which may include both DC and AC components. Additionally, **thermal noise** from the resistor is modeled as a stochastic term coupling to the junction’s voltage. 

### 1.1 Physical Equation

Let \(\phi(t)\) be the superconducting phase difference across the junction. The voltage across the junction is given by the Josephson relation:
\[
V(t) \;=\; \frac{\hbar}{2e}\,\frac{d\phi}{dt}.
\]
The circuit obeys Kirchhoff’s current law:
\[
I_{\text{bias}}(t)
\;=\;
I_c\,\sin(\phi(t))
\;+\;
\frac{V(t)}{R}
\;+\;
C \,\frac{dV(t)}{dt}
\;+\;
I_{\text{noise}}(t),
\]
where \(I_{\text{noise}}(t)\) represents thermal (Johnson–Nyquist) noise from the resistor. Substituting \(V(t) = \tfrac{\hbar}{2e}\dot{\phi}(t)\) yields a second‐order ODE in \(\phi\). We turn it into a **stochastic differential equation** (SDE) by modeling \(I_{\text{noise}}(t)\) as a random forcing, typically white Gaussian noise.

### 1.2 Dimensionless Form

We define a characteristic **plasma frequency** \(\omega_{p} = \sqrt{\tfrac{2eI_c}{\hbar C}}\) and time‐scale \(\tau = \omega_p\,t\). We also scale current by \(I_c\). The normalized system becomes:

\[
\begin{cases}
\frac{d\phi}{d\tau} = v,\\[5pt]
\frac{dv}{d\tau}
=
\gamma_{\mathrm{DC}} + \gamma_{\mathrm{AC}}\sin(\Omega\,\tau)
- \sin(\phi)
- \frac{1}{\beta_c}\,v
+ \eta(\tau),
\end{cases}
\]
where:

- \(\phi(\tau)\) is the dimensionless phase,
- \(v(\tau)\equiv \tfrac{d\phi}{d\tau}\) is the dimensionless voltage,
- \(\beta_c = \omega_p R C\) is the **Stewart–McCumber** parameter,
- \(\gamma_{\mathrm{DC}}, \gamma_{\mathrm{AC}}\) are dimensionless currents \(\tfrac{I_{\mathrm{DC}}}{I_c}\) and \(\tfrac{I_{\mathrm{AC}}}{I_c}\),
- \(\Omega = \tfrac{\omega_{\text{drive}}}{\omega_p}\) is the **dimensionless AC frequency**,
- \(\eta(\tau)\) is a Gaussian white noise of some amplitude reflecting the resistor’s thermal noise.

## 2. Numerical Integration: Heun Method

We solve the above SDE using a **Heun** (explicit trapezoid) scheme, which is a stochastic Runge–Kutta‐type integrator that attains **strong order** 1.0 for pathwise solutions. Each step from \(\tau_n\) to \(\tau_{n+1}\) proceeds:

1. **Predictor**: an explicit Euler step to guess the next state,
2. **Corrector**: re‐evaluate the drift at that guess and average with the old drift,
3. **Noise**: we add a random increment \(\sqrt{\Delta\tau}\,\mathcal{N}(0,1)\) in the equation for \(v\).

This approach is straightforward yet yields pathwise accuracy beyond Euler–Maruyama’s 0.5 order. In code, we do:

```python
# Pseudocode for one step of Heun
phi_n, v_n = phi[i], v[i]
dW = np.random.randn() * sqrt(dt)

# drift at old state
fphi_n = v_n
fv_n = gammaDC + gammaAC*sin(Omega*tau_n) - sin(phi_n) - (1/beta_c)*v_n

# predictor
phi_star = phi_n + fphi_n*dt
v_star   = v_n   + fv_n*dt + noiseAmp*dW

# drift at predicted state
fphi_star = v_star
fv_star   = gammaDC + gammaAC*sin(Omega*(tau_n+dt)) - sin(phi_star) - (1/beta_c)*v_star

# corrector
phi[i+1] = phi_n + 0.5*(fphi_n + fphi_star)*dt
v[i+1]   = v_n   + 0.5*(fv_n   + fv_star)*dt + noiseAmp*dW
```

In the actual code, we also include the **(1 - frac_4pi) * sin(phi_n) + frac_4pi * sin(phi_n/2)** for partial \(4\pi\)-periodic conduction.

## 3. Implementation Overview

1. **Dependencies**  
   - Python 3.7+  
   - NumPy for arrays, math functions  
   - Matplotlib for plotting  
   - (Optional) [tqdm](https://github.com/tqdm/tqdm) for the progress bar  
   - (Optional) [joblib](https://joblib.readthedocs.io/) for parallel runs  
   - (Optional) [Numba](https://numba.pydata.org/) for JIT compilation of the integrator

2. **Workflow**  
   - We parse the physical parameters \((R, C, T, I_c, f_{\text{drive}})\).  
   - Convert them to dimensionless forms \((\beta_c, \Omega, \dots)\).  
   - For each DC bias current \(\gamma_{\mathrm{DC}}\) in a user‐provided range, run multiple “ensemble” realizations of the SDE with random noise, discard an initial transient, then average the final dimensionless voltage.  
   - This yields an I–V curve that we can plot.

3. **Parallelization**  
   - If joblib is installed, we run each ensemble realization in parallel. Each run is independent because the random seeds differ.  
   - If not installed, we fall back to a serial loop.

4. **Numba Acceleration**  
   - We define a function `rsj_heun_run_numba` (decorated with `@njit`) to do the time stepping. The loop is compiled to native code.  
   - Gains can be large for bigger time spans or many ensemble runs.

## 4. Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/username/RCSJ_Simulation.git
   cd RCSJ_Simulation
   ```

2. (Optional) Create a virtual environment and install dependencies:
   ```bash
   python -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```

3. Run the main script:
   ```bash
   python rcsj_simulation.py
   ```
   This will print out some parameters (e.g., \(\beta_c\), \(\omega_p\), \(\Omega\)) and produce an I–V plot as dimensionless voltage vs. DC current.

## 5. Explanation of Key Code Segments

**`rsj_heun_run_numba(...)`**  
Implements the Heun method in dimensionless time steps. Inside the loop:
- `phi_n, v_n`: store current phase & voltage,
- `dW = np.random.randn() * sqrt(dt)`: single normal increment for the step,
- `fv_n`: drift function for voltage, capturing DC + AC + nonlinear sinusoidal conduction + resistor damping,
- predictor/corrector pattern for `(phi, v)`.

**`rcsj_iv_curve(...)`**  
- Loops over multiple DC bias points \(\gamma_{\mathrm{DC}}\).  
- At each DC bias, does `ensembleRuns` calls to `rsj_heun_run_numba`, storing final or steady average voltages.  
- Updates the noise amplitude each time, based on the updated average voltage (like in the original RCSJ code).

**Dimensionless Noise**  
- The noise amplitude depends on \(\mathrm{coth}\bigl(\tfrac{e v_0}{2k_B T}\bigr)\). This is the Johnson–Nyquist plus shot noise approximation for a given average dimensionless voltage \(v_0\).

## 6. Typical Results

When running `rcsj_simulation.py`, you should see:
- A textual progress bar (if tqdm is installed),
- For each DC bias, a short message or updated bar about the simulation,
- A final plot of \(\langle v\rangle / \Omega\) vs. \(\gamma_{\mathrm{DC}}\).  

Shapiro steps or other interesting structures in the dimensionless voltage may appear, depending on the AC amplitude \(\gamma_{\mathrm{AC}}\), frequency \(\Omega\), and temperature‐derived noise strength.

## 7. References

- K. K. Likharev, \emph{Dynamics of Josephson Junctions and Circuits}, Gordon and Breach, 1986.  
- R. F. Voss and R. A. Webb, \emph{Macroscopic Quantum Tunneling in Josephson Junctions}, Phys. Rev. Lett. 47, 265–268 (1981).  
- P. Hanggi and F. Marchesoni, \emph{Artificial Brownian motors: Controlling transport on the nanoscale}, Rev. Mod. Phys. 81, 387–442 (2009).  
- Rösler, Andreas. \emph{Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations}, SIAM J. Numer. Anal. 48, 3 (2010).  

---

**We hope you find this code and documentation helpful in studying the stochastic RCSJ model.** Feel free to open Issues or Pull Requests to improve or extend the code with additional integrators, higher‐order methods, or more advanced noise modeling.
