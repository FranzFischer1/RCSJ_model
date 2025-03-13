Below is an updated, more detailed README text suitable for a GitHub project. It:

1. Uses the standard `"$"` and `"$$"` LaTeX math delimiters with **no extra spaces** after the dollar sign.
2. Shows how the resistor noise formula converges to shot or Johnson noise in the corresponding limits.
3. Provides an explanation for the transient voltage and why we compute the mean steady voltage.
4. Includes the changes requested.

Feel free to adapt filenames, function names, or references to match your repository structure.

---

# RCSJ Simulation in Python

This repository provides a Python implementation of a **Resistively and Capacitively Shunted Josephson Junction (RCSJ)** model, including a thermal noise source from the shunt resistor and an optional \(4\pi\)-periodic current-phase relation. The code integrates a **stochastic differential equation** (SDE) for the junction’s phase using a Heun‐type (explicit trapezoidal) method in dimensionless form.

## 1  Physical Model and Noise Derivation

Consider a Josephson junction with critical current $I_{c}$. It is placed in parallel with a resistor $R$ and a capacitor $C$, all driven by a bias current $I_{\mathrm{bias}}(t)$. The voltage across the junction is:

$$V(t) = \frac{\hbar}{2e} \frac{d\phi(t)}{dt},$$

where \(\phi(t)\) is the superconducting phase difference. By Kirchhoff’s Current Law, the circuit equation becomes:

$$I_{\mathrm{bias}}(t)\;=\;I_{c}\,\sin\bigl(\phi(t)\bigr)\;+\;\frac{V(t)}{R}\;+\;C\,\frac{dV(t)}{dt}\;+\;I_{\mathrm{noise}}(t).$$

### 1.1  Thermal Noise: Johnson vs. Shot

At temperature \(T\), the resistor \(R\) produces **Johnson–Nyquist** noise. Its low-frequency spectral density in the simplest linear regime is:

$$S_{I}^\text{(Johnson)} = 4\,k_{B}\,T \,\frac{1}{R}, \quad(\text{A}^2/\text{Hz}),$$

where \(k_B\) is the Boltzmann constant. However, at higher voltages \(\lvert V\rvert\), the conduction can also exhibit **shot‐like** noise, approximated by:

$$S_{I}^\text{(shot)}(V) = 2\,e\,\frac{\lvert V\rvert}{R}.$$

A more **general** interpolation formula is:

$$S_{I}(V)\;\approx\;2\,e\,\frac{\lvert V\rvert}{R}\,\coth\!\Bigl(\frac{e\,\lvert V\rvert}{2\,k_{B}\,T}\Bigr).$$

- If \(\tfrac{e\,\lvert V\rvert}{2\,k_{B}\,T} \ll 1\), \(\coth(x)\approx\frac{1}{x}\), so \(S_{I}(V)\approx 4\,k_B\,T\,(1/R)\) => **Johnson** regime.
- If \(\tfrac{e\,\lvert V\rvert}{2\,k_{B}\,T} \gg 1\), \(\coth(x)\approx 1\), so \(S_{I}(V)\approx 2\,e\,\lvert V\rvert/R\) => **shot** regime.

Thus, for small voltages we recover Johnson noise, and for large voltages we approach shot noise.

### 1.2  Dimensionless RCSJ Equation

To eliminate physical constants, define:

- \(\omega_{p}=\sqrt{\tfrac{2\,e\,I_c}{\hbar\,C}}\): the **plasma frequency**,
- \(\tau=\omega_{p}\,t\): a dimensionless time,
- \(\gamma_{\mathrm{DC}}=\tfrac{I_{\mathrm{DC}}}{I_{c}},\;\;\gamma_{\mathrm{AC}}=\tfrac{I_{\mathrm{AC}}}{I_{c}}\): normalized DC/AC drives,
- \(\Omega=\tfrac{\omega_{\mathrm{drive}}}{\omega_p}\): dimensionless AC frequency,
- \(\beta_{c} = \omega_{p}\,R\,C\): the **Stewart–McCumber** parameter.

Let \(\phi(\tau)\) be the dimensionless phase, and define

$$v(\tau)=\frac{d\phi(\tau)}{d\tau}.$$

We can rewrite the circuit law as:

$$\frac{d\phi}{d\tau}=v,\quad\frac{dv}{d\tau}=\gamma_{\mathrm{DC}}+\gamma_{\mathrm{AC}}\sin(\Omega\,\tau)\;-\;I_{\mathrm{JJ}}(\phi)\;-\;\frac{v}{\beta_{c}}\;+\;\eta(\tau).$$

Here, the junction current \(I_{\mathrm{JJ}}(\phi)\) can be partially or wholly \(4\pi\)-periodic:

$$I_{\mathrm{JJ}}(\phi)\;=\;\bigl(1-\text{frac4pi}\bigr)\,\sin(\phi)\;+\;\bigl(\text{frac4pi}\bigr)\,\sin\!\bigl(\tfrac{\phi}{2}\bigr).$$

- If `frac4pi=0`, we recover the usual \(2\pi\)-periodic \(\sin(\phi)\).
- If `frac4pi=1`, we get \(\sin(\phi/2)\), doubling the period to \(4\pi\).

### 1.3  Adding the Noise and Updating its Amplitude

The random forcing \(\eta(\tau)\) models the resistor’s thermal (or shot) noise in dimensionless form. Suppose its amplitude depends on the junction’s average dimensionless voltage \(v_{0}\). Then:

$$\eta(\tau)\approx\text{noiseAmp}\times\mathcal{N}(0,1)\,\sqrt{d\tau},$$

where

$$S_{I}\approx2\,e\,\frac{v_{0}}{R}\,\coth\!\Bigl(\frac{e\,v_{0}}{2\,k_{B}\,T}\Bigr)\quad\Longrightarrow\quad\text{noiseAmp}=\sqrt{\frac{2\,S_{I}}{I_{c}^{2}\,\omega_{p}\,\Delta\tau}}.$$

In each time step, we compute the final/average dimensionless voltage \(v_{0}\), then use it to update `noiseAmp`. This ensures the noise transitions smoothly between Johnson or shot regimes depending on the operating point of the junction.

## 2  Heun Integrator in Python

We simulate over dimensionless timesteps \(\Delta\tau\). In each step:

1. Known \(\phi_{n},v_{n}\). Generate standard normal \(\xi\). Let \(dW=\xi\,\sqrt{\Delta\tau}\).
2. Drift at old state:

   $$f_{\phi}=v_{n},\quad f_{v}=\gamma_{\mathrm{DC}}+\gamma_{\mathrm{AC}}\sin(\Omega\,\tau_{n})-I_{\mathrm{JJ}}(\phi_{n})-\tfrac{v_{n}}{\beta_{c}}.$$

3. **Predictor** (Euler step):
   
   $$\phi_{\star}=\phi_{n}+f_{\phi}\,\Delta\tau,\quad v_{\star}=v_{n}+f_{v}\,\Delta\tau+(\text{noiseAmp})\,dW.$$

4. Evaluate drift at predicted state:

   $$f_{\phi,\star}=v_{\star},\quad f_{v,\star}=\gamma_{\mathrm{DC}}+\gamma_{\mathrm{AC}}\sin\bigl(\Omega(\tau_{n}+\Delta\tau)\bigr)-I_{\mathrm{JJ}}(\phi_{\star})-\tfrac{v_{\star}}{\beta_{c}}.$$

5. **Corrector**:

   $$\phi_{n+1}=\phi_{n}+\tfrac12\bigl(f_{\phi}+f_{\phi,\star}\bigr)\Delta\tau,\quad v_{n+1}=v_{n}+\tfrac12\bigl(f_{v}+f_{v,\star}\bigr)\Delta\tau+(\text{noiseAmp})\,dW.$$

## 3  Voltage Transient and Why We Take the Mean

When the system starts from some initial condition (e.g. \(\phi=0\), \(v=0\)), there is often a **transient** phase where the junction’s voltage quickly evolves from this initial state to a more statistically stationary regime. In practice, we:

- **Discard** some initial fraction of the simulation steps (e.g. 2/3) to remove transient effects.
- **Average** the final portion of the dimensionless voltage \(v(\tau)\) to obtain a representative “steady” or “long-time” average. This average is typically what we refer to as the measured or effective DC voltage. That is how we get an $I\!-\!V$ (or \(\gamma\!-\!v\)) characteristic.

## 4  Example Code Changes

Below is an updated snippet demonstrating how to incorporate the partial $4\pi$ conduction and the dynamic noise amplitude:

```python
def rsj_heun_run(
    gammaDC, gammaAC, Omega, beta_c,
    phi0, v0, tSim, dt,
    noiseAmp, frac_4pi
):
    import numpy as np
    N = int(np.floor(tSim / (Omega*dt))) if Omega != 0 else int(np.floor(tSim / dt))
    tau_array = np.arange(N+1)*dt
    phi_array = np.zeros(N+1)
    v_array   = np.zeros(N+1)

    phi_array[0] = phi0
    v_array[0]   = v0

    for i in range(N):
        phi_n = phi_array[i]
        v_n   = v_array[i]
        tau_n = tau_array[i]

        dW = np.random.randn() * np.sqrt(dt)

        # drift at old state
        # partial 4pi conduction:
        Ijj = (1.0 - frac_4pi)*np.sin(phi_n) + frac_4pi*np.sin(phi_n/2.0)
        fphi_n = v_n
        fv_n   = gammaDC + gammaAC*np.sin(Omega*tau_n) - Ijj - (1.0/beta_c)*v_n

        # predictor (Euler)
        phi_star = phi_n + fphi_n*dt
        v_star   = v_n   + fv_n*dt + noiseAmp*dW

        # drift at predicted state
        Ijj_star = (1.0 - frac_4pi)*np.sin(phi_star) + frac_4pi*np.sin(phi_star/2.0)
        fphi_star = v_star
        fv_star   = gammaDC + gammaAC*np.sin(Omega*(tau_n+dt)) - Ijj_star - (1.0/beta_c)*v_star

        # corrector
        phi_array[i+1] = phi_n + 0.5*(fphi_n + fphi_star)*dt
        v_array[i+1]   = v_n   + 0.5*(fv_n   + fv_star)*dt + noiseAmp*dW

    return tau_array, phi_array, v_array
```

---

We hope this detailed derivation clarifies **why** the noise formula transitions from Johnson to shot noise, **how** the dimensionless RCSJ SDE is formed (including the partial $4\pi$ conduction), and **why** we discard transients and take a mean to get the final DC voltage.
