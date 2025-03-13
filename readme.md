# RCSJ Simulation in Python

This repository provides a Python implementation of a **Resistively and Capacitively Shunted Josephson Junction (RCSJ)** model, including a thermal noise source from the shunt resistor, an RF source and an optional $4\pi$-periodic current-phase relation. The code integrates a **stochastic differential equation** (SDE) for the junction’s phase using a Heun‐type (explicit trapezoidal) method in dimensionless form. All  physical parameters for the JJ can be chosen. The scripts are subdivided in 4 parts: Functions, Parameters, Execution, Plotting. Just change the parameters section, if you just want to use it. If you want to compute shapiro maps or several IV-curves, write a for-loop around the execution line.

## Physical Model and Noise Derivation

Consider a Josephson junction with critical current $I_{c}$. It is placed in parallel with a resistor $R$ and a capacitor $C$, all driven by a bias current $I_{\mathrm{bias}}(t)$. The voltage across the junction is:

$$V(t) = \frac{\hbar}{2e} \frac{d\phi(t)}{dt},$$

where $\phi(t)$ is the superconducting phase difference. By Kirchhoff’s Current Law, the circuit equation becomes:

$$I_{\mathrm{bias}}(t) = I_{c} \sin\bigl(\phi(t)\bigr) + \frac{V(t)}{R} + C \frac{dV(t)}{dt} + I_{\mathrm{noise}}(t).$$

### Thermal Noise: Johnson vs. Shot

At temperature $T$, the resistor $R$ produces **Johnson–Nyquist** noise. Its low-frequency spectral density in the simplest linear regime is:

$$S_{I}^\text{(Johnson)} = 4 k_{B} T  \frac{1}{R}, \quad(\text{A}^2/\text{Hz}),$$

where $k_B$ is the Boltzmann constant. However, at higher voltages $\lvert V\rvert$, the conduction can also exhibit **shot‐like** noise, approximated by:

$$S_{I}^\text{(shot)}(V) = 2 e \frac{\lvert V\rvert}{R}.$$

A more **general** interpolation formula is:

$$S_{I}(V) \approx 2 e \frac{\lvert V\rvert}{R} \coth\!\Bigl(\frac{e \lvert V\rvert}{2 k_{B} T}\Bigr).$$

- If $\tfrac{e \lvert V\rvert}{2 k_{B} T} \ll 1$, $\coth(x)\approx\frac{1}{x}$, so $S_{I}(V)\approx 4 k_B T (1/R)$ => **Johnson** regime.
- If $\tfrac{e \lvert V\rvert}{2 k_{B} T} \gg 1$, $\coth(x)\approx 1$, so $S_{I}(V)\approx 2 e \lvert V\rvert/R$ => **shot** regime.

Thus, for small voltages we recover Johnson noise, and for large voltages we approach shot noise.

The random forcing $\eta(\tau)$ models the resistor’s thermal (or shot) noise in dimensionless form. Suppose its amplitude depends on the junction’s average dimensionless voltage $v_{0}$. Then:

$$\eta(\tau)\approx\text{noiseAmp}\times\mathcal{N}(0,1) \sqrt{d\tau},$$

where

$$S_{I}\approx2 e \frac{v_{0}}{R} \coth\!\Bigl(\frac{e v_{0}}{2 k_{B} T}\Bigr)\quad\Longrightarrow\quad\text{noiseAmp}=\sqrt{\frac{2 S_{I}}{I_{c}^{2} \omega_{p} \Delta\tau}}.$$

In each time step, we compute the final/average dimensionless voltage $v_{0}$, then use it to update `noiseAmp`. This ensures the noise transitions smoothly between Johnson or shot regimes depending on the operating point of the junction.

### Dimensionless RCSJ Equation

To eliminate physical constants, define:

- $\omega_{p}=\sqrt{\tfrac{2 e I_c}{\hbar C}}$: the **plasma frequency**,
- $\tau=\omega_{p} t$: a dimensionless time,
- $\gamma_{\mathrm{DC}}=\tfrac{I_{\mathrm{DC}}}{I_{c}},  \gamma_{\mathrm{AC}}=\tfrac{I_{\mathrm{AC}}}{I_{c}}$: normalized DC/AC drives,
- $\Omega=\tfrac{\omega_{\mathrm{drive}}}{\omega_p}$: dimensionless AC frequency,
- $\beta_{c} = \omega_{p} R C$: the **Stewart–McCumber** parameter.

Let $\phi(\tau)$ be the dimensionless phase, and define

$$v(\tau)=\frac{d\phi(\tau)}{d\tau}.$$

We can rewrite the circuit law as:

$$\frac{d\phi}{d\tau}=v,\quad\frac{dv}{d\tau}=\gamma_{\mathrm{DC}}+\gamma_{\mathrm{AC}}\sin(\Omega \tau) - I_{\mathrm{JJ}}(\phi) - \frac{v}{\beta_{c}} + \eta(\tau).$$

Here, the junction current $I_{\mathrm{JJ}}(\phi)$ can be partially or wholly $4\pi$-periodic:

$$I_{\mathrm{JJ}}(\phi) = \bigl(1-\text{frac4pi}\bigr) \sin(\phi) + \bigl(\text{frac4pi}\bigr) \sin\!\bigl(\tfrac{\phi}{2}\bigr).$$

- If `frac4pi=0`, we recover the usual $2\pi$-periodic $\sin(\phi)$.
- If `frac4pi=1`, we get $\sin(\phi/2)$, doubling the period to $4\pi$.

## Voltage Transient

When the system starts from some initial condition (e.g. $\phi=0$, $v=0$), there is often a **transient** phase where the junction’s voltage quickly evolves from this initial state to a more statistically stationary regime. Even in the "steady state", the Voltage is always oscillating due to the 2nd Josephson equation. In practice, we:

- **Discard** some initial fraction of the simulation steps (e.g. 2/3) to remove transient effects. You should play around with the time step, discard ratio and the integration time to find a good tradeoff. The integration time is coupled to the frequency of the RF drive, because this is the frequency at which the Voltage will oscillate at the shapiro steps. If you decrease the RF frequency, the integration time will be automatically chosen longer. If you do not use RF at all, put both the frequency and amplitude to 0.
- **Average** the final portion of the dimensionless voltage $v(\tau)$ to obtain a representative “steady” or “long-time” average. This average is typically what we refer to as the measured or effective DC voltage. That is how we get an $I\!-\!V$ (or $\gamma\!-\!v$) characteristic.

## Heun Integrator in Python

We simulate over dimensionless timesteps $\Delta\tau$. In each step:

Known $\phi_{n},v_{n}$. Generate standard normal $\xi$. Let $dW=\xi \sqrt{\Delta\tau}$.
Drift at old state:

   $$f_{\phi}=v_{n},\quad f_{v}=\gamma_{\mathrm{DC}}+\gamma_{\mathrm{AC}}\sin(\Omega \tau_{n})-I_{\mathrm{JJ}}(\phi_{n})-\tfrac{v_{n}}{\beta_{c}}.$$

**Predictor** (Euler step):
   
   $$\phi_{\star}=\phi_{n}+f_{\phi} \Delta\tau,\quad v_{\star}=v_{n}+f_{v} \Delta\tau+(\text{noiseAmp}) dW.$$

Evaluate drift at predicted state:

   $$f_{\phi,\star}=v_{\star},\quad f_{v,\star}=\gamma_{\mathrm{DC}}+\gamma_{\mathrm{AC}}\sin\bigl(\Omega(\tau_{n}+\Delta\tau)\bigr)-I_{\mathrm{JJ}}(\phi_{\star})-\tfrac{v_{\star}}{\beta_{c}}.$$

 **Corrector**:

   $$\phi_{n+1}=\phi_{n}+\tfrac12\bigl(f_{\phi}+f_{\phi,\star}\bigr)\Delta\tau,\quad v_{n+1}=v_{n}+\tfrac12\bigl(f_{v}+f_{v,\star}\bigr)\Delta\tau+(\text{noiseAmp}) dW.$$


