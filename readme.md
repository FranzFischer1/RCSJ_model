# RCSJ Simulation in Python

This repository provides a Python implementation of a **Resistively and Capacitively Shunted Josephson Junction (RCSJ)** model, including a thermal noise source from the shunt resistor and an optional $4\pi$-periodic current-phase relation. The code integrates a **stochastic differential equation** (SDE) for the junction’s phase, using a Heun‐type (explicit trapezoidal) method in dimensionless form.

## 1  Physical Model and Noise Derivation

Consider a Josephson junction with critical current $I_{c}$. It is placed in parallel with a resistor $R$ and a capacitor $C$, all driven by a current $I_{\mathrm{bias}}(t)$. The voltage across the junction is  
$$V (t) = \frac{\hbar}{2 e}\frac{d \phi (t)}{d t}$$

where $\phi (t)$ is the superconducting phase difference. By Kirchhoff’s Current Law, the circuit equation becomes  
$$I_{\mathrm{bias}} (t)= I_{c} \sin \bigl(\phi (t)\bigr)+ \frac{V (t)}{R}+ C \frac{d V (t)}{d t}+ I_{\mathrm{noise}} (t)$$

### 1.1  Thermal Noise Formula

The resistor $R$ at temperature $T$ generates **Johnson–Nyquist noise** with a spectral density  

$$S_{I} = 4 k_{B} T \frac{1}{R} \quad (\text{A}^2/\text{Hz})$$

where $k_{B}$ is Boltzmann’s constant. If additional shot noise or more advanced corrections apply, one may write a more general formula such as  

$$S_{I} (V ) \approx 2 e \frac{\lvert V \rvert}{R} \, \coth \!\Bigl(\frac{e\lvert V \rvert}{2 k_{B} T}\Bigr)$$

which **interpolates** between a linear Johnson regime at small $V$ and shot‐like behavior at large $V$. In dimensionless form, we factor out $I_{c}$ and the typical voltage scale.

### 1.2  Dimensionless RCSJ Equation

To eliminate physical constants, define:

- $\omega_{p} = \sqrt{\dfrac{2 e I_{c}}{\hbar C}}$ : the **plasma frequency**,
- $\tau = \omega_{p} t$ : a dimensionless time,
- $\gamma_{\mathrm{DC}} = \dfrac{I_{\mathrm{DC}}}{I_{c}}$ , $\gamma_{\mathrm{AC}} = \dfrac{I_{\mathrm{AC}}}{I_{c}}$ : normalized DC/AC drives,
- $\Omega = \dfrac{\omega_{\mathrm{drive}}}{\omega_{p}}$ : dimensionless AC frequency,
- $\beta_{c} = \omega_{p} \, R \, C$ : the Stewart–McCumber parameter.

Then, letting $\phi (\tau)$ be the dimensionless phase, define  

$$v (\tau) = \frac{d \phi (\tau)}{d \tau}$$

so that the circuit law becomes a system:

$$\frac{d \phi}{d \tau} = v \qquad
\frac{d v}{d \tau}
= \gamma_{\mathrm{DC}} + \gamma_{\mathrm{AC}} \,\sin \!\bigl(\Omega \,\tau\bigr)
- I_{\mathrm{JJ}}(\phi) - \frac{1}{\beta_{c}} \, v + \eta (\tau)$$

Here, the junction current $I_{\mathrm{JJ}}(\phi)$ can be purely **$2 \pi$**‐periodic (the usual \(\sin (\phi)\)) or partially $4 \pi$–periodic:

$$I_{\mathrm{JJ}}(\phi) = \bigl(1 - \text{frac\_4pi}\bigr)\,\sin\!\bigl(\phi\bigr)
\;+\;
\bigl(\text{frac\_4pi}\bigr)\,\sin\!\bigl(\tfrac{\phi}{2}\bigr).$$

If `frac_4pi = 0`, we recover the standard \(\sin(\phi)\). If `frac_4pi = 1`, we get \(\sin(\phi/2)\), effectively doubling the period to $4 \pi$.

### 1.3  Adding the Noise

The random forcing \(\eta (\tau)\) stems from the resistor’s thermal noise in dimensionless form. Suppose the amplitude depends on the average voltage, say

$$\eta (\tau) \approx \text{noiseAmp} \times \mathcal{N}(0,1)\,\sqrt{d\tau}$$

We can update `noiseAmp` each time if the system’s average dimensionless voltage $v_{0}$ changes significantly, e.g. from the interpolation formula
$$S_{I} \approx 2 e \,\frac{v_{0}}{R}\,\coth \!\Bigl(\frac{e \,v_{0}}{2 k_{B} T}\Bigr)\quad
\text{noiseAmp} = \sqrt{\frac{2 S_{I}}{I_{c}^{2} \,\omega_{p}\,\Delta \tau}}$$

## 2  Heun Integrator in Python

We solve the system in dimensionless time steps \(\Delta \tau\). In each step:

1. \(\phi_{n}, v_{n}\) are known. Generate a standard normal \(\xi\) and let $dW = \xi \sqrt{\Delta \tau}$.
2. Drift at old state:
   $$f_{\phi} = v_{n}, \quad
   f_{v} = \gamma_{\mathrm{DC}} + \gamma_{\mathrm{AC}}\sin (\Omega\, \tau_{n})
           - I_{\mathrm{JJ}}(\phi_{n})
           - \tfrac{1}{\beta_{c}} v_{n}.$$
   
4. **Predictor** (Euler):
   $$\phi_{\star} = \phi_{n} + f_{\phi}\,\Delta \tau,
   \quad
   v_{\star}   = v_{n} + f_{v}\,\Delta \tau + (\text{noiseAmp}) \, dW.
   $$
5. Evaluate drift at predicted state:
   $$
   f_{\phi,\star} = v_{\star},
   \quad
   f_{v,\star} = \gamma_{\mathrm{DC}} + \gamma_{\mathrm{AC}}\sin (\Omega(\tau_{n}+\Delta \tau))
                 - I_{\mathrm{JJ}}(\phi_{\star})
                 - \tfrac{1}{\beta_{c}} v_{\star}.
   $$
6. **Corrector**:
   $$
   \phi_{n+1} = \phi_{n} + \tfrac{1}{2}\bigl(f_{\phi} + f_{\phi,\star}\bigr)\,\Delta \tau,
   \quad
   v_{n+1}   = v_{n} + \tfrac{1}{2}\bigl(f_{v} + f_{v,\star}\bigr)\,\Delta \tau
               + (\text{noiseAmp}) \, dW.
   $$
