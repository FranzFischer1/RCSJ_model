# RCSJ Simulation in Python

This repository provides a Python implementation of a **Resistively and Capacitively Shunted Josephson Junction (RCSJ)** model. It accounts for thermal/shot noise from the shunt resistor, optionally includes a partially \(4\pi\)-periodic current-phase relation, and numerically integrates the resulting **stochastic differential equation (SDE)** for the junction’s phase using a Heun‐type (explicit trapezoidal) method in dimensionless form.

## 1  Physical Model and Noise Derivation

Consider a Josephson junction with critical current \$I_c\$. The junction is placed in parallel with a resistor \$R\$ and a capacitor \$C\$, all driven by a bias current \$I_{\mathrm{bias}}(t)\$. The voltage across the junction is 
\[
V(t) = \frac{\hbar}{2e}\,\frac{d\phi(t)}{dt},
\]
where \$\phi(t)\$ is the superconducting phase difference. By Kirchhoff’s Current Law, we have
\[
I_{\mathrm{bias}}(t)
= I_c \,\sin\bigl(\phi(t)\bigr)\;+\;\frac{V(t)}{R}\;+\;C\,\frac{dV(t)}{dt}\;+\;I_{\mathrm{noise}}(t).
\]

### 1.1  Noise Formula and Limits

The resistor \$R\$ at temperature \$T\$ typically generates **Johnson–Nyquist noise** with a low-frequency spectral density 
\[
S_I \;=\; 4\,k_B\,T\,\frac{1}{R}\quad(\text{A}^2/\text{Hz}),
\]
where \$k_B\$ is Boltzmann’s constant. However, when the junction’s voltage becomes large, **shot noise** may dominate. A common interpolation formula is:
\[
S_I(V)\;\approx\;2\,e\,\frac{|V|}{R}\;\coth\!\Bigl(\frac{e\,|V|}{2\,k_B\,T}\Bigr).
\]
- **Johnson‐like limit**: for small \$|V|\$, \$e|V|\ll k_B T\$, we have \$\coth(x)\approx 1/x\$ so
\[
S_I(V)\;\to\;4\,k_B\,T\,\frac{1}{R}\quad\text{(Johnson noise)}.
\]
- **Shot‐like limit**: for large \$|V|\$, \$e|V|\gg k_B T\$, \$\coth(x)\approx 1\$, so
\[
S_I(V)\;\to\;2\,e\,\frac{|V|}{R}\quad\text{(shot noise)}.
\]

In our simulation, we often define a dimensionless amplitude \(\text{noiseAmp}\) from \$S_I(V)\$ after scaling out \$I_c\$ and the plasma frequency, so the random forcing \$I_{\mathrm{noise}}(t)\$ is effectively 
\[
\text{noiseAmp}\times\mathcal{N}(0,1)\times\sqrt{\Delta \tau}\quad\text{in dimensionless units.}
\]

### 1.2  Dimensionless RCSJ Equation

We eliminate physical constants by defining:
- \$\omega_p = \sqrt{\tfrac{2\,e\,I_c}{\hbar\,C}}\$: the **plasma frequency**,
- \$\tau = \omega_p\,t\$: a dimensionless time,
- \$\gamma_{\mathrm{DC}} = \tfrac{I_{\mathrm{DC}}}{I_c}\$, \$\gamma_{\mathrm{AC}} = \tfrac{I_{\mathrm{AC}}}{I_c}\$: normalized DC/AC currents,
- \$\Omega = \tfrac{\omega_{\mathrm{drive}}}{\omega_p}\$: dimensionless AC drive frequency,
- \$\beta_c = \omega_p\,R\,C\$: the **Stewart–McCumber** parameter.

Let \$\phi(\tau)\$ be the dimensionless phase, and define
\[
v(\tau)\;=\;\frac{d\phi(\tau)}{d\tau}.
\]
Then the circuit equation becomes:

\[
\frac{d\phi}{d\tau}\;=\;v,\quad
\frac{dv}{d\tau}
\;=\;\gamma_{\mathrm{DC}}
+\gamma_{\mathrm{AC}}\sin\bigl(\Omega\,\tau\bigr)
-\,I_{\mathrm{JJ}}(\phi)
-\,\frac{v}{\beta_c}
+\,\eta(\tau),
\]
where \$\eta(\tau)\$ is the dimensionless random forcing from the resistor noise. If the junction has a mixture of \$2\pi\$ and \$4\pi\$ periodic current-phase relations, we can define
\[
I_{\mathrm{JJ}}(\phi)
\;=\;(1-\text{frac4pi})\,\sin(\phi)\;+\;\text{frac4pi}\,\sin\!\Bigl(\tfrac{\phi}{2}\Bigr).
\]
Setting \`frac4pi = 0\` recovers the standard \$2\pi\$ conduction (\$\sin(\phi)\$), while \`frac4pi = 1\` leads to \$\sin(\phi/2)\$, effectively doubling the period to \$4\pi\$.

### 1.3  Voltage Transient and Ensemble Mean

When we switch on a particular DC bias or noise amplitude, the dimensionless voltage \$v(\tau)\$ often undergoes a transient phase in which it moves from some initial condition to a quasi‐stationary state or a phase‐locked region. **We discard** an initial fraction of the simulation — e.g. the first \$30\%\$ or \$50\%\$ of time steps — to remove that transient. Then we average the final portion of \$v(\tau)\$ over time to get a representative “steady‐state” dimensionless voltage, or we run many “ensemble” simulations with different random seeds and compute the final time average. This yields the \emph{mean} dimensionless voltage for each DC bias setting, enabling an I–V (or \(\gamma_{\mathrm{DC}}\)–\(\langle v\rangle\)) curve that reflects the typical junction behavior.

## 2  Heun Integrator in Python

We solve the dimensionless SDE with time step \$\Delta\tau\$. For each step from \(\tau_n\) to \(\tau_{n+1}\):

1. We have \(\phi_n, v_n\). Generate a standard normal \(\xi\) and let \$dW = \xi\,\sqrt{\Delta\tau}\$.
2. **Drift at old state**:
   \[
   f_\phi = v_n,\quad
   f_v = \gamma_{\mathrm{DC}} + \gamma_{\mathrm{AC}}\,\sin(\Omega\,\tau_n)
          - I_{\mathrm{JJ}}(\phi_n)
          - \frac{v_n}{\beta_c}.
   \]
3. **Predictor** (Euler):
   \[
   \phi_\star = \phi_n + f_\phi\,\Delta\tau,\quad
   v_\star   = v_n + f_v\,\Delta\tau + (\text{noiseAmp})\,dW.
   \]
4. Evaluate drift at predicted state:
   \[
   f_{\phi,\star} = v_\star,\quad
   f_{v,\star} = \gamma_{\mathrm{DC}} + \gamma_{\mathrm{AC}}\sin(\Omega(\tau_n+\Delta\tau))
                 - I_{\mathrm{JJ}}(\phi_\star)
                 - \frac{v_\star}{\beta_c}.
   \]
5. **Corrector**:
   \[
   \phi_{n+1} = \phi_n + \frac{1}{2}\,(f_\phi + f_{\phi,\star})\,\Delta\tau,\quad
   v_{n+1}    = v_n + \frac{1}{2}\,(f_v + f_{v,\star})\,\Delta\tau + (\text{noiseAmp})\,dW.
   \]

This procedure is repeated for \$N\$ steps until we reach final time \(\tau = N\,\Delta\tau\). We then discard the transient portion and average or store the final dimensionless voltage.
