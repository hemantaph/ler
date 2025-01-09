### Optical depth for EPL+Shear lens

* $\phi_{\text{EPL}'}=P(SL|z_s, z_l, \theta_L)$: area of the double caustic.
    * $z_s$: source redshift
    * $z_l$: lens redshift
* $\theta_L\in \{\sigma, q, e_1,e_2, \gamma_1, \gamma_2, \gamma\}$
    * $\sigma$: velocity dispersion
    * $q$: axis-ratio
    * $e_1, e_2$: ellipticity
    * $\gamma_1, \gamma_2$: shear
    * $\gamma$: mass-density spectral index

The probability of strong lensing for a SIE lens is given by,

$$
\begin{equation}
\begin{split}
P(SL|z_s) &=  \int \frac{\phi_{\text{EPL}'}}{4\pi} \frac{dN(z_l)}{dV_c d\theta_L} \frac{dV_c}{dz_l} dz_l d\theta_L \\
P(SL|z_s) &=  \int \frac{\phi_{\text{EPL}'}}{4\pi} \big< n \big>_{\theta_L} P(\theta_L|z_l) \frac{dV_c}{dz_l} dz_l d\theta_L \\
P(SL|z_s) &=  \int \frac{\phi_{\text{EPL}'}}{4\pi} n_o P(\theta_L|z_l) \frac{dV_c}{dz_l} dz_l d\theta_L \\
& \;\;\;\;\; \text{$\phi_{\text{EPL}'}$ is caustic double area.}
\end{split} \tag{12a}
\end{equation}
$$

$$
\begin{equation}
\begin{split}
P(SL|z_s) &=  \int_0^{z_s} \Phi_{\text{EPL}'}(z_l) dz_l \\
& \text{where, $\Phi_{\text{EPL}'}(z_l)= \int \frac{\phi_{\text{EPL}'}}{4\pi} n_o P(\theta_L|z_l) \frac{dV_c}{dz_l} d\theta_L$ }. \\
& \;\;  \Phi_{\text{EPL}'}(z_l)= \left< \frac{\phi_{\text{EPL}'}}{4\pi} n_o \frac{dV_c}{dz_l} \right>_{\theta_L\in P(\theta_L|z_l)} \\
& \text{If $\theta_L$ is independent of $z_l$, then} \\
& \;\;  \Phi_{\text{EPL}'}(z_l)= \left< \frac{\phi_{\text{EPL}'}}{4\pi} n_o \frac{dV_c}{dz_l} \right>_{\theta_L\in P(\theta_L)}
\end{split} 
\end{equation}
$$