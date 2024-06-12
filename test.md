Given $d N^U_{obs}(z_s)$ is the number of detectable gravitational wave (GW) events from sources at red-shift $z_s$, then, rate of observing GW (number of such events happening per unit time) is given by,

$$
\begin{equation}
\begin{split}
\mathcal{R}_U &= \int_{z_{min}}^{z_{max}} \frac{d N^U_{obs}(z_s)}{d t} \\
&= \int_{z_{min}}^{z_{max}} \frac{d N^U_{obs}(z_s)}{d t \;dV_c} \frac{dV_c}{dz_s} dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d N^U_{obs}(z_s)}{\left(d t/(1+z_s)\right) \;dV_c}\frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d N^U_{obs}(z_s)}{d\tau \;dV_c}\frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s
\end{split} \tag{1} 
\end{equation}
$$

Note that if $d\tau=dt/(1+z_s)$ is considered the proper time and it can be converted from the time at detector frame $dt$ using the time-dilation factor $(1+z_s)$. Consequently, $\frac{d N^U_{obs}(z_s)}{d t \;dV_c}$ and $\frac{d N^U_{obs}(z_s)}{d \tau \;dV_c}$ are the observed merger rate density at detector-frame and source-frame respectively. We want to use the $R^U_{obs}(z_s)=\frac{d N^U_{obs}(z_s)}{d \tau \;dV_c}$ for our analysis as most observational papers and the output of theoretical predictions are in the source-frame. $4\pi \frac{dV_c}{dz_s}dz_s$ is considered a spherical-shell volume element in co-moving coordinates at red-shift $z_s$. So, the rate expression simplifies to integrating $\text{density} \times \text{time dilation effect} over the shell volume element.

Note: $\frac{dV_c}{dz_s}$ is the differential co-moving volume at red-shift $z_s$ and you can get the value by using `astropy` cosmology package for a given cosmology.

$$
\begin{equation}
\begin{split}
\mathcal{R}_U &= \int_{z_{min}}^{z_{max}} \frac{d N^U_{obs}(z_s)}{d \tau\;dV_c} \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
&= \int_{z_{min}}^{z_{max}} R^U_{obs}(z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s
\end{split} \tag{2}
\end{equation}
$$

Equation 2 is exactly same to the of equation A2 from [WIERDA et al. 2021](https://arxiv.org/pdf/2106.06303.pdf).

Now, we want to re-write it in terms intrinsic merger rate distribution $R(z_s)$. $R(z_s)$ is the merger rate density distribution regardless of whether the source is detectable or not.

$$
\begin{equation}
\begin{split}
\mathcal{R}_U &= \int_{z_{min}}^{z_{max}} R^U_{obs}(z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
&= \int_{z_{min}}^{z_{max}} R(z_s) P(obs|z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
\end{split} \tag{3.1}
\end{equation}
$$

$$
\begin{equation}
\begin{split}
\mathcal{R}_U &= \mathcal{N}^U \int_{z_{min}}^{z_{max}} P(z_s) P(obs|z_s) dz_s \\
\end{split} \tag{3.3}
\end{equation}
$$

$P(obs|z_s)$ is the probability of observing a GW event at red-shift $z_s$. Here, it is considered marginalized over all the gravitational wave parameters except the red-shift, i.e., $P(obs|z_s) = \int_{\theta} d\theta P(obs|z_s, \theta) P(\theta)$.

The merger rate density becomes the prior distribution $P(z_s)$ for the red-shift and it is normalized over the red-shift range. So the normalizing factor is given by,

$$
\begin{equation}
\begin{split}
\mathcal{N}^U = \int_{z_{min}}^{z_{max}} R(z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s 
\end{split} \tag{3.4}
\end{equation}
$$

Therefore,

$$
\begin{equation}
\begin{split}
P(z_s) = \frac{R(z_s)}{\mathcal{N}^U} \frac{1}{(1+z_s)} \frac{dV_c}{dz_s}
\end{split} \tag{3.5}
\end{equation}
$$

This is done because we want to sample the red-shift values from this prior distribution, which later will be used in the monte-carlo integration to estimate the rate.