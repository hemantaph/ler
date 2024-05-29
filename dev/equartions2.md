# Posterior overlap

$$
\begin{equation} 
\begin{split}
\mathcal{O}_U^L &= \frac{P(H_L|\{d_1,d_2\})}{P(H_U|\{d_1,d_2\})} \\ \nonumber
&= \frac{P(\{d_1,d_2\}|H_L)P(H_L)}{P(\{d_1,d_2\}|H_U)P(H_U)} \\ \nonumber
&= \mathcal{B}_U^L \mathcal{P}_U^L \\ \nonumber
\end{split}
\end{equation}
$$

$$
\begin{equation} 
\begin{split}
\mathcal{B}_U^L &= \frac{P(\{d_1,d_2\}|H_L)}{P(\{d_1,d_2\}|H_U)} \\ \nonumber
&= \frac{\int d\theta P(\{d_1,d_2\}|\theta,H_L)P(\theta|H_L)}{\int d\theta P(d_1|\theta,H_U)P(\theta|H_U) \int d\theta P(d_2|\theta,H_U)P(\theta|H_U)} \\
&= \frac{ \int d\theta \frac{P(\theta|d_1) P(d_1)}{P_{pe}(\theta)} \frac{P(\theta|d_2) P(d_2)}{P_{pe}(\theta)} P_{astro}(\theta|H_L) }
{\int d\theta \frac{P(\theta|d_1) P(d_1)}{P_{pe}(\theta)} P_{astro}(\theta|H_U) \int d\theta \frac{P(\theta|d_2) P(d_2)}{P_{pe}(\theta)} P_{astro}(\theta|H_U)} \\
&= \frac{ \int d\theta \frac{P(\theta|d_1)}{P_{pe}(\theta)} \frac{P(\theta|d_2)}{P_{pe}(\theta)} P_{astro}(\theta|H_L) }
{\int d\theta \frac{P(\theta|d_1)}{P_{pe}(\theta)} \frac{P(\theta|d_2)}{P_{pe}(\theta)} P_{astro}(\theta|H_U)} \\
\end{split}
\end{equation}
$$

* $\theta \in \{m_{1,1}, m_{2,1}, d_{L,1}, \iota_1, t_1, \alpha_1, \delta_1, m_{1,2}, m_{2,2}, d_{L,2}, \iota_2, t_2, \alpha_2, \delta_2\}$

* $\theta \in \{m_{1,1}: \text{mass of the heavier BH; first event}, m_{2,1}: \text{mass of the lighter BH; first event}, d_{L,1}: \text{luminosity distance; first event}, \iota_1: \text{inclination angle; first event}, t_1: \text{coalescence time; first event}, \alpha_1: \text{right ascension; first event}, \delta_1: \text{declination; first event}, m_{1,2}: \text{mass of the heavier BH; second event}, m_{2,2}: \text{mass of the lighter BH; second event}, d_{L,2}: \text{luminosity distance; second event}, \iota_2: \text{inclination angle; second event}, t_2: \text{coalescence time; second event}, \alpha_2: \text{right ascension; second event}, \delta_2: \text{declination; second event}\}$

Reparameterization:
  
* $\theta \in \{m_{1,1}, m_{2,1}, d_{L,1}, \iota_1, t_1, \alpha_1, \delta_1, m_{1,2}, m_{2,2}, d_{L,2}, \iota_2, t_2, \alpha_2, \delta_2\}$


* $\theta \in \{m_{1,1}, m_{2,1}, d_{L,1}, \iota_1, t_1, \alpha_1, \delta_1, m_{1,1}+\Delta m_{1,2}, m_{2,1}+\Delta m_{2,2}, d_{L,2}, \iota_2, t_1+\Delta t_2, \alpha_1+\Delta \alpha_2, \delta_1+\Delta \delta_2\}$

* $\theta \in \{m_{1,1}, m_{2,1}, d_{L,1}, \iota_1, t_1, \alpha_1, \delta_1, \Delta m_{1,2}, \Delta m_{2,2}, d_{L,2}, \iota_2, \Delta t_2, \Delta\alpha_2, \Delta\delta_2\}$

Let's tackle each part of the integrand separately.

$P_{astro}(\theta|H_L)$: astrophysical prior on the parameters with lensed hypothesis considered. This is for the detectable events. So, the selection effect is already taken into account. The distribution preserves the correlation between the parameters. But it doesn't take into account the correlation between the coalescence time and sky location. 

$$
\begin{equation} 
\begin{split}
P_{astro}(\theta|H_L) &= P_{astro}(m_{1,1}, m_{2,1}, d_{L,1}, \iota_1, t_1, \alpha_1, \delta_1, m_{1,2}, m_{2,2}, d_{L,2}, \iota_2, t_2, \alpha_2, \delta_2|H_L) \\ \nonumber
\\
&= P_{astro}(m_{1,2}, m_{2,2}, \alpha_2, \delta_2, \iota_2| m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, d_{L,2}, t_1, t_2, H_L) \\ \nonumber
&\;\;\;\;\;\;\;\;\;P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, d_{L,2}, t_1, t_2|H_L) \\ \nonumber
\\
&= P_{astro}(\Delta m_{1,2}, \Delta m_{2,2}, \Delta \alpha_2, \Delta \delta_2, \Delta \iota_{2}| m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, d_{L,2}, t_1, t_1+\Delta t_2, H_L) \\ \nonumber
&\;\;\;\;\;\;\;\;\;P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, d_{L,2}, t_1, t_1+\Delta t_2|H_L) \\ \nonumber
\\
&= P_{astro}(0,0,0,0,0| m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, d_{L,2}, t_1,\Delta t_2, H_L) \\ \nonumber
&\;\;\;\;\;\;\;\;\;P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, d_{L,2}, t_1, \Delta t_2|H_L) \\ \nonumber
\\
&= 1\times P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, d_{L,2}, t_1, \Delta t_2|H_L) \\ \nonumber
\\
&= P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, d_{L,2}, \Delta t_2|H_L) P(t_1|H_L) \\ \nonumber
\\
&= P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, d_{L,2}, \Delta t_2|H_L) P(t_1) \\ \nonumber
\\
&= P_{astro}(m_{1}, m_{2}, \alpha, \delta, \iota, d_{L,1}, d_{L,2}, \Delta t_2|H_L) P(t_1) \\ \nonumber
\end{split}
\end{equation}
$$

In the last step, I have considered $t_1$ is independent of the other parameters and lensing hypothesis. This prior is a uniform distribution over the observation time window. Since for the lensing events the parameters $[m_{1}, m_{2}, \alpha, \delta, \iota]$ are same for the images, I dropped the event index at the last step of the derivation.

$P(\theta|d_1)$: Posterior distribution of the parameters given the data from the first GW event.

$$
\begin{equation} 
\begin{split}
P(\theta|d_1) &= P(m_{1,1}, m_{2,1},\alpha_1, \delta_1, \iota_1, d_{L,1}, t_1|d_1) \\ \nonumber
&= P(m_{1,1}, m_{2,1},\alpha_1, \delta_1, \iota_1, d_{L,1}|t_1, d_1) P(t_1|d_1) \\ \nonumber
\end{split}
\end{equation}
$$

$P(t_1|d_1)$: Posterior distribution of the coalescence time of the first GW event given the data from the first GW event. This is a delta function and peaks at the coalescence time of the first GW event, say $T_1$.

Therefore,

$$
\begin{equation} 
\begin{split}
P(\theta|d_1) &= P(m_{1,1}, m_{2,1},\alpha_1, \delta_1, \iota_1, d_{L,1}|T_1, d_1) P(T_1) \\ \nonumber
\end{split}
\end{equation}
$$

Similarly, for the second event,

$$
\begin{equation} 
\begin{split}
P(\theta|d_2) &= P(m_{1,2}, m_{2,2},\alpha_2, \delta_2, \iota_2, d_{L,2}|T_2, d_2) P(T_2) \\ \nonumber
\end{split}
\end{equation}
$$

And, $P_{astro}(\theta|H_L)$ now becomes,

$$
\begin{equation} 
\begin{split}
P_{astro}(\theta|H_L) &= P_{astro}(m_{1}, m_{2}, \alpha, \delta, \iota, d_{L,1}, d_{L,2}, \Delta T_2|H_L) P(T_1) \\ \nonumber
\end{split}
\end{equation}
$$

Where, $\Delta T_2 = T_2 - T_1$.

Now, for the unlensed hypothesis,

$$
\begin{equation}
\begin{split}
P_{astro}(\theta|H_U) &= P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}, m_{1,2}, m_{2,2}, d_{L,2}, \iota_2, d_{L,2}, \Delta T_2|T_1, H_U) P(T_1) \\ \nonumber
\\
&= P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}|T_1, H_U) \\
& P_{astro}(m_{1,2}, m_{2,2}, d_{L,2}, \iota_2, d_{L,2}, \Delta T_2|T_1, H_U) P(T_1) \\ \nonumber
\\
&= P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}|H_U) \\
& P_{astro}(m_{1,2}, m_{2,2}, d_{L,2}, \iota_2, d_{L,2}, \Delta T_2|H_U) P(T_1) \\ \nonumber
\end{split}
\end{equation}
$$

$P(T_1)$ in $P_{astro}(\theta|H_L)$ and $P_{astro}(\theta|H_U)$ cancels out. Similarly, $P(T_1)$ and $P(T_2)$ in $P(\theta|d_1)$ and $P(\theta|d_2)$ cancels out between the numerator and denominator of the Bayes factor.

$P_{pe}(\theta)$: Prior distribution of the parameters which was used in the parameter estimation to get the posterior distribution of the parameters, $P(\theta|d_1)$ and $P(\theta|d_2)$.

$$
\begin{equation}
\begin{split}
P_{pe}(\theta) &= P_{pe}(m_{1}, m_{2}, d_{L}, \iota, t, \alpha, \delta)  \\ \nonumber
&= P_{pe}(m_{1}, m_{2}) P_{pe}(d_{L}) P_{pe}(\iota) P_{pe}(t) P_{pe}(\alpha) P_{pe}(\delta)  \\ \nonumber
\end{split}
\end{equation}
$$

$P_{pe}(t)$ will collapse to constants $P_{pe}(T_1)$ or $P_{pe}(T_2)$ in the numerator and denominator of the Bayes factor. And later, it will cancel out.

For $P_{pe}(\theta)$, I will consider the `bilby` prior for the masses, luminosity distance, inclination angle and sky location. It is summerized in the following table.

| Parameter | Unit                 | Sampling Method        | Range               |
|-----------|----------------------|------------------------|---------------------|
| $d_L$     | Mpc                  | Power-law              | [$d_{L,\text{min}}$,$d_{L,\text{min}}$] |
| $m_1,m_2$ | $\mathcal{M}_{\odot}$ | Uniform   | [$m_{\text{min}}$,$m_{\text{max}}$] |
| $ra$      | Radian               | Uniform                | [0, $2\pi$]         |
| $dec$     | Radian               | Cosine                 | [$-\pi/2$, $\pi/2$] |
| $\iota$   | Radian               | Sine                   | [0, $\pi$]          |
| $\psi$    | Radian               | Uniform                | [0, $\pi$]          |
| $\phi_c$  | Radian               | Uniform                | [0, $2\pi$]         |
| $t_c$     | Second               | Uniform                | [$t_{\text{min}}$, $t_{\text{max}}$] |


**Note**: For numerator the sampling will be done from a combine distribution of data1 and data2.

$P_{astro}(\theta|H_L) = \frac{P_{astro}(\theta|H_L)}{(P(\theta'|d_2) \cup P(\theta'|d_2))}\times (P(\theta'|d_1) \cup P(\theta'|d_2))$

$P_{astro}(\theta|H_L) = \frac{P_{astro}(m_{1}, m_{2}, \alpha, \delta, \iota, d_{L,1}, d_{L,2}, \Delta T_2|H_U)}{P_{comb}(m_{1}, m_{2}, \alpha, \delta, \iota, d_{L,1}, d_{L,2}|d_1, d_2)}P_{comb}(m_{1}, m_{2}, \alpha, \delta, \iota, d_{L,1}, d_{L,2}|d_1, d_2)$

## Bayes factor

$$
\begin{equation}
\begin{split}
\mathcal{B}_U^L &= \frac{ \int d\theta \frac{P(\theta|d_1)}{P_{pe}(\theta)} \frac{P(\theta|d_2)}{P_{pe}(\theta)} P_{astro}(\theta|H_L) }
{\int d\theta \frac{P(\theta|d_1)}{P_{pe}(\theta)} \frac{P(\theta|d_2)}{P_{pe}(\theta)} P_{astro}(\theta|H_U)} \\ \nonumber
=& \Bigg[\\
& \bigg< \\
&\;P(m_{1}, m_{2},\alpha, \delta, \iota, d_{L,1}|T_1, d_1)\times \\
&\;P(m_{1}, m_{2},\alpha, \delta, \iota, d_{L,2}|T_2, d_2) \times \\
&\; P_{astro}(m_{1}, m_{2}, \alpha, \delta, \iota, d_{L,1}, d_{L,2}, \Delta T_2|H_L) \div \\
& \bigg(P_{pe}(m_{1}, m_{2}) P_{pe}(\alpha) P_{pe}(\delta) P_{pe}(\iota) P_{pe}(d_{L})\bigg)^2 \div \\
&\; P_{comb}(m_{1}, m_{2}, \alpha, \delta, \iota, d_{L,1}, d_{L,2}|d_1, d_2) \\
&\bigg>_{ m_{1}, m_{2}, d_{L,1}, \iota, d_{L,2}\in P_{comb}(.|d_1, d_2)} \\
& \Bigg] \div \Bigg[\\
& \bigg<\\
&\; P_{astro}(m_{1,1}, m_{2,1}, \alpha_1, \delta_1, \iota_1, d_{L,1}|H_U)\times \\
&\; P_{astro}(m_{1,2}, m_{2,2}, d_{L,2}, \iota_2, d_{L,2}, \Delta T_2|H_U) \div \\
& \bigg(P_{pe}(m_{1}, m_{2}) P_{pe}(\alpha) P_{pe}(\delta) P_{pe}(\iota) P_{pe}(d_{L})\bigg)^2 \\
&\bigg>_{ m_{1,1}, m_{2,1},\alpha_1, \delta_1, \iota_1, d_{L,1}\in P(.|T_1,d_1),}\\&_{\; m_{1,2}, m_{2,2},\alpha_2, \delta_2, \iota_2, d_{L,2}\in P(.|T_2,d_2)} \\
& \Bigg]
\end{split}
\end{equation}
$$

Where, $<.>$ : average over the samples drawn from the priors.

Multidimensional KDE needs to be constructed both for sampling and pdf evaluation.

Priors:

* $P_{comb}(m_{1}, m_{2}, \alpha, \delta, \iota, d_{L,1}, d_{L,2}|d_1, d_2)$, $P(m_{1}, m_{2},\alpha, \delta, \iota, d_{L,1}|T_1, d_1)$, $P(m_{1}, m_{2},\alpha, \delta, \iota, d_{L,2}|T_2, d_2)$ are the considered priors for sampling the parameters.  

PDF:

<!-- * $P(m_{1}, m_{2},\alpha, \delta, \iota, d_{L,1}|T_1, d_1)$, $P(m_{1}, m_{2},\alpha, \delta, \iota, d_{L,2}|T_2, d_2)$, $P_{astro}(m_{1}, m_{2}, \alpha, \delta, \iota, d_{L,1}, d_{L,2}, \Delta T_2|H_L)$, $$ are the considered pdfs for evaluating the Bayes factor.

* $P(m_{1,1}, m_{2,1}, d_{L,1}, \iota_1|T_1, d_1)$, $P(\alpha_1, \delta_1|T_1, d_1)$, $P(m_{1,2}, m_{2,2}, d_{L,2}, \iota_2|T_2, d_2)$, $P(\alpha_2, \delta_2|T_2, d_2)$ are obtained from the posterior samples of the parameter estimation.

* $P_{pe}(m_{1}, m_{2})$, $P_{pe}(d_{L})$, $P_{pe}(\iota)$, $P_{pe}(\alpha, \delta)$ are the considered priors used for the parameter estimation. -->