Below is a clean “from-first-principles” explanation of what the handwritten note is doing, with the key goal:

> infer the hyper-parameter $\tau$ that controls the source redshift distribution $p(z_s\mid\tau)$, using multiple GW datasets $\{\vec d_i\}_{i=1}^{N}$, including the selection effect (i.e. conditioning on detection).

The notation follows the note, where the subscript “0” means Bilby PE done with a simplified prior $p_0(\theta)$.

---

## 1) Event-level inference (Bilby PE): what $p_0(\theta\mid d)$ means

For one event with data $\vec d$, Bilby returns posterior samples from

$$
p_0(\theta \mid \vec d)
=\frac{p_0(\vec d\mid \theta)\,p_0(\theta)}{p_0(\vec d)}\,,
$$

where

- $\theta$ = event parameters (masses, spins, distance or $z$, sky location, etc.)
- $p_0(\theta)$ = the Bilby sampling prior (simplified)
- $p_0(\vec d)$ = the evidence under that prior

The evidence is

$$
p_0(\vec d)=\int p_0(\vec d\mid\theta)\,p_0(\theta)\,d\theta.
$$

This “0” posterior is not yet a population inference result; it is just the per-event posterior under the PE prior.

---

## 2) Population model: what you want to infer

You want a population model in which the source redshift distribution is described by a hyper-parameter $\tau$ :

$$
p(z_s\mid\tau)
$$

In practice, the population model is over all event parameters :

$$
p(\theta\mid\tau),
$$

and $\theta$ contains $z_s$ (and everything else).

### Common factorization 

Often one uses

$$
p(\theta\mid\tau)=p(z_s\mid\tau)\,p(\lambda)\,,
$$

where $\lambda$ denotes the “other parameters” (masses/spins/orientation), possibly with their own hyper-parameters (but you can keep them fixed if the only target is $z_s$).

---

## 3) Where selection effects enter (the conceptual point)

The catalog you analyze is a set of detected events, not “all events that occurred”.

So the correct likelihood is always conditioned on detection:

$$
p(\vec d \mid \tau, {\rm det}).
$$

This introduces a normalization by the detection efficiency of the population, written as $\alpha(\tau)$ (some notes call it $\sigma(\tau)$ or $\beta(\tau)$).

---

## 4) Single-event likelihood with selection (the role of $\sigma$)

A convenient way to write the single-event contribution is

$$
p(\vec d_i \mid \tau, {\rm det}) =
\frac{
\int p(\vec d_i \mid \theta)\,p(\theta\mid\tau)\,d\theta
}{
\alpha(\tau)
}\,,
$$

where the selection function is

$$
\alpha(\tau)=
p({\rm det}\mid\tau)=
\int p({\rm det}\mid\theta)\,p(\theta\mid\tau)\,d\theta.
$$

Interpretation:

- Numerator: “support of the data under the population prior”
- Denominator: “fraction of the population that would be detected”

---

## 5) Multiple-event hierarchical posterior for $\tau$

For $N$ independent detected events $\{\vec d_i\}$, the hierarchical posterior is

$$
p(\tau \mid \{\vec d_i\}, {\rm det})
\propto
p(\tau)\,
\prod_{i=1}^{N}
p(\vec d_i \mid \tau, {\rm det}).
$$

Substituting the selection-corrected event likelihood gives the catalog likelihood

$$
\mathcal{L}(\tau)=
\prod_{i=1}^N
\frac{
\int p(\vec d_i \mid \theta)\,p(\theta\mid\tau)\,d\theta
}{
\alpha(\tau)
}.
$$

---

## 6) The key computational trick: reuse Bilby PE samples

Bilby gives

$$
p_0(\theta \mid \vec d_i)=
\frac{p(\vec d_i\mid \theta)\,p_0(\theta)}{p_0(\vec d_i)}.
$$

Rearranging,

$$
p(\vec d_i\mid \theta)=
\frac{p_0(\theta \mid \vec d_i)\,p_0(\vec d_i)}{p_0(\theta)}.
$$

Plug into the numerator integral:

$$
\int p(\vec d_i\mid \theta)\,p(\theta\mid\tau)\,d\theta=
p_0(\vec d_i)\,
\int p_0(\theta\mid\vec d_i)\,
\frac{p(\theta\mid\tau)}{p_0(\theta)}\,d\theta.
$$

Approximate the integral using posterior samples $\{\theta_{ik}\}_{k=1}^{n_i}\sim p_0(\theta\mid\vec d_i)$:

$$
\int p_0(\theta\mid\vec d_i)\,
\frac{p(\theta\mid\tau)}{p_0(\theta)}\,d\theta
\approx
\frac{1}{n_i}\sum_{k=1}^{n_i}
\frac{p(\theta_{ik}\mid\tau)}{p_0(\theta_{ik})}.
$$

Therefore, up to $\tau$-independent constants,

$$
p(\vec d_i \mid \tau, {\rm det})
\propto
\frac{1}{\alpha(\tau)}
\left\langle
\frac{p(\theta\mid\tau)}{p_0(\theta)}
\right\rangle_{\theta\sim p_0(\theta\mid \vec d_i)}.
$$

The Bilby evidences $p_0(\vec d_i)$ drop out because they do not depend on $\tau$.

---

## 7) Special case: redshift-only hyper-inference

If only the redshift part changes and everything else matches the PE prior, then

$$
\frac{p(\theta\mid\tau)}{p_0(\theta)}
\approx
\frac{p(z_s\mid\tau)}{p_0(z_s)}.
$$

So the event contribution becomes

$$
p(\vec d_i \mid \tau, {\rm det})
\propto
\frac{1}{\alpha(\tau)}
\left\langle
\frac{p(z_s\mid\tau)}{p_0(z_s)}
\right\rangle_{z_s\sim p_0(z_s\mid \vec d_i)}.
$$

and the catalog likelihood is

$$
\mathcal{L}(\tau)
\propto
\frac{1}{\alpha(\tau)^N}
\prod_{i=1}^N
\left\langle
\frac{p(z_s\mid\tau)}{p_0(z_s)}
\right\rangle_{z_s\sim p_0(z_s\mid \vec d_i)}.
$$

---

## 8) How to compute $\alpha(\tau)$ in practice (injections)

The selection function is

$$
\alpha(\tau)=\int p({\rm det}\mid\theta)\,p(\theta\mid\tau)\,d\theta.
$$

Monte Carlo with injections $\theta_j\sim p_{\rm inj}(\theta)$ gives

$$
\alpha(\tau)\approx
\frac{1}{N_{\rm inj}}
\sum_{j=1}^{N_{\rm inj}}
p({\rm det}\mid\theta_j)\,
\frac{p(\theta_j\mid\tau)}{p_{\rm inj}(\theta_j)}.
$$

---

## 9) Final one-line formula (what you implement)

The hierarchical posterior is

$$
p(\tau \mid \{\vec d_i\}, {\rm det})
\propto
p(\tau)\,
\frac{1}{\alpha(\tau)^N}
\prod_{i=1}^N
\left[
\frac{1}{n_i}\sum_{k=1}^{n_i}
\frac{p(\theta_{ik}\mid\tau)}{p_0(\theta_{ik})}
\right].
$$

This matches the structure in the handwritten note: reweight each event’s PE posterior samples by the ratio of the population prior to the PE prior, then divide by the selection function.
