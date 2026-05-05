# Tutorial: Computing the Source-Frame BBH Merger-Rate Density from SFR + Exponential Time Delay

This tutorial explains how to compute the **source-frame volumetric merger-rate density** of binary black holes (BBHs),

$$ 
\begin{equation}
\begin{split}
R_{\rm m}(z)\equiv \frac{dN_{\rm merge}}{dV_c\,dt_s},
\end{split}
\end{equation} 
$$

starting from the **cosmic star-formation rate density** (SFR) and a **time-delay distribution** between binary formation and merger.

The goal is that a Master's student can:
1. understand the mathematics, and  
2. implement it correctly in Python.

We follow the analytical formulation used in *Vitale et al., "Measuring the star formation rate with gravitational waves from binary black holes"*.

---

## 0) What you will compute (the target quantity)

### Source-frame merger-rate density

The source-frame merger-rate density is defined as

$$ 
\begin{equation}
\begin{split}
R_{\rm m}(z)\equiv \frac{dN_{\rm merge}}{dV_c\,dt_s}
\quad [{\rm Gpc^{-3}\,yr^{-1}}].
\end{split}
\end{equation} 
$$

- $dV_c$: comoving volume element  
- $dt_s$: **source-frame** time interval  
- $z$: merger redshift  

> This is an **astrophysical merger-rate density**, not a detector yield.

---

## 1) Ingredients of the model (physical picture)

BBHs in the field form from massive stellar binaries. The merger-rate density at redshift $z$ depends on:

1. **how many stars formed at earlier cosmic times** (the SFR history),  
2. **how efficiently that star formation produces BBH progenitors** (metallicity dependent),  
3. **how long it takes from formation to merger** (time delay).  

Mathematically, this becomes a **convolution** of the formation rate with a delay-time distribution.

---

## 2) Cosmology: mapping redshift $\leftrightarrow$ cosmic time

To connect **formation redshift** $z_f$ to **merger redshift** $z_m$, we use the cosmic time $t(z)$, interpreted as the **age of the Universe at redshift $z$**.

### 2.1 Expansion rate

Assuming flat $\Lambda$CDM, the Hubble rate is

$$ 
\begin{equation}
\begin{split}
H(z)=H_0\,E(z),
\qquad
E(z)=\sqrt{\Omega_m(1+z)^3+\Omega_\Lambda}.
\end{split}
\end{equation} 
$$

### 2.2 Differential time–redshift relation

Cosmic time decreases with redshift as

$$ 
\begin{equation}
\begin{split}
\frac{dt}{dz}=-\frac{1}{(1+z)H(z)}.
\end{split}
\end{equation} 
$$

We will often use the positive quantity

$$ 
\begin{equation}
\begin{split}
\left|\frac{dt}{dz}\right|=\frac{1}{(1+z)H(z)}.
\end{split}
\end{equation} 
$$

### 2.3 Delay time

If a BBH forms at redshift $z_f$ and merges at redshift $z_m$ (with $z_f>z_m$), the delay time is

$$ 
\begin{equation}
\begin{split}
t_d(z_m,z_f)=t(z_m)-t(z_f)\ge 0.
\end{split}
\end{equation} 
$$

This definition is unambiguous when $t(z)$ is the **cosmic age**.

---

## 3) Star-formation rate density $\psi(z)$

We use the Madau–Dickinson functional form

$$ 
\begin{equation}
\begin{split}
\psi(z)=\psi_0\,
\frac{(1+z)^\alpha}{1+\left(\frac{1+z}{C}\right)^\beta},
\end{split}
\end{equation} 
$$

with typical parameters:
- $\alpha=2.7$
- $\beta=5.6$
- $C=2.9$
- $\psi_0=0.015~M_\odot\,{\rm Mpc^{-3}\,yr^{-1}}$

Units:
- $\psi(z)$ is in $M_\odot\,{\rm Mpc^{-3}\,yr^{-1}}$.

---

## 4) Metallicity-dependent efficiency $\eta(z)$

BBHs form more efficiently at **low metallicity**, because stellar winds are weaker and massive stars retain more mass. We encode this using an efficiency factor $\eta(z)\in[0,1]$, interpreted as:

- $\eta(z)$ is the fraction of star formation occurring at metallicity below $0.1Z_\odot$.

### 4.1 Metallicity distribution (log-normal model)

Assume

$$ 
\begin{equation}
\begin{split}
\log_{10}\!\left(\frac{Z}{Z_\odot}\right)
\sim
\mathcal{N}\!\left(
\log_{10}\!\left(\frac{Z_{\rm mean}(z)}{Z_\odot}\right),
\sigma_Z^2
\right),
\end{split}
\end{equation} 
$$

with $\sigma_Z=0.5$ dex.

### 4.2 Efficiency as a cumulative probability

Then

$$ 
\begin{equation}
\begin{split}
\eta(z)=P(Z<0.1Z_\odot)
=
\frac{1}{2}\left[
1+\operatorname{erf}\left(
\frac{\log_{10}(0.1)-\log_{10}\!\left(\frac{Z_{\rm mean}(z)}{Z_\odot}\right)}
{\sqrt{2}\sigma_Z}
\right)\right].
\end{split}
\end{equation} 
$$

---

## 4.3 Mean metallicity evolution $Z_{\rm mean}(z)$ (Belczynski-style)

Vitale points to Belczynski et al. (2016) for $Z_{\rm mean}(z)$. The commonly used prescription is

$$ 
\begin{equation}
\begin{split}
\log_{10}\!\left(\frac{Z_{\rm mean}(z)}{Z_\odot}\right)
=
0.5+\log_{10}\!\left[
\frac{y(1-R)}{\rho_b}
\int_{z}^{20}
\frac{97.8\times 10^{10}\,\psi(z')}
{H_0\,E(z')\,(1+z')}
\,dz'
\right].
\end{split}
\end{equation} 
$$

with constants:
- $R=0.27$ (return fraction),
- $y=0.019$ (metal yield),
- $\rho_b=2.77\times 10^{11}\,\Omega_b\,h_0^2~M_\odot~{\rm Mpc^{-3}}$ (baryon density).

### Practical note on numerical implementation

Instead of hard-coding the factor $97.8\times 10^{10}$, the numerically safest approach is to compute the time element directly:

$$ 
\begin{equation}
\begin{split}
\frac{dz'}{H_0E(z')(1+z')}
=
\left|\frac{dt}{dz'}\right|\,dz',
\end{split}
\end{equation} 
$$

so the metallicity integral can be evaluated with consistent units as

$$ 
\begin{equation}
\begin{split}
\int_{z}^{z_{\max}}\psi(z')\left|\frac{dt}{dz'}\right|dz'.
\end{split}
\end{equation} 
$$

This preserves the same physics while avoiding unit-normalization mistakes.

---

## 5) BBH formation-rate density $R_f(z)$

Assume BBH formation follows SFR weighted by low-metallicity efficiency:

$$ 
\begin{equation}
\begin{split}
R_f(z)=A\,\eta(z)\,\psi(z).
\end{split}
\end{equation} 
$$

Here $A$ is an unknown proportionality constant, fixed by local-rate normalization.

---

## 6) Exponential time-delay distribution

We use an exponential time-delay kernel:

$$ 
\begin{equation}
\begin{split}
p(t_d|\tau)=\frac{1}{\tau}\exp(-t_d/\tau),
\qquad t_d\ge 0.
\end{split}
\end{equation} 
$$

It is normalized automatically:

$$ 
\begin{equation}
\begin{split}
\int_{0}^{\infty}p(t_d|\tau)\,dt_d=1.
\end{split}
\end{equation} 
$$

---

## 7) Merger-rate density $R_m(z)$: the key equation

A BBH that merges at $z_m$ must have formed earlier at $z_f\ge z_m$. The merger-rate density is

$$ 
\begin{equation}
\begin{split}
R_m(z_m)
=
\int_{z_m}^{\infty}
dz_f\,
\left|\frac{dt}{dz_f}\right|
R_f(z_f)\,
p\!\left(t_d(z_m,z_f)\mid \tau\right).
\end{split}
\end{equation} 
$$

Substituting
$R_f(z_f)=A\,\eta(z_f)\psi(z_f)$ and the exponential kernel gives

$$ 
\begin{equation}
\begin{split}
R_m(z_m)
=
A\int_{z_m}^{\infty}
dz_f\,
\left|\frac{dt}{dz_f}\right|
\eta(z_f)\psi(z_f)\,
\frac{1}{\tau}\exp\!\left[-\frac{t(z_m)-t(z_f)}{\tau}\right].
\end{split}
\end{equation} 
$$

This is the **main formula** to implement.

---

## 8) Normalization using the local BBH merger rate (GWTC-4)

The unknown constant $A$ is fixed by imposing:

$$ 
\begin{equation}
\begin{split}
R_m(0)=R_0^{\rm BBH}.
\end{split}
\end{equation} 
$$

A representative local BBH rate from GWTC-4 is

$$ 
\begin{equation}
\begin{split}
R_0^{\rm BBH}\approx 19~{\rm Gpc^{-3}\,yr^{-1}}.
\end{split}
\end{equation} 
$$

### 8.1 Compute an unnormalized shape integral

Define the “shape” integral (with $A=1$):

$$ 
\begin{equation}
\begin{split}
I_0(\tau)=
\int_{0}^{\infty}
dz_f\,
\left|\frac{dt}{dz_f}\right|
\eta(z_f)\psi(z_f)\,
\frac{1}{\tau}\exp\!\left[-\frac{t(0)-t(z_f)}{\tau}\right].
\end{split}
\end{equation} 
$$

### 8.2 Solve for $A$

Then

$$ 
\begin{equation}
\begin{split}
A=\frac{R_0^{\rm BBH}}{I_0(\tau)}.
\end{split}
\end{equation} 
$$

Now $R_m(z)$ is fully specified.

---

## 9) Step-by-step computational algorithm (what to code)

### Step 1: Choose cosmology
Choose $(H_0,\Omega_m,\Omega_\Lambda)$, or use a built-in cosmology such as Planck18 in `astropy`.

### Step 2: Create a redshift grid
Choose $z_{\max}$ (e.g. $10$ or $20$) and construct
$z=[0,\Delta z,\dots,z_{\max}]$.

### Step 3: Compute $t(z)$ and $\left|\frac{dt}{dz}\right|$
Evaluate $t(z)$ (cosmic age) and

$$ 
\begin{equation}
\begin{split}
\left|\frac{dt}{dz}\right|=\frac{1}{(1+z)H(z)}.
\end{split}
\end{equation} 
$$

### Step 4: Compute $\psi(z)$
Evaluate the Madau–Dickinson SFR $\psi(z)$ on the grid.

### Step 5: Compute $\eta(z)$
1. compute $Z_{\rm mean}(z)$ using the Belczynski-style integral,
2. compute $\eta(z)$ using the error-function CDF formula.

### Step 6: Choose $\tau$
Pick a time-delay scale (example: $\tau=1~{\rm Gyr}$).

### Step 7: Compute the unnormalized $R_m(z)$
For each merger redshift $z_m$, evaluate the integral over $z_f\in[z_m,z_{\max}]$.

### Step 8: Normalize using GWTC-4
Enforce $R_m(0)=R_0^{\rm BBH}$ to solve for $A$.

---

## 10) Final output and detector-frame conversion

### Source-frame output
You obtain $R_m(z)$ in ${\rm Gpc^{-3}\,yr^{-1}}$.

### Detector-frame redshift rate density (all sky)

The detector-frame redshift distribution (events per observed time per redshift) is

$$ 
\begin{equation}
\begin{split}
\frac{d\dot N}{dz}
=
\frac{R_m(z)}{1+z}\,\frac{dV_c}{dz}.
\end{split}
\end{equation} 
$$

Here:
- $(1+z)^{-1}$ converts source-frame time to observer-frame time,
- $\frac{dV_c}{dz}$ is the comoving volume element over the full sky.

---

## 11) Practical notes (numerical and scientific)

1. **Upper limit of the formation redshift integral**  
   Replace $\infty$ by a finite $z_{\max}$. In practice $z_{\max}\sim 10$–$20$ is sufficient.

2. **Unit consistency**  
   If $\tau$ is in Gyr, then $t(z)$ must be in Gyr.  
   Using `astropy` avoids common unit conversion errors.

3. **Interpretation of $\tau$**  
   Small $\tau$: mergers closely track the SFR peak.  
   Large $\tau$: the merger rate shifts to lower redshift.

4. **Normalization**  
   Fixing $R_m(0)=R_0^{\rm BBH}$ sets the overall amplitude; the rest is the redshift evolution.

---

## 12) Summary of the full mathematical model (one box)

$$ 
\begin{equation}
\begin{split}
R_m(z_m)
=
A\int_{z_m}^{z_{\max}}
dz_f\,
\left|\frac{dt}{dz_f}\right|
\eta(z_f)\psi(z_f)\,
\frac{1}{\tau}\exp\!\left[-\frac{t(z_m)-t(z_f)}{\tau}\right].
\end{split}
\end{equation} 
$$

with

$$ 
\begin{equation}
\begin{split}
A=\frac{R_0^{\rm BBH}}{R_m(0)\big|_{A=1}},
\qquad
R_0^{\rm BBH}\approx 19~{\rm Gpc^{-3}\,yr^{-1}}.
\end{split}
\end{equation} 
$$

This completes the theory needed to implement the BBH merger-rate density in Python.
