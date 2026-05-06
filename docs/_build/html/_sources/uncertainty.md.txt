# Taking care of uncertainty

## SIS case

$$ 
\begin{split}
P(z_s)
= \frac{1}{{\cal N}_{\rm U}}\,
\frac{R_{\rm U}(z_s)}{1+z_s}\,
\frac{dV_c}{dz_s} \, ,
\end{split}
$$

$$ 
\begin{split}
{\cal N}_{\rm L}
&= {\cal N}_{\rm U}\int P({\rm SL}\mid z_s)\,
\frac{1}{{\cal N}_{\rm U}}\,
\frac{R_{\rm U}(z_s)}{1+z_s}\,
\frac{dV_c}{dz_s}\,dz_s \, \\
&= {\cal N}_{\rm U}\int P({\rm SL}\mid \sigma, z_L, z_s)\,
P(\sigma\mid z_L) \, P(z_L\mid z_s)
P(z_s) \, dz_s \\
&= {\cal N}_{\rm U}\int \frac{\sigma^{\rm SIS}_{\rm SL}}{4\pi}\,
P(\sigma\mid z_L) \, P(z_L\mid z_s)
P(z_s) \, dz_s \\
\end{split}
$$

New Normalization :

$$ 
\begin{split}
{\cal N}_{\rm L}
= \frac{1}{{\cal N}_{\rm U}}\,
\left\langle
\frac{\sigma^{\rm SIS}_{\rm SL}}{4\pi}\,
\frac{P_{\rm new}(\sigma\mid z_L) \, P_{\rm new}(z_L\mid z_s)
P_{\rm new}(z_s)}{P_{\rm old}(\sigma\mid z_L) \, P_{\rm old}(z_L\mid z_s)
P_{\rm old}(z_s)}
\right\rangle_{\substack{
\sigma\sim P_{\rm old}(\sigma\mid z_L) \\
z_L\sim P_{\rm old}(z_L\mid z_s) \\
z_s\sim P_{\rm old}(z_s)
}} 
\end{split}
$$

