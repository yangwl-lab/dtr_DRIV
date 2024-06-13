$$
\hat{\Lambda}\{t ; \theta, \alpha\}=\int_0^t \frac{\sum_i \exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)\left\{d N_i(s)-D_i(s)\theta d s-\alpha^\top L d s\right\}}{\sum_i \exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)}
$$


$$
\hat{\Lambda}\{t_m ; \theta, \alpha\} - \hat{\Lambda}\{t_{m-1} ; \theta, \alpha\}=\int_{t_{m-1}}^{t_m} \frac{\sum_i \exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)\left\{d N_i(s)-D_i(s)\theta d s-\alpha^\top L d s\right\}}{\sum_i \exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)}
$$

## Notations

$$
K_i = \exp\left(\int_0^{t_{m-1}}D_i(u)\theta du\right)
$$

$$
A = \sum_j K_j \{1 - D_j(s)\}Y_j(s)
$$

$$
B = \sum_j K_j  D_j(s)Y_j(s)
$$

In the interval $[t_{m-1}, t_m]$, All of $K_i$, $A$ and $B$ are constant.

## Integral

$$
\begin{aligned}
    &\hat{\Lambda}\{t_m ; \theta, \alpha\} - \hat{\Lambda}\{t_{m-1} ; \theta, \alpha\}\\
    =& \int_{t_{m-1}}^{t_m}\frac{\sum_i \exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)d N_i(s)}{\sum_i \exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)}\\
     &-\sum_i\int_{t_{m-1}}^{t_m}\frac{K_i Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\} \exp\left(\int_{t_{m-1}}^{s-}D_i(u)\theta du\right) ds}{\exp\left\{ \theta (s - t_{m-1})\right\}B + A}
\end{aligned}
$$

For the second item, it can be decomposed into
$$
\begin{aligned}
    &\int_{t_{m-1}}^{t_m}\frac{K_i Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}D_i(s) \exp\left\{ \theta (s - t_{m-1})\right\} ds}{\exp\left\{ \theta (s - t_{m-1})\right\}B + A}\\
    &+\int_{t_{m-1}}^{t_m}\frac{K_i Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}\{1 - D_i(s)\}  ds}{\exp\left\{ \theta (s - t_{m-1})\right\}B + A}\\
    =& \frac{K_i Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}D_i(s)}{\theta B} \log \left[\frac{\exp\left\{ \theta (t_m - t_{m-1})\right\}B + A}{B+A}\right]\\
    &+\frac{K_i Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}\{1-D_i(s)\}}{\theta A} \log\left[ \frac{B+A}{B + A\exp\left\{ \theta (t_{m-1} - t_{m})\right\}}\right]
\end{aligned}
$$

Find the derivative at x
$$
\begin{aligned}
    &\frac{\partial \left[\hat{\Lambda}\{t_m ; \theta, \alpha\} - \hat{\Lambda}\{t_{m-1} ; \theta, \alpha\}\right]}{\partial \theta}\\
    =& \sum_i\int_{t_{m-1}}^{t_m}\frac{ \left\{\int_0^{s-} D_i(u) d u\right\}\exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)\left[\sum_j \exp \left\{\int_0^{s-} D_j(u) d u \theta\right\} Y_j(s)\right] - \exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)\left[\sum_j \left\{\int_0^{s-} D_j(u) d u\right\}\exp \left\{\int_0^{s-} D_j(u) d u \theta\right\} Y_j(s)\right]d N_i(s)}{\left[\sum_j \exp \left\{\int_0^{s-} D_j(u) d u \theta\right\} Y_j(s)\right]^2}\\
    &+ \frac{\left[K_i' Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}D_i(s) + K_iY_iD_i(s)\right]\theta B -K_i Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}D_i(s)(B + \theta B')}{(\theta B)^2}\log \left[\frac{\exp\left\{ \theta (t_m - t_{m-1})\right\}B + A}{B+A}\right]\\
    &+ \frac{K_i Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}D_i(s)}{\theta B} \times\frac{B+A}{\exp\left\{ \theta (t_m - t_{m-1})\right\}B + A}\times \frac{\left[(t_m - t_{m-1})\exp\left\{ \theta (t_m - t_{m-1})\right\}B + \exp\left\{ \theta (t_m - t_{m-1})\right\}B' + A\right](A+B) - \left[\exp\left\{ \theta (t_m - t_{m-1})\right\}B + A\right] (B' + A')}{(B+A)^2}\\
    &+\frac{\left[K_i' Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}\{1-D_i(s)\} \right]\theta A -K_i Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}D_i(s)(A + \theta A')}{(\theta A)^2}\log\left[ \frac{B+A}{B + A\exp\left\{ \theta (t_{m-1} - t_{m})\right\}}\right]\\
    &+\frac{K_i Y_i(s)\left\{D_i(s)\theta + \alpha^\top L_i\right\}\{1-D_i(s)\}}{\theta A}\times \frac{B + A\exp\left\{ \theta (t_{m-1} - t_{m})\right\}}{B+A}\times \frac{(B'+A')\left[B + A\exp\left\{ \theta (t_{m-1} - t_{m})\right\}\right] - (B+A)\left[B' + A'\exp\left\{ \theta (t_{m-1} - t_{m})\right\} + (t_{m-1} - t_m)A\exp\left\{ \theta (t_{m-1} - t_{m})\right\}\right]}{\left[B + A\exp\left\{ \theta (t_{m-1} - t_{m})\right\}\right]^2}
\end{aligned}
$$



## Another way

$$
d\hat{\Lambda}(t_m;\theta, \alpha) = \frac{\sum_i \exp \left\{\int_0^{t_m-} D_i(u) d u \theta\right\} Y_i(t_m)\left\{d N_i(t_m)-D_i(t_m)\theta dS(t_m)-\alpha^\top L dS(t_m)\right\}}{\sum_i \exp \left\{\int_0^{t_m-} D_i(u) d u \theta\right\} Y_i(t_m)}
$$

where,
$$
dS(t_m) = t_m - t_{m-1}
$$