# Baseline Hazard

$$
\hat{\Lambda}\{t ; \theta, \alpha\}=\int_0^t \frac{\sum_i \exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)\left\{d N_i(s)-D_i(s)\theta d s-\alpha^\top L d s\right\}}{\sum_i \exp \left\{\int_0^{s-} D_i(u) d u \theta\right\} Y_i(s)}
$$

+ Model misspecification: using cross-fitting to estimate the baseline hazard where I made it a step function
+ Found the previous code used an approximate version of baseline hazard, where I made the denominator a step function, but the simulation is good
+ Approximated the baseline hazard using step function, but simulation is wierd. $\hat{\alpha}$, average of 1000 repetition is not equal to 0.25, whether for $N = 800$ or $N = 1600$.
  + Checked the coding, changed the baseline hazard into the true function, $0.25t$, the estimation is accurate.
  + There are some trends, since when I thought of the proof that the convergence may be related to the $max_k\{T^{(k+1)} - T^{(k)}\}$, I changed the grid from 0.01 to the original value, that the time points estimation was conducted is changed from about 200 to 800 (for $N = 800$). the estimation is better, but it is still not equal to $0.25$.

## Proof details

$$
\hat{\Lambda}\{t ; \theta, \alpha\} - \Delta(t)\bar{\Lambda}(t;\theta, \alpha) = \hat{\Lambda}\{t_m ; \theta, \alpha\} - \Delta(t)\bar{\Lambda}(t_m;\theta, \alpha) - \left[\bar{\Lambda}(t;\theta, \alpha) - \bar{\Lambda}(t_m;\theta, \alpha)\right]
$$

$$
\left|\hat{\Lambda}\{\tau ; \theta, \alpha\} - \bar{\Lambda}(\tau;\theta, \alpha)\right|
$$

$$
\int_0^\tau \left|\hat{\Lambda}\{t ; \theta, \alpha\} - \bar{\Lambda}(t;\theta, \alpha)\right|dte^{\theta_{\max} t}
$$

## Computation details

$$
-\sum_i\int_0^\tau Z_i^c\exp \left\{\int_0^{t-} D_i(u) d u \theta\right\} Y_i(t)d\hat{\Lambda}(t;\theta, \alpha)
$$