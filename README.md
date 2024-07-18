# Computing-Large-Deviation-Functionals-under-constraints

In this repository, there are the C-scripts we use to compute the large deviation functions in "On Magnetic Models in Wavefunction Ensembles", published in [Entropy 2023, 25, 564](https://www.mdpi.com/1099-4300/25/4/564) covering theorem 3.
The corresponding preprint version "On Schrödingerist Quantum Thermodynamics" is available on [arXiv:2208.07688](http://arxiv.org/abs/2208.07688), where this is referred to as theorem 2.

We suggest to refer to the arXiv version since we prefer its organization and there are some corrections to conceptual sentences (proofs are unchanged).
Here we refer to the arXiv version for the enumeration.

The repository contains all the necessary header files to run the main script Large_Deviations_SCW.c, which computes the large deviation functionals (157) (or (161)) and (162). These are rate functionals of NOT-identically distributed independent Gaussian random variables; therefore, no explicit formula is possible. Moreover, the computed probabilities in (146) depend on a parameter \(x\), which sets the set of wavefunctions over which the probability integral (144) is computed. The computational techniques are developed for Gaussian random variables because of how wavefunction ensembles are defined, but they do not restrict to the Gaussian case.





## Computational Appendix F

**Figures:**
- **Figure 1**: [Geometry of the Domain D](figures/figure_1.png) - This figure shows the geometry of the domain D.See the papers for how to obtain it.
- **Figure 2**: [Domain and Constraint Set, by Simulation](figures/SCW_D_and_G_plot.pdf) - This figure illustrates the domain and constraint set by sampling.
- **Figure 3**: [I-functions of x](figures/SCWIplot.pdf) - This figure plots the I-functions of x.



**Description:**

After implementing the functions appearing in Math Appendix D in a program, we employed a sampling scheme to locate the constraint set. The idea is to choose at random a ray emanating from point $'P'$ of [figure_1](figures/figure_1.png), then a point on the ray staying within the domain $D$. (We also tried a ray emanating from the origin, which gave similar results.)

By experiment, we discovered that for $x$ near $0$ near zero, the region $G$ shrank to nearly a line close to line $A$ in [figure_1](figures/figure_1.png), while the constraint set for generating $I_1$ lays almost at the left endpoint, $'P'$. Hence, we needed sampling schemes that could be biased to prefer points near the endpoints of a given interval $[a, b]$ of the real axis. Such schemes are given by:

- To bias near $b$, choose a point $s$ by the scheme:
  
  $$s = a + \log(z \eta u + 1)/\eta$$
  
  $$z = (e^{\eta(b-a)} - 1)/\eta$$

- And to bias near $a$:

  $$s = a - \log(1 - z \eta u \exp(-\eta(b - a)))/\eta$$

  with the same expression for $z$.

Plugging in a uniform random variable (produced by the system RNG) for $u$ yields the same scheme for $s \in [a, b]$. With $\eta = 0$ the sampling is uniform on the interval; large positive $\eta$ yields bias. We will write: $s = biased-sample\(a, b\)$ for the random sample.

Our sampling scheme in the region $D$ was the following:

1. Choose $\alpha: \alpha = biased-sample\(-x/(1 + x), -x/(1 - \epsilon)\)$;
2. Choose $\theta_1: \theta_1 = biased-sample\(-1/(2x), 5/(1 - x)\)$;
3. Let $\theta_2 = \alpha [\theta_1 + 1/(2x)] $.

If the selected $(\theta_1, \theta_2)$ pass all the tests to lie in $D$, accept the values; otherwise, reject them.

To locate the constraint set, we sampled many pairs $(\theta_1, \theta_2)$ as above and evaluated the two partial derivatives of $c$, keeping and plotting the points that had the required signs. The results—see [figure_2](figures/figure_2.png); parameters in the figure were $x = 0.7$ and $\epsilon = 0.3$—show that $G$ lies in the upper half of the plane, and is disjoint from the horizontal axis. [Figure_3](figures/figure_1.png) shows curves of the I-functions obtained by sampling for a few values of $x$ ($0.1$ to $0.7$, in increments of $0.1$). We used $10$ billion samples at each $x$-value.

Unfortunately, we were unable to estimate the I-functions for $x < 0.1$ because few or no sampled points fell in the constraint set (even with various choices of bias), for reason indicated earlier.
