# Lennard-Jones Ground States

### Problem Description

The 2D system of interest consists of $N$ particles within a square box under periodic boundary conditions interacting with the Lennard-Jones potential,
$$
V_{\text{LJ}}(r)=4\epsilon[(\frac{\sigma}{r})^{12}-(\frac{\sigma}{r})^6],
$$
whose graph looks like[^1]

![img](https://upload.wikimedia.org/wikipedia/commons/e/e7/Graph_of_Lennard-Jones_potential.png)

with the minimum located at $r_m=2^{1/6}\sigma$. The goal is to find the ground state of the system.

The ground state is approached by cooling down the system. For any specific particle configuration characterized by $r^N$, the total energy $U$ is given by the pair decomposable form[^2]
$$
U(r^N)=\sum_{i>j=1}^NV_{\text{LJ}}(|\boldsymbol{r}_i-\boldsymbol{r}_j|).
$$

Beginning with a random initial particle configuration $r^N_0$ and a large initial temperature $T_0$, we try moving one particle by a random displacement each time and calculate the associated energy change $\Delta E$. If $\Delta E < 0$, we always accept the move; else we accept the move with a Boltzmann factor probability, i.e., $\exp(-\Delta E / T)$, where $T$ is the temperature. The temperature $T$ should be gradually lowered after each equilibrium at a particular temperature. Eventually, as $T$ approaches zero, fewer and fewer uphill moves will be accepted, and the ground state should be achieved.

One specific ground state can be utilized to verify the code. With the box length set as $1$, at the large $N$ limit, the triangular lattice constant should be
$$
a = 1 / \left(\sqrt{\frac{\sqrt{3}}{2}N}\right) = \sqrt{\frac{2}{\sqrt{3}N}}.
$$
If $a = r_m$, i.e.
$$
\sigma = \frac{2^{1/3}3^{-1/4}}{\sqrt{N}},
$$
then the triangular lattice should definitely be the ground state configuration, which can be compared with the numerical result.


[^1]: By <a href="https://en.wikipedia.org/wiki/User:TimeStep89" class="extiw" title="w:User:TimeStep89">TimeStep89</a> - I computed the data and created the plot with software consistent with other plots for the "Lennard-Jones potential" Wikipedia article., <a href="https://creativecommons.org/licenses/by/4.0" title="Creative Commons Attribution 4.0">CC BY 4.0</a>, <a href="https://commons.wikimedia.org/w/index.php?curid=133883513">Link</a>
[^2]: Which is only an approximation. (David Chandler， *Introduction To Modern Statistical Mechanics*)
