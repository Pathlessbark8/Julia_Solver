# Adjoint Equations for Maxwell's Equations

## Frequency Domain
Basic equations:

$\nabla \times E = -j\omega\mu_r\mu_0H$

$\nabla \times H = j\omega\varepsilon_r\varepsilon_0E + \sigma E - J$

Altered:

$\nabla \times\frac{1}{\mu_r}\nabla \times E = -j\omega\mu_0\nabla \times H$

$-j\omega\mu_0\nabla \times H = +\omega^2\varepsilon_r\mu_0\varepsilon_0E -j\omega\mu_0 \sigma E +j\omega\mu_0 J$

Second order:

$\nabla \times\frac{1}{\mu_r}\nabla \times E -\omega^2\varepsilon_r\mu_0\varepsilon_0E +j\omega\mu_0 \sigma E = j\omega\mu_0 J$

Some objective function:

$F(\varepsilon_r, \mu_r, E)$

Derivatives:

$
\frac{dF(\varepsilon_r, \mu_r, E)}{d\varepsilon_r}=
\frac{\partial F(\varepsilon_r, \mu_r, E)}{\partial\varepsilon_r} + 
\frac{\partial F(\varepsilon_r, \mu_r, E)}{\partial E}
\frac{dE}{d\varepsilon_r}
$

Implicit derivative:

$
L[E; \varepsilon_r, \mu_r, ]=\nabla \times\frac{1}{\mu_r}\nabla \times E -\omega^2\varepsilon_r\mu_0\varepsilon_0E +j\omega\mu_0 \sigma E - j\omega\mu_0 J
$

$
\frac{d}{d\varepsilon_r}L[ E; \varepsilon_r, \mu_r] = 
\frac{\partial L}{\partial \varepsilon_r} + 
\frac{\partial L}{\partial E}\frac{dE}{d\varepsilon_r}=0
$

$
\frac{\partial L}{\partial \varepsilon_r} = -\omega^2\mu_0\varepsilon_0E
$

$
\frac{\partial L}{\partial E}\frac{\partial E}{\partial \varepsilon_r} = L[\frac{\partial E}{\partial \varepsilon_r}]
$ (L is linear operator on E)

$
L[\frac{dE}{d\varepsilon_r}] = \omega^2\mu_0\varepsilon_0E
$

$
L[\frac{dE}{d\sigma}] = -j\omega\mu_0E
$

Where adjoint comes in:

$
\frac{dF(\varepsilon_r, \mu_r, E)}{d\varepsilon_r}=
\frac{\partial F(\varepsilon_r, \mu_r, E)}{\partial\varepsilon_r} + 
\frac{\partial F(\varepsilon_r, \mu_r, E)}{\partial E}
L^{-1}[\omega^2\mu_0\varepsilon_0E]
$

$
\lambda = \frac{\partial F(\varepsilon_r, \mu_r, E)}{\partial E}
L^{-1}
$

$
L^H\lambda = \frac{\partial F(\varepsilon_r, \mu_r, E)}{\partial E}^H
$

$
\frac{dF(\varepsilon_r, \mu_r, E)}{d\varepsilon_r}=
\frac{\partial F(\varepsilon_r, \mu_r, E)}{\partial\varepsilon_r} + 
\lambda^H\omega^2\mu_0\varepsilon_0E^H
$

$
\frac{dF(\varepsilon_r, \mu_r, E)}{d\sigma}=
\frac{\partial F(\varepsilon_r, \mu_r, E)}{\partial\sigma} - 
\lambda^Hj\omega\mu_0E^H
$

## Integral Form

$
\int W\cdot(\nabla \times\frac{1}{\mu_r}\nabla \times E -\omega^2\varepsilon_r\mu_0\varepsilon_0E +j\omega\mu_0 \sigma E)dV = \int W\cdot j\omega\mu_0 J dV
$

$
\int W\cdot\nabla \times\frac{1}{\mu_r}\nabla \times E dV
+\int W\cdot (\omega^2\varepsilon_r\mu_0\varepsilon_0E +j\omega\mu_0 \sigma E)dV = \int W\cdot j\omega\mu_0 J dV
$

$
\int W\cdot\nabla \times\frac{1}{\mu_r}\nabla \times E dV
= \int \nabla \times W\cdot\frac{1}{\mu_r}\nabla \times E dV
-\oint (W \times \frac{1}{\mu_r}\nabla \times E) \cdot dS
$

Now this is just a number:

$
I[E,W;\varepsilon_r, \mu_r, \sigma]=\int \frac{1}{\mu_r}\nabla \times W\cdot\nabla \times E dV
-\oint (W \times \frac{1}{\mu_r}\nabla \times E) \cdot dS
-\int  (\omega^2\varepsilon_r\mu_0\varepsilon_0-j\omega\mu_0\sigma )W\cdot EdV - \int j\omega\mu_0  W\cdot J dV
$

Vector form ($\frac{\partial}{\partial W}$):

$
L[E;\varepsilon_r, \mu_r, \sigma]=\int \frac{1}{\mu_r}(\nabla \times N)\cdot(\nabla \times N)E dV
-\oint (Ni \times \frac{1}{\mu_r}\nabla \times N_j)E \cdot dS
-\int  (\omega^2\varepsilon_r\mu_0\varepsilon_0-j\omega\mu_0\sigma )(N_i\cdot N_j)EdV - \int j\omega\mu_0  N_i\cdot J dV
$

$
\frac{d}{d\varepsilon_r}(L[E; \varepsilon_r, \mu_r, \sigma]) = 
\frac{\partial L}{\partial \varepsilon_r} + 
\frac{\partial L}{\partial E}\frac{dE}{d\varepsilon_r}=0
$

### $\frac{\partial L}{\partial E}$=?

$
\frac{\partial L}{\partial E}=\int \frac{1}{\mu_r}(\nabla \times N)\cdot(\nabla \times N) dV
-\oint (Ni \times \frac{1}{\mu_r}\nabla \times N_j) \cdot dS
-\int  (\omega^2\varepsilon_r\mu_0\varepsilon_0-j\omega\mu_0\sigma )(N_i\cdot N_j)dV=K
$

### $\frac{\partial L}{\partial \varepsilon_r}$=?
$
\frac{\partial L}{\partial \varepsilon_r}=-\int  (\omega^2\mu_0\varepsilon_0)(N_i\cdot N_j)EdV
$
($\in \mathbb{R}^{N\times P}$)
(N size of E, P size of epsilon)

### $\frac{\partial L}{\partial \sigma}$=?
$
\frac{\partial L}{\partial \sigma}=\int  (j\omega\mu_0)(N_i\cdot N_j)EdV
$


### $\frac{\partial L}{\partial \mu_r}$=?
$
\frac{\partial L}{\partial \mu_r}=-\int \frac{1}{\mu_r^2}(\nabla \times N_i)\cdot(\nabla \times N_j)E dV
+\oint (Ni \times \frac{1}{\mu_r^2}\nabla \times N_j)E \cdot dS
$

Each of these derivatives of parameters is a vector of derivatives for each parameter varied. 
Efficiency comes from using only local parts of integrals that parameters affect, and adjoint of objective function. 

Parallel implementaiton would require a vecscatter of electric field and adjoint parameters.
N -> M*P

# Continuous Adjoint

$\nabla \times E = - \mu\frac{\partial H}{\partial t}$

$\nabla \times H =  \epsilon\frac{\partial E}{\partial t} + J$

Optimization problem:

$\min_p \int g(E,H,\epsilon, \mu)dt$

s.t.

$ \nabla \times E = - \mu\frac{\partial H}{\partial t}$

$\nabla \times H =  \epsilon\frac{\partial E}{\partial t} + J$

$$
L(E,H,\epsilon, \mu, \lambda_1, \lambda_2) = \int g(E,H,\epsilon, \mu) 
+ \lambda_1(\nabla \times E + \mu\frac{\partial H}{\partial t}) 
+ \lambda_2(\nabla \times H - \epsilon\frac{\partial E}{\partial t} - J) dt
$$

Stationary wrt $\lambda_1, \lambda_2$, solve ODE for E and H.

Stationary wrt $E,H$.

$$
0=\frac{\partial L}{\partial E}=  \int \frac{\partial g}{\partial E}dt
+ \int \frac{\partial}{\partial E}\lambda_1(\nabla \times E) dt
+ \int \frac{\partial}{\partial E}\lambda_2(- \epsilon\frac{\partial E}{\partial t}) dt
$$

$$
0=  \int \frac{\partial g}{\partial E}dt
+ \int -\nabla \times \lambda_1 dt
+ \int - \epsilon\frac{\partial \lambda_2}{\partial t} dt
$$

$$
0=\frac{\partial L}{\partial H}=  \int \frac{\partial g}{\partial H}dt
+ \int \frac{\partial}{\partial H}\lambda_1(\mu\frac{\partial H}{\partial t}) dt
+ \int \frac{\partial}{\partial H}\lambda_2(\nabla \times H) dt
$$

$$
0=\frac{\partial L}{\partial H}=  \int \frac{\partial g}{\partial H}
+ \mu\frac{\lambda_1}{\partial t}
-\nabla \times \lambda_2 dt
$$

Adjoint ODEs:

$$
\frac{\partial g}{\partial E}=\nabla \times \lambda_1 + \epsilon\frac{\partial \lambda_2}{\partial t}
$$

$$
-\frac{\partial g}{\partial H}=\mu\frac{\lambda_1}{\partial t}
-\nabla \times \lambda_2
$$




