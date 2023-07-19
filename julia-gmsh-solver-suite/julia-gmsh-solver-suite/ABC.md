# Vector ABC Condition

## First Order condition:

$\hat{n} \times \nabla \times E = -jk \hat{n}\times\hat{n}\times E$

FEM Left hand side integral:

$\int \psi \cdot \nabla \times \nabla \times E - k^2 \psi \cdot E dV$

$=\int \nabla \times \psi \cdot  \nabla \times E - k^2 \psi \cdot E dV - \int (\psi \times \nabla \times E) \cdot \hat{n}dS$

$=\int \nabla \times \psi \cdot  \nabla \times E - k^2 \psi \cdot E dV - \int \psi \cdot ( \hat{n} \times \nabla \times E) dS$

Substituting in ABC condition:

$=\int \nabla \times \psi \cdot  \nabla \times E - k^2 \psi \cdot E dV - \int \psi \cdot (-jk \hat{n} \times \hat{n} \times E) dS$

$=\int \nabla \times \psi \cdot  \nabla \times E - k^2 \psi \cdot E dV - \int jk (\hat{n} \times \psi) \cdot ( \hat{n} \times E) dS$

## Second Order Condition

$\hat{r} \times \nabla \times E = -jk \hat{r}\times\hat{r}\times E + \frac{r}{2(jkr+1)} \times \{\nabla \times [\hat{r}(\nabla\times E)_r] + (s-1)\nabla_t (\nabla\cdot E_t) + (2-s)jk\nabla_t E_r\}$

Jianming-Jin says set s=2 to preserve symmetry (s=1/2 optimal):

$\hat{r} \times \nabla \times E = -jk \hat{r}\times\hat{r}\times E + \frac{r}{2(jkr+1)} \times \{\nabla \times [\hat{r}(\nabla\times E)_r] + \nabla_t (\nabla\cdot E_t) \}$

If we have zero divergence?

$\hat{r} \times \nabla \times E = -jk \hat{r}\times\hat{r}\times E + \frac{r}{2(jkr+1)} \times \{\nabla \times [\hat{r}(\hat{r}\cdot(\nabla\times E))] \}$

Let $\hat{n}=\hat{r}$.

$I_{ABC} = - \int (\psi \times \nabla \times E) \cdot \hat{r}dS$

$I_{ABC} = - \int \psi \cdot ( \hat{r} \times \nabla \times E) dS$

$I_{ABC} = - \int \psi \cdot ( -jk \hat{r}\times\hat{r}\times E + \frac{r}{2(jkr+1)} \times \{\nabla \times [\hat{r}(\hat{r}\cdot(\nabla\times E))] \}) dS$

$I_{ABC} = - \int \psi \cdot ( -jk \hat{r}\times\hat{r}\times E)dS - \int \psi \cdot (\frac{r}{2(jkr+1)} \times \{\nabla \times [\hat{r}(\hat{r}\cdot(\nabla\times E))] \})) dS$

$I_{ABC} = - \int jk (\hat{n} \times \psi) \cdot ( \hat{n} \times E)dS - \int \psi \cdot (\frac{r}{2(jkr+1)} \times \{\nabla \times [\hat{r}(\hat{r}\cdot(\nabla\times E))] \}) dS$

$I_{ABC} = - \int jk (\hat{n} \times \psi) \cdot ( \hat{n} \times E)dS - \int (\frac{r}{2(jkr+1)})(\hat{r}\cdot \nabla \times \psi)(\hat{r}\cdot \nabla \times E) dS$


$K = \int \nabla \times N_i \cdot \nabla \times N_j - k^2 N_i \cdot N_j dV$

$K_{abc} = \int jk (\hat n \times N_i) \cdot (\hat n \times N_j) + \beta (\hat n \cdot \nabla \times N_i)(\hat n \cdot \nabla \times N_j) dV$

$M = \int k^2 N_i \cdot N_j dV$