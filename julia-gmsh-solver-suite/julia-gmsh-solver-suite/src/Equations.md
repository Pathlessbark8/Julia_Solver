# Derivation of the Helmholtz Equation for 2D TM Problems

Starting from:

$j\omega \varepsilon \vec{E} - \nabla \times \vec{H} = -\vec{J}$, $\qquad \qquad$ $j\omega \mu \vec{H} + \nabla \times \vec{E} = \vec{0}$

The vector wave equations is obtained by solving for $\nabla \times \vec{H}$:

$\nabla \times \vec{H} = -\frac{1}{j \omega \mu}\nabla \times \nabla \times \vec{E}$

Giving:

$j\omega \varepsilon \vec{E} + \frac{1}{j \omega \mu}\nabla \times \nabla \times \vec{E} = -\vec{J}$

Multiplying through by $j\omega \mu$:

$j^2\omega^2 \varepsilon \mu \vec{E} + \nabla \times \nabla \times \vec{E} = -j \omega \mu \vec{J}$

or 

$-k^2\vec{E} +  \nabla \times \nabla \times \vec{E} = -j \omega \mu \vec{J}$

For 2D TM problems with no variation in the $z$ direction, you can work out the fact that the z-component of the curl-curl operator is:

$\left(\nabla \times \nabla \times \vec{E}\right)_z = -\partial_{xx} E_z - \partial_{yy}E_z = -\nabla^2 E_z$ (note the negative sign).

Therefore the scalar Helmholtz equation of interest is:

$-k^2E_z -\nabla^2 E_z = -j \omega \mu J_z$

or:

$k^2E_z +\nabla^2 E_z = j \omega \mu J_z$

# Derivation of signs for Galerkin testing of scalar Helmholtz equation

If we test the scalar Helmholtz equation with $\psi$ then:

$\int \psi k^2E_z d\vec{x} + \int \psi \nabla^2 E_z d\vec{x} = j \omega \mu \int \psi J_z d \vec{x}$

The vector identity $\nabla \cdot \psi \nabla E_z = \psi \nabla^2 E_z + \nabla \psi \cdot \nabla E_z$ is useful here (effectively the product rule). We can re-write as:

$\psi \nabla^2 E_z = \nabla \cdot \psi \nabla E_z - \nabla \psi \cdot \nabla E_z$

which when substituted into the tested equations results in:

$\int \psi k^2E_z d\vec{x} - \int \nabla \psi \cdot \nabla E_z d\vec{x} + \int \nabla \cdot \left(\psi \nabla E_z\right) d\vec{x} = j \omega \mu \int \psi J_z d \vec{x}$

and the divergence theorem results in

$\int \psi k^2E_z d\vec{x} - \int \nabla \psi \cdot \nabla E_z d\vec{x} + \oint \hat{n} \cdot \left(\psi \nabla E_z\right) d\vec{x} = j \omega \mu \int \psi J_z d \vec{x}$

From the perspective of implementing the volumetric terms, (the surface term can be ignored if FEM basis functions enforcing continuity are enforced) we have:

$ k^2 \cal{M} E_z - \cal{S}E_z = \cal{M}J_z$

where the stiffness matrix $\cal{S} = \cal{S}_{xx} + \cal{S}_{yy}$ with entries equal to $\left(\cal{S}_{xx}\right)_{ij} = \int (\partial_x b_i) (\partial_x b_j) d\vec{x}$. Here we assume that $k$ is a constant on a given element (does not vary locally).


# Derivation of the Helmholtz Equation for 2D TE Problems

The vector wave equation we start with is:

$-k^2\vec{E} +  \nabla \times \nabla \times \vec{E} = -j \omega \mu \vec{J}$

When we use a curl-confirming vector basis function expansion, we will expand $\vec{E}$ as a sum of scalar coefficients times vector basis functions $\vec{\psi}$. Testing will then involve integrating the dot product: 

$\vec{\psi} \cdot \nabla \times (\vec{F}) = \nabla \cdot (\vec{F} \times \vec{\psi}) + (\nabla \times \vec{\psi}) \cdot \vec{F}$

such that

$\vec{\psi} \cdot \nabla \times (\nabla \times \vec{E}) = \nabla \times ((\nabla \times \vec{E}) \times \vec{\psi} ) + (\nabla \times \vec{\psi}) \cdot (\nabla \times \vec{E})$

It follows that the tested vector wave equation can be written as

$-k^2 \int \vec{\psi} \cdot \vec{E} d\vec{x} + \int (\nabla \times \vec{\psi})\cdot(\nabla \times \vec{E}) d\vec{x} + \oint \hat{n} \times ((\nabla \times \vec{E}) \times \vec{\psi}) = -j\omega\mu \int \vec{\psi}\cdot \vec{J} d\vec{x}$

As the curl is given by:

$\nabla \times \vec{E} = \left|\begin{matrix} \hat{x} & \hat{y} & \hat{z} \\ \partial_x & \partial_y & \partial_z \\ E_x & E_y & E_z \end{matrix}\right| = \begin{bmatrix} \partial_y E_z  - \partial_z E_y \\ \partial_z E_x - \partial_x E_z \\ \partial_x E_y - \partial_y E_x \end{bmatrix}$

The TE problem with no variation in $z$ results in the curl:

$\nabla \times \vec{E} = \begin{bmatrix} \partial_y E_z \\ - \partial_x E_z \\ \partial_x E_y - \partial_y E_x \end{bmatrix}$

As the field equations uncouple so that we can solve for $E_x$ and $E_y$ independently of $E_z$, we take the curl:

$\nabla \times \vec{E} = \begin{bmatrix} 0 \\ 0 \\ \partial_x E_y - \partial_y E_x \end{bmatrix}$

Such that 

$(\nabla \times \vec{\psi})\cdot(\nabla \times \vec{E}) = \begin{bmatrix} 0 \\ 0 \\ \partial_x \psi_y - \partial_y \psi_x \end{bmatrix}^{T} \begin{bmatrix} 0 \\ 0 \\ \partial_x E_y - \partial_y E_x \end{bmatrix} = (\partial_x \psi_y)(\partial_x E_y) + (\partial_y \psi_x)(\partial_y E_x) - (\partial_x \psi_y) (\partial_y E_x) - (\partial_y \psi_x)(\partial_x E_y) $

However, Gmsh provides a way to evaluate both the vector basis function $\vec{\psi}$ and its curl directly. Which implies we only need to integrate the product of the $z$-components of the basis functions to make progress.