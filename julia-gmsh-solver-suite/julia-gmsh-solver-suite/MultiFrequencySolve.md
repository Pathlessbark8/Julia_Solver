# Multi-frequency solve Math

$K(\omega)u(\omega)=b(\omega)$

## First Derivative
$\frac{d}{d\omega}[K\cdot u]=\frac{d}{d\omega}b$

$\frac{dK}{d\omega}u + K \frac{du}{d\omega}=\frac{db}{d\omega}$

$K \frac{du}{d\omega}=\frac{db}{d\omega} - \frac{dK}{d\omega}u$

## Second Derivative

$\frac{d}{d\omega}[\frac{dK}{d\omega}u + K \frac{du}{d\omega}]=\frac{d}{d\omega}[\frac{db}{d\omega}]$

$\frac{d^2K}{d\omega^2}u + 2\frac{dK}{d\omega}\frac{du}{d\omega} + K \frac{d^2u}{d\omega^2}=\frac{d^2b}{d\omega^2}$

$K \frac{d^2u}{d\omega^2}=\frac{d^2b}{d\omega^2} - \frac{d^2K}{d\omega^2}u - 2\frac{dK}{d\omega}\frac{du}{d\omega}$

## Derivative of K

$K_{ij} = \int \nabla\times N_i \cdot \nabla\times N_j - (\omega^2\varepsilon\mu-j\omega\mu\sigma ) N_i \cdot N_j dV$

$\frac{dK_{ij}}{d\omega} = \int -(2\omega\varepsilon\mu-j\mu\sigma ) N_i \cdot N_j dV$

$\frac{d^2K_{ij}}{d\omega^2} = \int - (2\varepsilon\mu) N_i \cdot N_j dV$

## Taylor Series of u

$u(\omega) = u(\omega_0) + \frac{du}{d\omega}d\omega + \frac{1}{2}\frac{d^2u}{d\omega^2}d\omega^2 + O(d\omega^3)$

$d\omega = \omega - \omega_0$
