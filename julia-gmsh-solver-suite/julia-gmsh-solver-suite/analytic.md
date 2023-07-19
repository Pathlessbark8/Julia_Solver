# TE Field from Mz Point Source

$F(\vec r) = \frac{\epsilon}{4\pi}\int M(\vec r')\frac{e^{-jkR}}{R}dV'$

$R=|\vec r - \vec r'|$

$M(r') = \hat z\delta(\vec r')$

$F_z = \frac{\epsilon}{4\pi}\frac{e^{-jkR}}{R}$

$\vec E = -\frac{1}{\epsilon} \nabla \times \vec F$

$\vec E = -\frac{1}{4 \pi} \nabla \times \hat z f(x,y)$

$f(x,y) = \frac{e^{-jkR}}{R}$

$\vec E = -\frac{1}{4 \pi} [\hat x \frac{\partial f(x,y)}{\partial y} - \hat y \frac{\partial f(x,y)}{\partial x}]$

$\vec E = -\frac{1}{4 \pi} [\hat x \frac{\partial}{\partial y}\frac{e^{-jkR}}{R} - \hat y \frac{\partial}{\partial x}\frac{e^{-jkR}}{R}]$

$\frac{\partial}{\partial R}\frac{e^{-jkR}}{R} = \frac{e^{-jkR}}{R}\frac{jkR+1}{R}$

$\frac{\partial R}{\partial x} = \frac{x-x'}{R}$

$\frac{\partial R}{\partial y} = \frac{y-y'}{R}$

$\vec E = -\frac{1}{4 \pi} [\hat x \frac{e^{-jkR}}{R}\frac{jkR+1}{R}\frac{y-y'}{R} - \hat y \frac{e^{-jkR}}{R}\frac{jkR+1}{R}\frac{x-x'}{R}]$

$\vec E = -\frac{1}{4 \pi} \frac{e^{-jkR}(jkR+1)}{R^3}[\hat x (y-y') - \hat y (x-x')]$

