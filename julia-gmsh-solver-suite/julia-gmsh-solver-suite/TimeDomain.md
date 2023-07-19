# Time Domain Equations

## Maxwell Stuff

Base equations:
$\nabla \times E = -\mu \frac{\partial H}{\partial t} - M$

$\nabla \times H = \varepsilon \frac{\partial E}{\partial t}  + \sigma E + J$

Scaling:

$\frac{1}{\mu_r}\nabla \times E = -\mu_0 \frac{\partial H}{\partial t} - \frac{1}{\mu_r}M$

$-\mu_0\nabla \times H = -\mu_0\varepsilon \frac{\partial E}{\partial t}  -\mu_0\sigma E -\mu_0J$

Take Curl of (1), time derivative of (2):

$\nabla \times\frac{1}{\mu_r}\nabla \times E = -\mu_0 \nabla \times\frac{\partial H}{\partial t} M$

$-\mu_0\nabla \times \frac{\partial H}{\partial t} = -\mu_0\varepsilon \frac{\partial^2 E}{\partial t^2}  -\mu_0\sigma \frac{\partial E}{\partial t} -\mu_0\frac{\partial J}{\partial t}$

Create curl curl:

$\nabla \times\frac{1}{\mu_r}\nabla \times E +\mu_0\varepsilon \frac{\partial^2 E}{\partial t^2}  +\mu_0\sigma \frac{\partial E}{\partial t} = -\mu_0\frac{\partial J}{\partial t} - \nabla \times\frac{1}{\mu_r}M$

Combine $\varepsilon_0\mu_0=\frac{1}{c_0^2}$

$\nabla \times\frac{1}{\mu_r}\nabla \times E + \frac{1}{c_0^2}\varepsilon_r \frac{\partial^2 E}{\partial t^2}  +\frac{\eta_0}{c_0}\sigma \frac{\partial E}{\partial t} = -\mu_0\frac{\partial J}{\partial t} - \nabla \times\frac{1}{\mu_r}M$

With $\frac{\eta_0}{c_0} = \frac{\sqrt{\mu_0/\varepsilon_0}}{c_0}= \sqrt{\mu_0/\varepsilon_0}\sqrt{\varepsilon_0\mu_0} = \mu_0$

## Time Stepping

# Discrete Equations

$\frac{1}{c_0^2}T_\varepsilon \frac{\partial^2 e}{\partial t^2}  +\frac{\eta_0}{c_0}T_\sigma \frac{\partial e}{\partial t} + S_\mu e = -f(t)$

$T_\varepsilon \frac{e^{n+1} - 2e^n + e^{n-1}}{\Delta t^2c_0^2}  + \eta_0T_\sigma \frac{e^{n+1} - e^{n-1}}{2\Delta tc_0}  + S_\mu (\beta e^{n+1} + (1-2\beta) e^{n} + \beta e^{n-1}) = -(\beta f^{n+1} + (1-2\beta) f^{n} + \beta f^{n-1})$

Let $\tau = \Delta t c_0$

$T_\varepsilon \frac{e^{n+1} - 2e^n + e^{n-1}}{\tau^2}  + \eta_0T_\sigma \frac{e^{n+1} - e^{n-1}}{2\tau}  + S_\mu (\beta e^{n+1} + (1-2\beta) e^{n} + \beta e^{n-1}) = -(\beta f^{n+1} + (1-2\beta) f^{n} + \beta f^{n-1})$

Find $e^{n+1}$:

$[\frac{1}{\tau^2}T_\varepsilon + \frac{\eta_0}{2\tau}T_\sigma + \beta S]e^{n+1}=[\frac{2}{\tau^2}T_\varepsilon - (1-2\beta)S]e^n - [\frac{1}{\tau^2}T_\varepsilon - \frac{\eta_0}{2\tau}T_\sigma + \beta S]e^{n-1} -(\beta f^{n+1} + (1-2\beta) f^{n} + \beta f^{n-1})$



