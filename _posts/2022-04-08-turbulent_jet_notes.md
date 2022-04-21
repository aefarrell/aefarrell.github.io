---
title: "Turbulent Jets"
toc: true
toc_label: "contents"
toc_sticky: true
comments: true
categories:
  - notes
tags:
  - turbulent jets
  - chemical releases
tagline: "Notes on turbulent jets and velocity profiles."
header:
  overlay_image: /images/chuttersnap-unsplash-header.jpg
  caption: "Photo by [CHUTTERSNAP](https://unsplash.com/@chuttersnap?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText) on [Unsplash](https://unsplash.com/)"
---

# Turbulent Jets

In [a previous post](https://aefarrell.github.io/2021/04/10/turbulent_jet_example/) I worked through a chemical release modeled as a turbulent jet and while I mentioned there were several ways modeling the jet, I didn't go into any of them. I'm taking the opportunity here to collect my notes on turbulent jets, some different ways of modeling the jets, and the relative performance of each approach.

## Observations on Turbulent Jets

We are considering a submerged circular jet, issuing from a surface, with the coordinate system centered on the jet. Since it is circular, the natural coordinate system is cylindrical with a downstream distance *z*, radial distance *r*, and angular coordinate *&theta;*. The jet is fully turbulent when the Reynolds number, $Re \gt 2000$, where the Reynolds number is calculated with respect to the initial jet velocity and jet diameter

$$ Re = { \rho_j v_0 d_0 \over \mu_j } $$

We are also considering the case where the densities of the two fluids are similar, where we take "similar" to mean 
$$ \frac{1}{4} \le { \rho_{a} \over \rho_{j} } \le 4 $$

Where subscript *a* indicates the ambient fluid and *j* the jet. For much the experimental data the jet and ambient fluid are the same fluid, e.g. a jet of air into air or water into water.

Turbulent jets expand by entraining ambient fluid, tracing out a cone defined by a jet angle $\alpha \approx 15-25^\circ$. The mixing layer penetrates into the jet forming the potential core, inside is pure jet material and outside is mixed. After approximately 6 diameters the region is fully developed.

![image.png](/images/turbulent_jet_notes_files/att1.png)

Empirical approximations of the velocity profile are often given with respect to this jet angle or, equivalently, the slope of the line (i.e. $\tan \frac{\alpha}{2}$). A related way of parameterizing the jet is in terms of a width parameter *b*. Typically this is the width of the velocity profile at half-height $b_{1/2}$ (though not always). With a constant jet angle and a self-similar velocity profile the width is directly proportional to the downstream distance $ b_{1/2} = \tan \left( \frac{\alpha_{1/2} }{2} \right) z = c z$.

Where the value of *c* is can be found in the literature

|   c           | Reference                     |
|:-------------:|:------------------------------|
| 0.082 - 0.097 | Garde, 2010[^garde]           |
| 0.0848        | Bird, 2007[^bird]             |
| 0.10          | Rajaratnam, 1974[^rajaratnam] |


[^garde]: R. J. Garde, *Turbulent Flows, 3rd Ed.*, New AcademicScience, London (2010)

[^bird]: B. R. Bird, W. E. Stewart, E. N. Lightfoot, *Transport Phenomena, 2nd Ed.*, John Wiley & Sons, New York (2007)

[^rajaratnam]: N. Rajaratnam, *Turbulent Jets*, Elsevier, Amsterdam (1974)

At this point it is common to introduce a variable $ \xi = {r \over b_{1/2} }$ or $ \xi = {r \over z }$ where we are taking advantage of the fact that $ b_{1/2} \propto z $. This is a *scaled* radial distance, using the width at half-height as a characteristic length. It is important to keep track of which definition of &xi; is being used as they differ by a scaling factor. The reason for this change of variables is the observation that the *shape* of the velocity profile is the same at any downstream point, it is merely scaled down in height and wider as one travels downstream. That is $ { \bar{v}_z \over \bar{v}_{max} } = f \left( \xi \right) $ is the same for all downstream distances (in the region where the jet is fully developed).

Another important observation is that the center-line velocity, the max velocity in the jet, decays with the inverse of the downstream distance, i.e.

$$ \bar{v}_{max} \propto z^{-1} $$

Putting those two observations together we expect the velocity profile to have the form

$$ \bar{v}_z = { \mathrm{const} \over z } f \left( \xi \right)$$

## Modeling Turbulent Jets

To set up our system we consider the case of a jet coming out of a point on an infinite surface into a quiescent medium, and that the jet and medium have the same density. This is a major simplification, but it makes the math easier to deal with. The coordinate system is centered at this point and all momentum in the jet ultimately comes from the origin.

The boundary conditions for the problem are:
1. at the center-line, *r=0*, the velocity is entirely in the z-directly
1. at the center-line, *r=0*, the velocity in the z-direction is at a maximum
1. as the radius increases, *r &rarr; &infin;* , the velocity in the z-direction goes to zero

### Time Averaged Values

Since we are concerned with turbulent flow, we can employ [Reynolds decomposition](https://en.wikipedia.org/wiki/Reynolds_decomposition) to transform the velocities like so

$$ v_z = \bar{v}_z + v^{\prime}_{z} $$
$$ v_r = \bar{v}_r + v^{\prime}_{r} $$

where $\bar{v}$ is the time-smoothed velocity and $v^{\prime}$ is an instantaneous deviation such that $\bar{v^{\prime} } = 0$ and the time-averaging operator follows the [Reynolds criteria](https://en.wikipedia.org/wiki/Reynolds_operator).



### Equations of Motion

The equations of motion in terms of time-smoothed velocities are
$$ \rho {D \mathbf{\bar{v} } \over D t } = - \nabla \bar{p} - \nabla \cdot \mathbf{ \bar{\tau} } + \rho \mathbf{g}  $$

Where $\mathbf{ \bar{\tau} }$ is the turbulent stress and includes the [Reynolds stresses](https://en.wikipedia.org/wiki/Reynolds_stress).

With the z component, in cylindrical coordinates[^bird]

$$ \rho \left( {\partial \over \partial t} \bar{v}_z + \bar{v}_r {\partial \bar{v}_z \over \partial r} + {\bar{v}_\theta \over r} {\partial \bar{v}_z \over \partial \theta} + \bar{v}_z {\partial \bar{v}_z \over \partial z} \right) 
\\ = - {\partial \bar{p} \over \partial z} - {1 \over r} {\partial \left( r \bar{\tau}_{rz} \right) \over \partial r } - {1 \over r} {\partial \bar{\tau}_{\theta z} \over \partial \theta } - {\partial \bar{\tau}_{z z} \over \partial z } + \rho g_z$$

Making the assumptions:
1. Zero pressure gradient ( ${\partial p \over \partial z} = 0$ )
1. Steady state ( ${\partial \over \partial t} \left( \cdots \right) = 0$ )
1. Axisymmetric ( ${\partial \over \partial \theta} \left( \cdots \right) = 0$ )
1. Effect of gravity can be neglected ( $ \rho g_z \approx 0$ )
1. Within the jet $ \mid v_z \mid \gg \mid v_r \mid $ and, by boundarly layer approximation, $ \bar{\tau}_{z z} $ can be neglected[^boundarylayer]

The equations of motion, in the z direction, simplifies to

$$ \bar{v}_r {\partial \bar{v}_z \over \partial r} + \bar{v}_z {\partial \bar{v}_z \over \partial z} = - {1 \over \rho r} {\partial \left( r \bar{\tau}_{rz} \right) \over \partial r } $$


[^boundarylayer]: The boundary layer approximation is that 
    $$ { \partial^2 \bar{v}_z \over \partial z^2 } \ll { \partial^2 \bar{v}_z \over \partial r^2 } $$
    and if we suppose that
    $$ \bar{\tau}_{z z} \propto { \partial \bar{v}_z \over \partial z } $$
    and 
    $$ \bar{\tau}_{r z} \propto {\partial \bar{v}_z \over \partial r} $$
    we find 
    $$ { \partial \bar{\tau}_{z z} \over \partial z } \propto {\partial^2 \bar{v}_z \over \partial z^2} \ll { \partial r \bar{\tau}_{r z} \over \partial r} \propto {\partial^2 \bar{v}_z \over \partial r^2} $$
    and thus we can assume the free turbulence is dominated by 
    $$ \bar{\tau}_{r z} $$
    and
    $$ { \partial \bar{\tau}_{z z} \over \partial z} \approx 0 $$

### Equation of Continuity

The continuity equation in terms of time-smoothed velocities is

$$ {\partial \rho \over \partial t} + \nabla \cdot \rho \mathbf{ \bar{v} } = 0 $$

In cylindrical coordinates[^bird]

$$ {\partial \rho \over \partial t} + {1 \over r} {\partial \rho r \bar{v}_r \over \partial r} + {1 \over r} { \partial \rho \bar{v}_\theta \over \partial \theta} + {\partial \rho \bar{v}_z \over \partial z} = 0 $$

Making the assumptions:
1. Steady state ( ${\partial \over \partial t} \left( \cdots \right) = 0$ )
1. Axisymmetric ( ${\partial \over \partial \theta} \left( \cdots \right) = 0$ )
1. Incompressible ( ${\partial \rho \over \partial z} = {\partial \rho \over \partial r} = {\partial \rho \over \partial \theta} = 0 $ )

The equation of continuity simplifies to

$$ {1 \over r} {\partial r \bar{v}_r \over \partial r} + {\partial \bar{v}_z \over \partial z} = 0 $$

### Stokes Stream Function

To simplify things down to working with one dependent variable we introduce a [Stokes stream function](https://en.wikipedia.org/wiki/Stokes_stream_functionhttps://en.wikipedia.org/wiki/Stokes_stream_function) $\psi$ defined such that

$$ \bar{v}_z = -{1 \over r} {\partial \psi \over \partial r} $$

and 

$$ \bar{v}_r = {1 \over r} {\partial \psi \over \partial z} $$

This definition ensures that the equation of continuity is satisfied. Suppose that $ \psi = k z F\left(\xi\right) $, where *F* is a unitless function of $\xi = \frac{r}{z}$ and $k$ is a constant with units $ [[ \mathrm{length} ]]^2 \times [[ \mathrm{time} ]]^{-1} $, then

$$ \bar{v}_z = -{1 \over r} {\partial \xi \over \partial r} {\partial \psi \over \partial \xi}
= -{1 \over r} {1 \over z} {k z F^{\prime} } \\
\\= -{k \over z} {F^{\prime} \over \xi} = { \mathrm{const} \over z } f \left( \xi \right)$$

Which matches what we expect from the empirical observations (which is why we supposed that form of the stream function in the first place). We can use this definition to work out some other useful terms

$$ {\partial \bar{v}_z \over \partial z} = {k \over z^2} F^{\prime \prime} $$

$$ {\partial \bar{v}_z \over \partial r} = -{k \over z^2} \left( { F^{\prime \prime} \over \xi} - { F^{\prime} \over \xi^2} \right) $$

$$ \bar{v}_r = { k \over z } \left( { F \over \xi } - F^{\prime} \right) $$

Substituting these back into the equation of motion, in the z direction, leads to

$$ \left( k \over z \right)^2 \left[ { F F^{\prime \prime} \over \xi } - {F F^{\prime} \over \xi^2} + { \left( F^{\prime} \right)^2 \over \xi } \right]  = {1 \over \rho} {\partial \over \partial r} \left( r \bar{\tau}_{rz} \right) $$

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = {1 \over \rho} {\partial \over \partial r} \left( r \bar{\tau}_{rz} \right) $$

Which is suggestive of the overall approach to follow: find an expression for the right hand side of this differential equation, integrate both sides with respect to *&xi;*, and solve for *F(&xi;)*

### Boundary Conditions

The initial boundary conditions of the problem were that:
1. $\bar{v}_r = 0$ at *r=0*
1. ${\partial \bar{v}_z \over \partial r} = 0$ at *r=0* (i.e. the velocity is at a maximum)
1. $\bar{v}_z \to 0$ as *r &rarr; &infin;* (i.e. the velocity decays to zero)

In terms of *F* and *&xi;* these become:
1. ${F \over \xi} - F^{\prime} = 0$ at *&xi;=0*, which implies *F=0* at *&xi;=0*
1. $F^{\prime \prime} - {F^{\prime} \over \xi}  = 0$ at *&xi;=0*
1. ${F^{\prime} \over \xi} \to 0$ as *&xi; &rarr; &infin;*

### Momentum Balance

To determine the constant *k* we use a momentum balance: the momentum flux, *J*, in the z direction is constant. Initially the momentum flux is

$$ J = \rho v_0^2 A_0 = \rho v_0^2 {\pi \over 4} d_0^2$$

and at some point *z* downstream of the origin we have

$$ J = \int_{0}^{2\pi} \int_{0}^{\infty} \rho \bar{v}_z^2 r dr d\theta \\
     = 2 \pi \rho \int_{0}^{\infty} \bar{v}_{z,max}^2 \left( \bar{v}_z \over \bar{v}_{z,max} \right)^2 r dr \\
     = 2 \pi \rho \bar{v}_{z,max}^2 \int_{0}^{\infty} \left( \bar{v}_z \over \bar{v}_{z,max} \right)^2 r dr \\
     = 2 \pi \rho k^2 \int_{0}^{\infty} \left( f\left( \xi \right) \right)^2 \xi d \xi $$
     
Taking the integral to be *I*, and equating the initial momentum flux with the momentum flux at point *z*

$$ J = \rho v_0^2 {\pi \over 4} d_0^2 = 2 \pi \rho k^2 I $$

$$ k = \sqrt{1 \over 8 I } v_0 d_0 $$

**Note**  I've played a little fast and loose with the definition of $\bar{v}_z$ in that I am implicitly assuming $f(\xi) = {-F^{\prime}(\xi) \over \xi}$ which isn't strictly true, there can be scaling factor. In practice all of these are collected together into one constant so it doesn't matter, but that is something to be aware of as the definition of *k* here is really $k\times \mathrm{const}$ where $\mathrm{const} = {-F^{\prime}(\xi) \over \xi} \div f(\xi) $
{: .notice}

## Prandtl Mixing Length

The Prandtl mixing length model makes the assumption that momentum transfer occurs over some "mixing length" *l* such that

$$ \bar{\tau}_{rz} = -\rho l^2 \left| {\partial \bar{v}_z \over \partial r} \right| \left( {\partial \bar{v}_z \over \partial r} \right)$$

We suppose that the mixing length is proportional to the width of the velocity profile $b_{1/2}$, the characteristic length for the velocity profile, which we know is proportional to the downstream distance *z*

$$ l \propto b_{1/2} \propto z$$

$$ l = c z $$

Where $c$ is some unitless constant. Making the observation that $ {\partial \bar{v}_z \over \partial r} < 0 $ we can make the simplification

$$ \bar{\tau}_{rz} = \rho c^2 z^2 \left( {\partial \bar{v}_z \over \partial r} \right)^2$$

### Setting up the ODE

Recall that the equation of motion in the z direction is (in terms of the unitless function *F*) is

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = {1 \over \rho} {\partial \over \partial r} \left( r \bar{\tau}_{rz} \right) $$

Substituting the expression for $ \bar{\tau}_{rz} $ we have

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = c^2 z^2 {\partial \over \partial r} \left( r \left( \partial \bar{v}_{z} \over \partial r \right)^2 \right) $$

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = c^2 z^2 \left( \left( \partial \bar{v}_{z} \over \partial r \right)^2 + 2r \left( \partial \bar{v}_{z} \over \partial r \right) \left( \partial^2 \bar{v}_{z} \over \partial r^2 \right) \right) $$

Substituting in the expressions for $ {\partial \bar{v}_z \over \partial r} $ and $ {\partial^2 \bar{v}_z \over \partial r^2} $ we arrive at

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = c^2 \left(k \over z \right)^2 \left( 1 \over \xi \right)\left( F^{\prime \prime} - { F^{\prime} \over \xi } \right) \left(  2 F^{\prime \prime \prime} - 3 { F^{\prime \prime} \over \xi } + { F^{\prime} \over \xi^2 }\right) $$

$$ { d \over d \xi } \left( F F^{\prime} \over \xi \right) = c^2 { d \over d \xi } \left( 1 \over \xi \right)\left( F^{\prime \prime} - { F^{\prime} \over \xi } \right)^2 $$

Integrating both sides

$$ \left( F F^{\prime} \over \xi \right) = c^2 \left( 1 \over \xi \right)\left( F^{\prime \prime} - { F^{\prime} \over \xi } \right)^2 + \mathrm{const}$$

By applying the boundary conditions we find the constant of integration is zero, thus

$$ F F^{\prime} = c^2 \left( F^{\prime \prime} - { F^{\prime} \over \xi } \right)^2 $$

Making the substitution $ \phi = a^{-1} \xi $ where $ a = c^{2/3} $^[consts]

$$ F F^{\prime} = \left( F^{\prime \prime} - { F^{\prime} \over \phi } \right)^2 $$

$$ F^{\prime \prime} = { F^{\prime} \over \phi } + \sqrt{ F F^{\prime} } $$

Which is in a form that can be solved numerically.

[^consts]: parameterizing this in terms of *a* is the standard way of presenting the problem

### Solving the ODE

We can solve the ODE and perform the integral needed for the momentum balance at the same time. First we define a vector *u* such that:

$$ \mathbf{u} =  \begin{bmatrix} u_{1} \\ u_{2} \end{bmatrix} = \begin{bmatrix} F \\ F^{\prime} \end{bmatrix}$$

The ODE then becomes:

$$ {d \mathbf{u} \over dt } = \begin{bmatrix} F^{\prime} \\ F^{\prime \prime} \end{bmatrix}  = \begin{bmatrix} u_{2} \\ \frac{ u_{2} }{t} + \sqrt{ u_{1} u_{2} } \end{bmatrix} $$

Which has a singularity at *t=0*, but one that can be easily dealt with by setting the initial value of the derivatives to[^init]

$$ {d \mathbf{u} \over dt }_{t=0} = \begin{bmatrix} 0 \\ -1 \end{bmatrix}$$

Putting that together, the ODE can be integrated easily.

**Note**  Because of how $\bar{v}_z$ and $\bar{v}_r$ were defined $ -{ F^{\prime} \over \phi } \ge 0 $, i.e. $ { F^{\prime} \over \phi } \le 0 $. For the signs to work out, $ F \le 0 $ and $ F^{\prime} \le 0 $ (since $F F^{\prime} \ge 0$).
{: .notice}


[^init]: From the boundary conditions we know *F'(0) = 0* but what about *F''*? Taking the ratio 
    $$ { \bar{v}_z \over \bar{v}_{z,max} }_{r=0} = - {F^{\prime} \over \phi }_{\phi=0} = 1 $$
    we find ${F^{\prime} \over \phi } = -1$ at *&phi; = 0* and, from the boundary conditions, 
    $$ F^{\prime \prime} = {F^{\prime} \over \phi } $$
    at *&phi; = 0*, therefore *F''(0) = -1*


```julia
using StaticArrays
using DifferentialEquations: ODEProblem, Tsit5, solve, TerminateSteadyState

function sys(u,p,t)
    u₁, u₂ = u[1], u[2]
    if t > 0.0
        du₁ = u₂
        du₂ = u₂/t + √(max((u₁*u₂),0))
    else
        du₁ = 0.0
        du₂ = -1.0
    end
    
    return SA[du₁; du₂]
end

u0    = SA[0.0; 0.0]
tspan = (0.0, 6.0)
prob  = ODEProblem(sys, u0, tspan)
sol   = solve(prob, Tsit5(), dtmax=0.1, cb=TerminateSteadyState())

print(sol.retcode)
```

    Success


![svg](/images/turbulent_jet_notes_files/output_13_0.svg)


Using the solution in terms of &phi; we can write a function *f(&xi;)*


```julia
function f_pml(ξ; a=0.066)
    ϕ = abs(ξ)/a
    
    if ϕ >0
        F, F′ = sol(ϕ)
        f = -F′/ϕ
        f  = max(f, 0)
    else
        f = 1
    end
    
    return f
end
```


### Comparison with Tollmien

The classic treatment of the Prandtl mixing length model is from Tollmien[^tollmien] in which, instead of solving numerically in the way shown above, the ODE is further transformed and a series expansion which used to generate a table of results. More often than not it is these tabulated values, or similar ones[^rajaratnam2], that are presented as the solution to the model.

We can easily compare the result here with the tabulated values and verify for ourselves that we have indeed solved the right differential equation. Though by solving numerically in this way we can control the level of precision and easily generate smooth interpolations. In my opinion, this makes using the ODE solution far more convenient than the tabulated values.

[^tollmien]: W. Tollmien, *Zeitschrift f&uuml;r angewandte Mathematik und Mechanik*, **6**, 468-478 (1926), reprinted and translated in [NACA-TM-1085](https://ntrs.nasa.gov/search?reportNumber=NACA-TM-1085)

[^rajaratnam2]: see (Rajaratnam 1974)[^rajaratnam] page 39 which has a table with &phi; spaced every 0.1, though it has an error at &phi;=1 the value of ${F^{\prime} \over \phi }$ should be 0.606 but is given as 0.505 (presumably a typo)


![svg](/images/turbulent_jet_notes_files/output_19_0.svg)


### Width at Half Height

The width at half height, $b_{1/2}$, is an important parameter and often velocity profiles are scaled relative to this. To compare different models on a fair basis, it is a good idea to determine what the model parameters are relative to $b_{1/2}$. Then each model can be scaled to the same $b_{1/2}$ and compared, apples-to-apples.

In this case we don't have a closed form for the velocity profile so we need to solve for &phi; such that *f(&phi;)=0* numerically.


```julia
using Roots: find_zero

ϕ_half = find_zero( ϕ -> f_pml(ϕ; a=1)-0.5, (1, 1.25))
```


    1.2277667062444657


and we then write the model parameter *a* in terms of $b_{1/2}$

$$ a = \frac{1}{ \phi_{1/2} } \frac{ b_{1/2} }{z} = 0.814 \frac{ b_{1/2} }{z} $$

Using a default value for $ { b_{1/2} \over z } = 0.0848 $ we arrive at


```julia
b_half = 0.0848

a = b_half/ϕ_half
```


    0.06906849613098656


Several sources have tabulated values for *a*

| a     | Reference                     |
|:-----:|:------------------------------|
| 0.063 | Tollmien, 1926[^tollmien]     |
| 0.066 | Rajaratnam, 1974[^rajaratnam] |

and the result of this notebook compares with those

### Velocity Profile

Now that we have completed the integration we can calculate the parameter *k*, using the equation derived from the momentum balance

$$ k = \sqrt{1 \over 8 I } v_0 d_0 $$

with the value of the integral coming directly from the ode solver


```julia
using NumericalIntegration: integrate

ϕ, F′ = sol.t, sol[2,:]

# trim any unphysical values
F′[F′.>0] .= 0.0

function integrand(ϕ, F′)
    if ϕ>0
        return F′^2/ϕ
    else
        return 0
    end
end

I = integrate(ϕ, integrand.(ϕ, F′))
I = a^2 * I
```


    0.002573069044757039


Allowing us to write the velocity profile as

$$ \bar{v}_z = 6.97 { v_0 d_0 \over z} f(\xi) $$

## Eddy Viscosity

The eddy viscosity model makes the assumption that the turbulent shear stress depends on the rate of strain in a manner that is analogous to laminar flow, with the constant of proportionality being the *eddy viscosity* &epsilon;:

$$ \bar{\tau}_{rz} = - \rho \varepsilon {\partial \bar{v}_z \over \partial r}$$


### Setting up the ODE

Recall that the equation of motion in the z direction is (in terms of the unitless function *F*)

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = {1 \over \rho} {\partial \over \partial r} \left( r \bar{\tau}_{rz} \right) $$

Substituting the expression for $ \bar{\tau}_{rz} $ we have

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = - \varepsilon {\partial \over \partial r} \left( r \left( \partial \bar{v}_{z} \over \partial r \right) \right) $$

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = - \varepsilon \left( \left( \partial \bar{v}_{rz} \over \partial r \right) + r \left( \partial^2 \bar{v}_{rz} \over \partial r^2 \right) \right) $$

Substituting in the expressions for $ {\partial \bar{v}_z \over \partial r} $ and $ {\partial^2 \bar{v}_z \over \partial r^2} $ we arrive at

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = { k \varepsilon \over z^2} \left( F^{\prime \prime \prime} - { F^{\prime \prime} \over \xi } + { F^{\prime} \over \xi^2 }  \right) $$

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = { k \varepsilon \over z^2} { d \over d \xi } \left( F^{\prime \prime} - { F^{\prime} \over \xi } \right) $$

at this point we note that *k* and &epsilon; have the same units of $ [[ \mathrm{length} ]]^2 \times [[ \mathrm{time} ]]^{-1} $ and are independent of *z* and &xi;, so we propose that $ \varepsilon = c k $ where *c* is some unknown constant of proportionality.

$$ \left( k \over z \right)^2 { d \over d \xi } \left( F F^{\prime} \over \xi \right) = c \left( k \over z \right)^2 { d \over d \xi } \left( F^{\prime \prime} - { F^{\prime} \over \xi } \right) $$

$$ { d \over d \xi } \left( F F^{\prime} \over \xi \right) = c { d \over d \xi } \left( F^{\prime \prime} - { F^{\prime} \over \xi } \right) $$

Integrating both sides

$$ { F F^{\prime} \over \xi } = c \left( F^{\prime \prime} - { F^{\prime} \over \xi } \right) + \mathrm{const}$$

By applying the boundary conditions we find the constant of integration is zero, thus

$$ F F^{\prime} = c \left( \xi F^{\prime \prime} -  F^{\prime} \right) $$

$${ d \over d \xi } \left( \frac{1}{2} F^2 \right) = c { d \over d \xi } \left( \xi F^{\prime} -  2 F \right) $$

Integrating both sides

$$ \frac{1}{2} F^2 = c \left( \xi F^{\prime} -  2 F \right) + \mathrm{const}$$

By applying the boundary conditions we find the constant of integration is zero, thus

$$ c \xi F^{\prime} = \frac{1}{2} F^2 -  2c F $$

Which is separable

$$ \int { d \xi \over \xi} = \int { c \over {\frac{1}{2} F^2 -  2c F} } dF $$

Integrating one last time

$$ \log \left( C_1 \xi \right) = \frac{1}{2} \log \left( F \over F + 4 c \right) $$

Where *C<sub>1</sub>* is an undetermined constant of integration. Re-arranging and solving for *F* we arrive at

$$ F\left( \xi \right) = { 4 c C_1 \xi^2 \over {1 - C_1 \xi^2 } } $$

A common substitution is $ C_1 = - \left( C_2 \over 2 \right)^2 $ then

$$ F\left( \xi \right) = { - c \left( C_2 \xi \right)^2 \over {1 + \frac{1}{4} \left( C_2 \xi \right)^2 } } $$

What we need, for the velocity profile, is the first derivative of *F*, which is

$$ F^{\prime}\left( \xi \right) = { - 2 c C_2^2 \xi \over \left( 1 + \left( C_2 \xi \over 2 \right)^2 \right)^2 } $$

and finally

$$ \bar{v}_z = -{k \over z} {F^{\prime} \over \xi} \\
 = -{k \over z} { 1 \over \xi }{ - 2 c C_2^2 \xi \over \left( 1 + \left( C_2 \xi \over 2 \right)^2 \right)^2 } \\
 = {2 \varepsilon C_2^2 \over z} \left( 1 + \left( C_2 \xi \over 2 \right)^2 \right)^{-2} $$
 
$$ f \left( \xi \right) = { \bar{v}_z \over \bar{v}_{z,max} } =  \left( 1 + \left( C_2 \xi \over 2 \right)^2 \right)^{-2} $$


```julia
f_ev(ξ; C₂=15.1) = ( 1 + (C₂*ξ/2)^2 )^-2
```


### Width at Half Height

Since we have a convenient closed form for the velocity profile, we can calculate what the parameter $C_2$ is in terms of the width at half height rather easily.

$$ f(\xi) =  \left( 1 + \left( C_2 \xi \over 2 \right)^2 \right)^{-2} $$

$$ \frac{1}{2} =  \left( 1 + \left( {C_2 \over 2} { b_{1/2} \over z }\right)^2 \right)^{-2} $$

$$ C_2 = 2 \sqrt{\sqrt{2}-1} \frac{z}{ b_{1/2} } $$

using the same parameterization as above we get


```julia
C₂ = 2*√(√(2)-1)/b_half
```


    15.179109738339214


### Velocity Profile

Returning to the momentum balance, we need to solve the integral:

$$ I = \int_{0}^{\infty} f\left( \xi \right)^2 \xi d \xi\\
= \int_{0}^{\infty} \xi \left( 1 + \left( C_2 \xi \over 2 \right)^2 \right)^{-4} d\xi $$

Which [can be integrated](https://www.wolframalpha.com/input?i2d=true&i=Integrate%5B%CE%BE*Power%5B%5C%2840%291%2BDivide%5B1%2C4%5DPower%5Bc*%CE%BE%2C2%5D%5C%2841%29%2C-4%5D%2C%7B%CE%BE%2C0%2C%E2%88%9E%7D%5D) to give
$$ I = {2 \over 3} C_2^{-2} $$

and finally

$$ k = \sqrt{ 3 \over 16 } C_2 v_0 d_0 $$

with the velocity profile as

$$ \bar{v}_z = 6.57 { v_0 d_0 \over z} \left( 1 + 57.6 \xi^2 \right)^{-2} $$


## Empirical Velocity Profiles

Perhaps the most widely used turbulent jet model is simply an empiirical gaussian fit to the data. These are easy to use -- no solving of ODEs required -- and fitting them to data is relatively straight forward. There is no real theoretical basis that I am aware of, merely based on the observation that a gaussian function fits the velocity profile well.

$$ f \left( \xi \right) = \exp \left( -c \xi^2 \right) $$

Where *c* is a parameter determined by fitting to a dataset.


```julia
f_emp(ξ; c=72) = exp(-c*ξ^2)
```


### Width at Half Height

Since we have a convenient closed form for the velocity profile, we can calculate what the parameter $c$ is in terms of the width at half height rather easily

$$ f(\xi) =  \exp \left( -c \xi^2 \right) $$

$$ \frac{1}{2} =  \exp \left( - c  \left( \frac{ b_{1/2} }{z} \right)^2 \right) $$

$$ c = \ln \left( 2 \right) \left( \frac{z}{ b_{1/2} } \right)^2 $$

using the same parameterization as above we get


```julia
c = log(2)/b_half^2
```


    96.39039423504045


### Velocity Profile

Returning to the momentum balance, we need to solve the integral:

$$ I = \int_{0}^{\infty} f\left( \xi \right)^2 \xi d \xi\\
= \int_{0}^{\infty} \xi \exp \left( -2 c \xi^2 \right) d\xi $$

Which [can be integrated](https://www.wolframalpha.com/input?i2d=true&i=Integrate%5B%CE%BE*Exp%5B-2*c*Power%5B%CE%BE%2C2%5D%5D%2C%7B%CE%BE%2C0%2C%E2%88%9E%7D%5D) to give

$$ I = {1 \over 4 c} $$

and finally

$$ k = \sqrt{ c \over 2 } v_0 d_0 $$


with the velocity profile as

$$ \bar{v}_z = 6.94 { v_0 d_0 \over z} \exp\left( -193.0 \xi^2 \right) $$


## Comparing the Models

At this point two models of velocity were derived using different models of the free turbulent stress and one purely empirical model was introduced. Each of these models uses a different set of parameters, and have different strengths and weaknesses in terms of usability. To compare them like-for-like we can scale each to the same width at half height, which is shown below along with some measured data[^pope]

We can also calculate a Mean Square Error (MSE) and evaluate which model is a better fit to the observed velocity profile.


[^pope]: Stephen B. Pope, *Turbulent Flows*, Cambridge University Press, Cambridge (2000) Points captured from a figure using [WebPlotDigitizer](https://automeris.io/WebPlotDigitizer)


![svg](/images/turbulent_jet_notes_files/output_46_0.svg)


    Prandtl Mixing Length Model MSE 0.00051
    Eddy Viscosity Model        MSE 0.00092
    Gaussian (empirical) Model  MSE 0.00054


Interestingly the Prandtl mixing length model works the best, though the gaussian fit is close enough as to be essentially the same given this data set. Which is convenient as a gaussian fit is easier to work with. The eddy viscosity model is the easiest to derive, however it clearly does not work as well for the outer parts of the jet.

The above approach, the one you will most likely see in the literature, compares each model scaled to the same height and width. Which is sensible if one is planning on fitting data, and allowing that the height and width to be free parameters. However we know, from the analysis above, that the height of each model is dependent upon the width, so might be instructive to look at how that plays out in practice.

Suppose we are looking at a velocity profile far enough downstream to be in the fully developed flow, say $ z = 7 d_0 $


![svg](/images/turbulent_jet_notes_files/output_48_0.svg)


Note that in the region near the center-line the three models are no longer particularly close to one another and the eddy viscosity and prandtl mixing length models have changed places. Relative to the predicted $v_{max}$ the the eddy viscosity model stays high when compared to the prandtl mixing length model, however the eddy viscosity model predicts a lower $v_{max}$ such that the effect is entirely reversed.

It's also worth noting that the gaussian fit and the prandtl mixing length model track one another reasonably well. I have a gaussian fit *of* the Tollmien tabulated results used in some papers when a smooth interpolation of the intermediate values is required and this suggests that may not be a bad idea. Though, to me, just solving the ode is easier, on a modern machine it takes milliseconds or less, and a good ode package like `DifferentialEquations.jl` provides a higher-order interpolation for free.

This comparison has been done with each of the model parameters set based on a shared width. However there are as many different ways of arriving at the model parameters as there are datasets to fit against. There is a wide spread in tabulated values in the literature and so the predictions of two independently arrived at models can be quite different due all of these factors coming together.

## Where to go from here

All of this work was to determine the *velocity* field, which is not necessarily what anyone cares about. In a release scenario, for example, it is concentration that is most relevant. For a heat transfer application, perhaps, you may care about the temperature field instead. However, with the velocity field the concentrations, temperatures, total entrained flow, etc. can be easily derived.


For a complete listing of code used to generate data and figures, please see the [corresponding julia notebook](https://nbviewer.org/github/aefarrell/aefarrell.github.io/blob/main/_notebooks/2022-04-08-turbulent_jet_notes.ipynb)
{: .notice--info}


---
