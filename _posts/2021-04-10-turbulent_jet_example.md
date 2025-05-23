---
title: "Turbulent Jet Example - Acetylene Leak"
last_modified_at: 2023-12-26
toc: true
toc_label: "Contents"
toc_sticky: true
comments: true
categories:
  - examples
tags:
  - chemical releases
  - dispersion modelling
  - hazard screening
  - turbulent jets
tagline: "Estimating the explosive mass"
header:
  overlay_image: /images/chuttersnap-unsplash-header.jpg
  overlay_filter: 0.25
  caption: "Photo by [CHUTTERSNAP](https://unsplash.com/@chuttersnap?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText) on [Unsplash](https://unsplash.com/)"
---


# Turbulent Jet Example

In previous examples I discussed release scenarios involving vapour clouds spreading over a large area, carried by the wind. In those examples the momentum of the jet of fluid was not very important relative to the ambient wind conditions and could be ignored. In this example I am looking at the opposite extreme, a release from a pressure vessel inside a building where the momentum of the jet dominates.

## The Scenario

Consider, for an example, a leak from an acetylene cylinder inside a large building, such as in a warehouse or shop. We imagine, for convenience, that the air within the building is quiescent. For the sake of an example suppose the leak is a 1/4 in. hole, similar in diameter to a typical acetylene hose, and that the operating pressure at that point is 15psig<a href="#fn-1" class="sidenote-number"></a><span class="sidenote" id="fn-1">From CGA G-1 2009 the safe operating pressure of an acetylene system</span> We are interested in exploring the concentration distribution as the acetylene jets into the air and mixes, with our reference concentration of interest being half the LEL of 2.5%(vol).


```julia
using Unitful: @u_str, ustrip

inch = ustrip(u"m", 1u"inch") # unit conversion inch->m
psi = ustrip(u"Pa", 1u"psi")  # unit conversion psi->Pa

p₂ = 14.7psi   # atmospheric pressure, Pa absolute
T₂ = 25+273.15 # ambient temperature, K

d  = 0.25inch  # diameter of the hole, m
p₁ = 15psi+p₂  # pressure of the acetylene, Pa absolute
T₁ = T₂        # the release temperature, K
```

We can look up some properties of acetylene in Perry's<a href="#fn-2" class="sidenote-number"></a><span class="sidenote" id="fn-2">[Poling *et al*](#poling-2007) "Physical and Chemical Data".</span>


```julia
# universal gas constant, J/mol/K
R = 8.31446261815324 

# ideal gas density, kg/m³
ρ(p,T;MW) = (p*MW)/(R*T)/1000

# gas viscosity correlation, Pa*s
μ(T;C) = (C[1]*T^(C[2]))/(1+(C[3]/T)+(C[4]/T^2)) 

# Properties of Acetylene
MWⱼ = 26.037 # molar mass, kg/kmol
LEL = 0.025  # Lower explosive limit, vol/vol
k   = 1.26   # ratio cp/cv at 15C
μⱼ  = μ(T₁;C=[1.2025e-6,0.4952,291.4,0])
ρ₁  = ρ(p₁,T₁;MW=MWⱼ)

# Properties of Air
MWₐ = 28.960  # molar mass, kg/kmol
ρ₂  = ρ(p₂,T₂;MW=MWₐ)
```

## The Release Rate

We can model the release as a gas jet<a href="#fn-3" class="sidenote-number"></a><span class="sidenote" id="fn-3">[AIChE/CCPS](#ccps-1999) *Guidelines for Consequence Analysis*, 29.</span> where the gas is ideal and the expansion through the jet is an isentropic process<a href="#fn-4" class="sidenote-number"></a>

<span class="sidenote" id="fn-4">[AIChE/CCPS](#ccps-1999) *Guidelines for Consequence Analysis* has a mistake in equation 2.16, the version given here is correct</span>

$$ G = \rho u = c_d \sqrt{ \rho_1 p_1 \left( 2 k \over k-1 \right) \left[ \left(p_2 \over p_1\right)^{2 \over k} - \left(p_2 \over p_1\right)^{k+1 \over k} \right]} $$


for non-choked flow and

$$ G = c_d \sqrt{ \rho_1 p_1 k \left( 2 \over k+1 \right)^{k+1 \over k-1} } $$

for choked flow, which occurs when

$$ \left(p_2 \over p_1 \right) \lt \left( 2 \over k+1 \right)^{k \over k-1} $$

Where *G* is the mass velocity of acetylene discharged through the hole (in kg/m²/s), *c<sub>d</sub>* is the discharge coefficient which can be assumed to be 0.61<a href="#fn-5" class="sidenote-number"></a><span class="sidenote" id="fn-5">[AIChE/CCPS](#ccps-1999) *Guidelines for Consequence Analysis*, 30.</span>, and the rest are as defined earlier. I am assuming, here, that the hole is circular for simplicity.



```julia
(p₂/p₁) < (2/(k+1))^(k/(k-1)) 
```




    true



Therefore the flow is choked and


```julia
c_d = 0.61

G = c_d * √(ρ₁*p₁*k*(2/(k+1))^((k+1)/(k-1)) )
```




    267.1556913840265



The density at the orifice is reduced, through the expansion and, for an isentropic process, is related to the pressure by

$$ {\rho_o \over \rho_1} = \left( p_o \over p_1 \right)^{1 \over k} $$

Where subscript *o* indicates at the orifice. At this point, after the expansion $p_o = p_2$ and

$$ \rho_o = \rho_1 \left( p_o \over p_1 \right)^{1 \over k} $$


```julia
ρₒ = ρ₁*(p₂/p₁)^(1/k)
```




    1.2307940295609565



The velocity at the orifice, i.e. after the gas has expanded, is then

$$ u_o = {G \over \rho_o} $$


```julia
uₒ = G/ρₒ
```




    217.05962571115586



## Jet Behavior

To model the concentration profile I am going to assume a turbulent jet, from a circular hole, mixing with air. In this case the density of air and acetylene are similar and so a simple turbulent jet model is appropriate. If there was a significant difference in densities then a density correction would be needed, however for many applications "close" means a ratio of ambient to jet densities between<a href="#fn-6" class="sidenote-number"></a><span class="sidenote" id="fn-6">[Poleshaw and Golub](#poleshaw-2023) "Jets".</span>

$$ \frac{1}{4} \le { \rho_{a} \over \rho_{j} } \le 4 $$

Where subscript *a* indicates the ambient fluid and *j* the jet.

Circular turbulent jets expand by entraining ambient fluid, tracing out a cone defined by a jet angle $\alpha \approx 15-25^\circ$. The mixing layer penetrates into the jet forming the potential cone, inside is pure jet material and outside is mixed. After approximately 6 hole diameters the region is fully developed.<a href="#fn-7" class="sidenote-number"></a><span class="sidenote" id="fn-7">[Revill](#revill-1992) "Jet mixing".</span>


![image.png](/images/turbulent_jet_example_files/att1.png)

Empirical approximations of the velocity, and concentration, profiles are often given with respect to this jet angle or, equivalently, the slope of line (i.e. $\tan \frac{\alpha}{2}$)

Another important factor is the Reynolds number, the jet is fully turbulent when $Re \gt 2000$, where the Reynolds number is calculated with respect to the initial jet velocity and jet diameter (i.e. the hole diameter)

$$ Re = { \rho u d \over \mu } = { G d \over \mu }$$


```julia
0.25 < (ρ₂/ρₒ) < 4
```




    true




```julia
Re = G*d/μⱼ

Re > 2000
```




    true



The densities are within the appropriate range and the flow is fully turbulent, so the turbulent jet model requirements are satisfied.


### Velocity and Concentration distributions

There are many different empirical velocity distributions as well as velocity distributions derived from theories of turbulent mixing available in various references. Mostly of the same general type (gaussian), but parametrized slightly differently. However, in my experience, there are far fewer concentration distributions available, this is not too critical due to an interesting result in turbulent mass transfer for jets<a href="#fn-8" class="sidenote-number"></a><span class="sidenote" id="fn-8">[Bird, Stewart, and Lightfoot](#bird-2007) *Transport Phenomena*, 416.</span>


$$ { C \over C_{max} } = \left( v_z \over v_{z,max} \right)^{Sc_t} $$

That is, at a given distance *z* away from the hole, the concentration profile is the velocity profile raised to the power $Sc_t$ -- the turbulent Schmidt number. Experimentally this is approximately 0.7. Note also that $C_{max}$ and $v_{z,max}$ are taken at the centerline. Physically this means that the concentration profile, at a given downstream distance, is wider than the velocity distribution; concentration expands more.

A similar way of capturing the same phenomenon that is often seen with empirical velocity distributions is to define a width parameter $b$ and note that the equivalent width for the concentration profile is $1.17b$<a href="#fn-9" class="sidenote-number"></a><span class="sidenote" id="fn-9">[Kaye, Khan, and Testik](#kaye-2018) "Environmental Fluid Mechanics".</span> and substitute in accordingly.

In this example I am using the empirical concentration given in [Lees'](#lees-1996)<a href="#fn-10" class="sidenote-number"></a><span class="sidenote" id="fn-10">[Lees](#lees-1996) *Loss Prevention*, 15/140.</span> for simplicity

$$ {C \over C_0 } = k_2 \left( d_h \over z \right) \left( \rho_z \over \rho_o \right)^{0.5} \exp \left( - \left( k_3 r \over z \right)^2 \right) $$

Note also the ratio of densities, the density $\rho_z$ is the density of the jet at some distance *z* and it is common to conservatively take this as $\rho_a$.

The parameters $k_2$ and $k_3$ are empirically derived for the particular jet and $k_2$ is a function of Reynolds number below $ Re \lt 20000 $ <a href="#fn-11" class="sidenote-number"></a><span class="sidenote" id="fn-11">[Long](#long-1963) "Estimation of the Extent of Hazard Areas Around a Vent".</span>. The conservative values suggested are *6* and *5* respectively.


```julia
function C(r, z; C₀=1.0, k₂=6, k₃=5, d=d, ρz=ρ₂, ρₒ=ρₒ)
    C = C₀ * k₂ * (d/z) * √(ρz/ρₒ) * exp(-(k₃*r/z)^2)
end
```

    
![svg](/images/turbulent_jet_example_files/output_19_0.svg)
    



At this point it is worth pointing out that the model of the jet is independent of the discharge rate. The concentration profile is only a function of the hole diameter and the fluid density. The velocity in the jet, and the amount of air entrained in the jet, do depend strongly on the initial discharge rate but in such a way that the concentration does not. As the jet velocity increases proportionally more air is entrained and the concentration profile remains constant.



    
![svg](/images/turbulent_jet_example_files/output_21_0.svg)


## Explosive Mass

Now that we have a model of the jet, showing the concentration of acetylene, the most relevant parameter we would want to know is the explosive mass such that some blast modeling could be done.

The most obvious way to do this is to integrate over the jet, using cylindrical coordinates for convenience

$$ m_e = \int \rho C(r,z) dV = 2\pi \rho_o \int_{0}^{\infty} \int_{0}^{\infty} C(r,z) r dr dz $$

Except that we define the explosive mass to be the volume where $ C > \frac{1}{2} LEL $. A lazy way to do this is to define a function that equals $C$ if it is $ \gt \frac{1}{2} LEL $ and zero otherwise.

The potential core region is poorly described by this model, and the closer to the origin of the jet the more un-physical the results: giving concentrations greater than 100% and being undefined completely at the origin. One way of hand waving this away is to chop off any concentrations above 100%.


```julia
function igrd(v; lim=0.5*LEL)
    r, z = v
    
    if z>0
        c = C(r,z)
        c = c<lim ? 0 : min(1,c)
    else
        c = 0
    end

    return r*c
end
```



Integrating over some plausible bounds, taken by looking at the plots above, gives the volume of acetylene.


```julia
using HCubature: hcubature

I, err = hcubature(igrd, [0, 0], [0.25, 2.0], atol = 1e-8)
```




    (0.0008207940258726464, 9.999922827914883e-9)



Which can be plugged into the equation to calculate the final explosive mass.


```julia
mₑ = 2*π*ρₒ*I
```




    0.006347452155224944



To give a sense of how much this is, the explosive mass is equivalent to ~1s of discharge at the steady state discharge rate.


```julia
m = G*(π/4)*d^2

mₑ/m
```




    0.7502356087241902



## Conclusions

Turbulent jet mixing is a much simpler model for estimating releases, especially when using empirical models, compared to models for plumes influences by buoyancy and wind. There are much fewer parameters that need to be estimated.

One big weakness to the model as presented here is that it does not take into account the enclosed space. If the assumption is that the warehouse is large and ignition sources are numerous then that likely doesn't matter, the acetylene leak will ignite before it has a chance to accumulate. However it will grossly underestimate the potential explosive mass that could develop as the acetylene disperses through the air of warehouse, since the model presumes the ambient air has no acetylene in it and is effectively infinite in extent.

This limitation would, for me, motivate exploring more detailed models of gas build up in enclosed spaces



For a complete listing of code used to generate data and figures, please see the [corresponding julia notebook](https://github.com/aefarrell/aefarrell.github.io/blob/main/_notebooks/2021-04-10-turbulent_jet_example.ipynb)
{: .notice--info}

## References

+ <a name="ccps-1999">AIChE/CCPS</a>. *Guidelines for Consequence Analysis of Chemical Releases.* New York: American Institute of Chemical Engineers, 1999.
+ <a name="bird-2007">Bird</a>, R. Byron, Warren E. Stewart, and Edwin N. Lightfoot. *Transport Phenomena*. 2nd ed. Hoboken: John Wiley & Sons, 2007. [archive](https://archive.org/details/transportphenome0000bird_n8h5)
+ <a name="kaye-2018">Kaye</a>, Nigel B., Abdul A. Khan, and Firat Y. Testik. "Environmental Fluid Mechanics" in *Handbook of Environmental Engineering*. Edited by Myer Kutz. New York: John Wiley & Sons, 2018.
+ <a name="lees-1996">Lees</a>, Frank P. *Loss Prevention in the Process Industries*. 2nd ed. Oxford: Butterworth-Heinemann, 1996.
+ <a name="long-1963">Long</a>, V.D. "Estimation of the Extent of Hazard Areas Around a Vent," *Second Symposium On Chemical Process Hazards* (1963): 6-14
+ <a name="poleshaw-2013">Poleshaw</a>, Yury V., and V.V. Golub. "Jets" in *[Thermopedia](https://www.thermopedia.com/content/903/)* 2013. [doi:10.1615/AtoZ.j.jets](https://dx.doi.org/10.1615/AtoZ.j.jets)
+ <a name="poling-2007">Poling</a>, Bruce E., George H. Thomson, Daniel G. Friend, Richard L. Rowley, and W. Vincent Wilding. "Physical and Chemical Data" in *Perry's Chemical Engineers' Handbook*. 8th ed. Edited by Don W. Green. New York: McGraw Hill, 2007.
+ <a name="revill-1992">Revill</a>, B. K. "Jet mixing" in *Mixing in the Process Industries*. 2nd ed. Edited by N. Harnby, M.F. Edwards, and A. W. Nienow. Oxford: Butterworth-Heinemann, 1992.
