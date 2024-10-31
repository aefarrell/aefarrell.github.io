---
title: "Relief Valve Sizing with Real Gases"
toc: true
toc_label: "Contents"
toc_sticky: true
comments: true
categories:
  - examples
tags:
  - pressure relief
  - gas flow
  - equations of state
tagline: "compressible flow using equations of state"
header:
  overlay_image: /images/adiabatic_compressible_flow_files/pipes-unsplash-header.jpg
  overlay_filter: 0.25
  caption: "Photo by [Ricardo Gomez Angel](https://unsplash.com/es/@rgaleriacom?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText) on [Unsplash](https://unsplash.com/)"
---

# Relief Valve Sizing with Real Gases

Very often, in chemical engineering, the line between problems one can solve one's self, with code or more often than not in a spreadsheet, and problems that are solved with a piece of commercial software is when ideal fluid assumptions break down. Relief valve sizing presents a typical example: if the fluid is (approximately) an ideal gas then sizing is simple and often done in a spreadsheet. When this isn't the case, if the compressiblity is <0.8 or >1.1 <a href="#fn-1" class="sidenote-number"></a><span class="sidenote" id="fn-1">*[API 520](#API-2020)*, 68.</span>, then one typically has to turn to some commercial software. Models of real fluids are complicated and extracting the relevant thermodynamic properties from them can be quite tedious when doing it all from scratch.

The Julia ecosystem has several libraries that help bridge this gap, exposing a wide array of equations of state for real fluids as well as robust libraries for the integration and optimization methods needed to solve real fluid problems. This post walks through how to size a relief device, in gas service, starting from an ideal gas and working through various methods for real gases using equations of state.

## Sizing a Pressure Relief Valve

The general idea for sizing a relief valve is to determine the minimum area required such that the mass flow through the valve equals the mass flow required for the governing release case via the relation<a href="#fn-2" class="sidenote-number"></a><span class="sidenote" id="fn-2">*[API 520](#API-2020)*, equation B.5</span>


$$ W = K A G $$


where $W$ is the mass flow-rate, $K$ is a capacity correction, $A$ is the flow area, and $G$ the theoretical (friction-less) mass flux through the valve. Valves are sized based on the theoretical flow area.

The general process is as follows:

1. Determine the governing release rate, $W$
2. Determine the capacity correction, $K$
3. Calculate the mass flux through the valve, $G$
4. Calculate the theoretical flow area $A = \frac{W}{K G}$
5. Select the appropriate valve with a flow area $>A$.

Standards, such as [API 520](#API-2020), give equations that combine steps 3 and 4 and absorb unit-conversions into the constants, so that the equation is in a more convenient form, but this is what is happening under the hood.

The complications creep in through calculating $G$, it is path-dependent and is a function of the equation of state for the fluid. For gas releases the relief device is typically treated as an isentropic nozzle, the assumption being that the flow-rate through the valve is typically large enough that any heat transfer can be neglected.

<figure>
<img src="/images/relief_valve_sizing_files/pressure_relief.png" />
<figcaption> A hypothetical pressure relief device, connected to a pressure reservoir (1) and discharging into the atmosphere (2).</figcaption>
</figure>

Consider the differential form of the [mechanical energy balance](https://en.wikipedia.org/wiki/Bernoulli%27s_principle), along a streamline from the stagnation point, in the vessel, through the valve and out into the atmosphere, assuming no elevation change and no friction

$$ dP + \rho u du = 0 $$

$$ u du = -\frac{dP}{\rho} = - vdP $$

Integrating from the stagnation point to the throat of the nozzle gives

$$ \frac{1}{2} u_t^2 = - \int_{P_1}^{P_t} v dP $$

Where the velocity at the stagnation point, $u_1=0$. Putting this in terms of the mass flux $u = v G$

$$ \frac{1}{2} v_t^2 G_t^2 = - \int_{P_1}^{P_t} v dP $$

$$ G_t = \frac{1}{v_t} \sqrt{-2 \int_{P_1}^{P_t} v dP} = \rho_t \sqrt{2 \int_{P_t}^{P_1} v dP} $$

This integral cannot be solved directly at this point as the conditions at the throat of the nozzle are undefined. So solving this requires simultaneously solving for the nozzle conditions, $P_t, T_t$.

If we specify that the streamline follows an isentropic path, then we can construct a constrained maximization problem: the nozzle conditions are the $P_t$ and $T_t$ which maximizes $G_t$ where the integration is taken along an isentropic path.

### Choked Flow

In the case where flow is choked, i.e. the flow in the nozzle reaches sonic velocity, the maximum $G_t$ occurs at the sonic velocity with a pressure $P_t > P_2$. This can allow for the direct calculation of the mass flux as $G_t = \rho_t c_t$, where $c_t$ is the sonic velocity at the throat. No integration required.

## A Motivating Example

Consider the release of ethane from a vessel at 200 bar and 400 K, for the sake of simplicity assume the release is directly into the atmosphere at 1 bar and 288.15 K (15°C) (the flow is going to be choked, so this doesn't actually matter).

```julia
using Unitful.jl
```

```julia
begin
# the vessel properties
    P₁ = 200u"bar"
    T₁ = 400u"K"

# the ambient properties
    P₂ = 1u"bar"
    T₂ = 288.15u"K"
end
```

We can use `Clapeyron.jl` to initialize a few example equations of state for ethane. In this case I'm going to use an ideal gas model (`ReidIdeal` is an ideal gas model that also includes correlations for the ideal gas heat capacity), a cubic equation of state (volume translated Peng Robinson), and an empirical Helmholtz model (GERG-2008).

```julia
using Clapeyron
```

```julia
begin
# assorted equations of state for ethane
    ig_ethane = ReidIdeal(["ethane"])
    vtpr_ethane = VTPR(["ethane"]; idealmodel = ReidIdeal)
    gerg_ethane = GERG2008(["ethane"])
end
```

At system conditions ethane is a super critical fluid, with the temperature and pressure above the critical point, which can be modelled as a dense gas.

## The Ideal Gas Case

Considering the choked flow case, we know that $G_t = \frac{c_t}{v_t}$ and, for an ideal gas, the sonic velocity is given by<a href="#fn-3" class="sidenote-number"></a><span class="sidenote" id="fn-3">[Tilton](#tilton-2008), "Fluid and Particle Dynamics," equation 6-113.</span>

$$ c = \sqrt{ {k R T} \over M} = \sqrt{k P v} $$

Combining these we have

$$ G_t = \frac{c_t}{v_t} = { \sqrt{ k P_t v_t } \over v_t } = \sqrt{ k P_t \over v_t } $$

It can be shown that, along an isentropic path defined by $P v^k = \mathrm{const}$, the critical pressure ratio is<a href="#fn-4" class="sidenote-number"></a><span class="sidenote" id="fn-4">[Tilton](#tilton-2008), "Fluid and Particle Dynamics," equation 6-119.</span>

$$ {P_t \over P_1} = { P_{chk} \over P_1 } = \left(2 \over {k+1} \right)^{k \over {k-1} } $$

Which allows us to write

$$ P_t = P_1 {P_t \over P_1} = P_1 \left(2 \over {k+1} \right)^{k \over {k-1} } $$

and (using $P_1 v_1^k = P_t v_t^k$)

$$ v_t = v_1 \left(P_1 \over P_t \right)^{1 \over k} = v_1 \left(2 \over {k+1} \right)^{-1 \over {k-1} } $$

Substituting back into the equation for $G_t$<a href="#fn-5" class="sidenote-number"></a><span class="sidenote" id="fn-5">*[API 520](#API-2020)*, equation B.21</span>

$$ G_t = \sqrt{ k \frac{P_1}{v_1} \left(2 \over {k+1} \right)^{k+1 \over {k-1} } } $$

or, to put it in terms of density

$$ G_t = \sqrt{ k P_1 \rho_1 \left(2 \over {k+1} \right)^{k+1 \over {k-1} } } $$

where $k$, the isentropic expansion factor for an ideal gas, is the ratio of heat capacities

$$ k = { c_{p,ig} \over c_{v,ig} } $$

{% capture inlinenote-1 %}
This is the basis of [API 520 Part 1](#API-2020) equation 9 where the following substitutions is made:
	
$$ \rho = { {P M} \over {Z R T} } $$

and the constant $R$ and some unit conversions are rolled up into the constant 0.03948 in the expression for $C$

$$ R = 8.314 { {\mathrm{m^3} \cdot \mathrm{Pa} } \over {\mathrm{mol} \cdot \mathrm{K} } } = 8,314 { {\mathrm{m^3} \cdot \mathrm{Pa} } \over {\mathrm{kmol} \cdot \mathrm{K} } } = 8,314 { {\mathrm{kg} \cdot \mathrm{m^2} } \over { \mathrm{kmol} \cdot \mathrm{s^2} \cdot \mathrm{K} } } 	$$

$$ {1 \over \sqrt{R} } \left[ { \sqrt{ \mathrm{kmol} \cdot \mathrm{K} } \cdot \mathrm{s} } \over { \sqrt{\mathrm{kg} } \cdot \mathrm{m} } \right] \times 3600 \left[ \mathrm{s} \over \mathrm{h} \right] \times 10^{-6} \left[ \mathrm{m^2} \over \mathrm{mm^2} \right] \times 10^3 \left[ \mathrm{Pa} \over \mathrm{kPa} \right] \\= 0.03948 \left[ \sqrt{\mathrm{kmol} \cdot \mathrm{kg} \cdot \mathrm{K} } \over { \mathrm{h} \cdot \mathrm{mm^2} \cdot \mathrm{kPa} } \right] $$

{% endcapture %}

<div class="notice">
  {{ inlinenote-1 | markdownify }}
</div>

We can use `Clapeyron.jl` to calculate $k$ at any given temperature, using correlations for the ideal gas heat capacity.

```julia
function isentropic_expansion_factor(model::IdealModel, P, T; z=[1.0])
    cₚ_ig = isobaric_heat_capacity(model, P, T; phase=:vapor)
    cᵥ_ig = isochoric_heat_capacity(model, P, T; phase=:vapor)
    return cₚ_ig/cᵥ_ig
end
```

From which we calculate k= 1.146.

We can check our work by comparing with the tabulated values. At 15°C and 1 atm we calculate k= 1.193 which is the same as the tabulated value of 1.19<a href="#fn-6" class="sidenote-number"></a><span class="sidenote" id="fn-6">*[API 520](#API-2020)*, 70.</span> (given at 15°C and 1 atm).


```julia
# this is a hack, ideal models in Clapeyron do not return a molar weight
# and so cannot return a mass density
Clapeyron.mw(model::IdealModel) = Clapeyron.mw(vtpr_ethane)
```

```julia
function mass_flux_choked(model, P, T; z=[1.0])
    k = isentropic_expansion_factor(model, P, T; z=z)
    ρ = mass_density(model, P, T, z; phase=:vapor)
    Gₜ² = k*P*ρ*(2/(k+1))^((k+1)/(k-1))
    return √(Gₜ²)
end
```

The theoretical mass flux is then 38359 kg m^-2 s^-1

The ideal gas model, when the flow is choked, calculates the mass flux directly without needing to calculate the actual conditions at the nozzle. These can be calculated easily as well<a href="#fn-7" class="sidenote-number"></a><span class="sidenote" id="fn-7">[Tilton](#tilton-2008), "Fluid and Particle Dynamics," equations 6-119 and 6-120.</span>.

```julia
nozzle_pressure_ideal(P, T, k) = P*(2/(k+1))^(k/(k-1))
```

```julia
nozzle_temperature_ideal(P, T, k) = T*(2/(k+1))
```

The pressure at the nozzle is 115 bar the temperature at the nozzle is 373 K, which is above the critical point. The fluid super-critical and choked when leaving the PSV.


## The Isentropic Expansion Factor

At the vessel conditions, the VTPR model of ethane gives a compressibility factor of 0.672 (GERG-2008 model gives a similar value of 0.69), well below 0.8 and therefore outside the range where the ideal gas model is expected to work well.

In this case an alternative method is to calculate what the *effective* isentropic expansion factor would be, for the real gas, assuming that the real fluid obeys

$$ P_1 v_1^n = P_t v_t^n $$

where $n$ is a constant.

The derivation of $n$ follows from the definition of the speed of sound in a gas<a href="#fn-8" class="sidenote-number"></a><span class="sidenote" id="fn-8">[Tilton](#tilton-2008), "Fluid and Particle Dynamics," 6-22; [Gmehling *et al.*](#gmehling-2012), *Chemical Thermodynamics*, 113</span>

$$ c = \sqrt{ \left( {\partial P} \over {\partial \rho} \right)_S} =\sqrt{ -v^2 \left( {\partial P} \over {\partial v} \right)_S} $$

The constant entropy partial derivative can be re-written to eliminate entropy<a href="#fn-9" class="sidenote-number"></a><span class="sidenote" id="fn-9">[Gmehling *et al.*](#gmehling-2012), *Chemical Thermodynamics*, 660</span>

$$ \left( {\partial P} \over {\partial v} \right)_S = { { \left( {\partial S} \over {\partial T} \right)_P \left( {\partial P} \over {\partial T} \right)_v } \over { \left( {\partial S} \over {\partial T} \right)_v \left( {\partial v} \over {\partial T} \right)_P } } $$

Using the relations<a href="#fn-10" class="sidenote-number"></a><span class="sidenote" id="fn-10">[Gmehling *et al.*](#gmehling-2012), *Chemical Thermodynamics*, equations C.21 and C.8 (respectively)</span>

$$ c_p = T \left( {\partial S} \over {\partial T} \right)_P $$

and

$$ c_v =  T \left( {\partial S} \over {\partial T} \right)_V $$

we get

$$ \left( {\partial P} \over {\partial v} \right)_S = {c_p \over c_v} \left( {\partial P} \over {\partial v} \right)_T $$

and the sonic velocity is then

$$ c = \sqrt{ -v^2 {c_p \over c_v} \left( {\partial P} \over {\partial v} \right)_T } $$

equating this to the ideal gas case, $c = \sqrt{n P v}$, and solving for $n$ gives<a href="#fn-11" class="sidenote-number"></a><span class="sidenote" id="fn-11">*[API 520](#API-2020)*, equation B.13</span>

$$ n = -\frac{v}{P} {c_p \over c_v} \left( {\partial P} \over {\partial v} \right)_T $$

where $n$ has been used to distinguish it from $k$ (the ideal gas case). This is the version of $n$ presented in the typical references, such as [API 520](#API-2020). The derivation, however, hints at a useful shortcut to calculating $n$ that does not require digging into the internals of `Clapeyron.jl` to retrieve partial derivatives:

$$ n = { c^2 \over {P v} } = { {\rho c^2} \over P } $$

The remainder of the calculations are identical as the ideal gas case, simply substituting $n$ wherever $k$ appears. Unfortunately $n$ is not actually constant and depends on the temperature and pressure, which are not actually known in the nozzle, so the temperature and pressure at the stagnation point are often used instead.

```julia
function isentropic_expansion_factor(model, P, T; z=[1.0])
    ρ = mass_density(model, P, T, z)
    c = speed_of_sound(model, P, T, z)
    n = ρ*c^2/P
    return n
end
```

Using effective isentropic expansion factors from the VTPR equation of state, the theoretical mass flux is 57811 kg m^-2 s^-1 ( 59321 kg m^-2 s^-1 from GERG-2008 ). This is quite a bit larger than the ideal case, indicating that the ideal gas law leads to a significantly over-sized PRV, 51.0% larger when compared to sizing done using the isentropic expansion factor and the VTPR equation of state.

<figure>
<img src="/images/relief_valve_sizing_files/figure1.svg" />
<figcaption>The isentropic expansion factor for ethane at 400K, calculated for a range of stagnation pressures.</figcaption>
</figure>

The isentropic expansion factor method works best when $n$ is approximately constant over the isentropic path. As the above figure shows, this breaks down in ethane for pressures greater than ~100 bar. It also shows that the different equations of state start to diverge greatly further into the super-critical regime.


## Solving the Choked Flow Energy Balance

An alternative approach to direct integration, and one I have seen more often in older references, is to perform an energy balance over the isentropic path and, assuming the flow is choked, solve for sonic velocity in the nozzle<a href="#fn-12" class="sidenote-number"></a><span class="sidenote" id="fn-12">[Crowl *et al*](#crowl-2008), "Process Safety," 23-55; [Gmehling *et al.*](#gmehling-2012), *Chemical Thermodynamics*, 603.</span>. Consider an energy balance starting at the stagnation point, (1), and following an isentropic path to immediately after the throat of the nozzle (t).

$$ h_1 = h_t + \frac{1}{2} c_t^2 $$

Where $c_t$ is the speed of sound at nozzle, a function of $P_t$ and $T_t$. The procedure is then to solve the system of equations given by the energy balance and the entropy balance, $s_1 = s_t$, for $P_t$ and $T_t$, then the theoretical mass flux is given by

$$ G_t = \rho_t c_t $$

There are a few ways this could be done, but one simple way is to divide the problem into two:

+ Define the isentropic path, i.e. find the isentropic temperature for a given pressure *P*

```julia
using Roots, ForwardDiff
```

```julia
function isentropic_temperature(model, P, s, T0; z=[1.0])
	return find_zero( (x) -> s - entropy(model, P, x, z), T0)
end
```

+ Use the energy balance to solve for the pressure, following the isentropic path.

```julia
# this is a hack, Clapeyron's molecular weight function 
# does not have units
# using multiple dispatch, units can be added as needed
molecular_weight(model,z,::Number) = Clapeyron.molecular_weight(model,z)
```

```julia
molecular_weight(model,z,::Quantity) = Clapeyron.molecular_weight(model,z)*1u"kg"
```

```julia
function choked_energy_balance(model, P, s, T0; z=[1.0])
    T = isentropic_temperature(model, P, s, T0; z=z)
    c = speed_of_sound(model, P, T, z)
    h = enthalpy(model, P, T, z)
    Mw = molecular_weight(model,z,P)
    return h/Mw + 0.5*c^2
end
```

```julia
function choked_energy_balance_pressure(model, s, h, P0, T0; z=[1.0])
    return find_zero( (p) -> h - choked_energy_balance(model, p, s, T0), P0)
end
```

Calculating the mass flux then merely wraps these two steps in some set-up:
1. use the isentropic expansion factor to generate some initial guesses
2. calculate the initial entropy and specific enthalpy
3. solve for the pressure at the nozzle using the energy balance
4. find the density and sonic velocity at nozzle conditions
5. calculate the theoretical mass flux

```julia
function mass_flux_choked_energy_balance(model, P, T; z=[1.0])
    # estimate nozzle conditions using the expansion factor
    n = isentropic_expansion_factor(model, P, T; z=z)
    T0 = nozzle_temperature_ideal(P, T, n)
    P0 = nozzle_pressure_ideal(P, T, n)

    # calculate the entropy and specific enthalpy at initial conditions
    s₁ = entropy(model, P, T)
    h₁ = enthalpy(model, P, T)/molecular_weight(model, z, P)

    # find the pressure and temperature such that the enthalpy
    # and energy balance
    Pₜ = choked_energy_balance_pressure(model, s₁, h₁, P0, T0; z=z)
    Tₜ = isentropic_temperature(model, Pₜ, s₁, T0; z=z)

    # mass flux is the sonic velocity at nozzle conditions
    cₜ = speed_of_sound(model, Pₜ, Tₜ, z)
    ρₜ = mass_density(model, Pₜ, Tₜ, z)
    return ρₜ*cₜ
end
```

Solving the choked flow energy balance, using VTPR equation of state, the theoretical mass flux is 54629 kg m^-2 s^-1 ( 54353 kg m^-2 s^-1 from GERG-2008 ). This is also quite a bit larger than the ideal case, 42.0% larger. Though the values for the two equations of state are closer, indicating that this method is less sensitive to model choice.

An alternative way of doing this would be to solve for $P_t$ and $T_t$ simultaneously. However, one advantage of doing it as independent steps is that the isentropic path is available, outside of calculating the mass flux. Which can be useful for other calculations, or simply to plot it on a phase diagram as is shown below.

<figure>
<img src="/images/relief_valve_sizing_files/figure2.svg" />
<figcaption>The isentropic paths for the ideal gas, effective isentropic factor, and true isentropic path methods.</figcaption>
</figure>

In this case the ideal gas method and the isentropic expansion factor method bracket the more exact method of solving the energy balance directly.

As it is written, this method would need to be modified to allow for non-choked flow. This is done by eliminating the assumption $u_t = c_t$ and instead finding the conditions which maximize $G_t$ (while still satisfying the energy balance). This will arrive at the same solution, in the case of choked flow, but with a little more effort.


## Direct Integration

Direct integration is the method most commonly recommended today, as it is entirely general. It can be used to solve all flow conditions from liquids to gases as well as two-phase mixtures. As a reminder, this method constitutes finding the $P_t$ and $T_t$ that maximize the mass flux given by

$$ G_t = \rho_t \sqrt{2 \int_{P_t}^{P_1} v dP} $$

Since the isentropic path was already defined, we can easily generate an expression for the specific volume as a function of pressure, $v(P)$

```julia
function isentropic_specific_volume(model, P, s, T0; z=[1.0])
    T = isentropic_temperature(model, P, s, T0; z=z)
    V = volume(model, P, T, z)
    MW = molecular_weight(model, z, P)
    return V/MW
end
```

The function we wish to optimize is the numerical integral, from $P_t$ to $P_1$, where $P_t$ is the optimization variable. In this case it is $-G$ since the standard form for optimization is minimizing.

Below I define a factory function since the objective function is a function of only one variable, $P_t$, but requires a lot of parameters, this simplifies the problem set-up.

```julia
using QuadGK
```

```julia
function objective_factory(model, P, T; z=[1.0])
    # estimate nozzle conditions using the expansion factor
    n = isentropic_expansion_factor(model, P, T; z=z)
    T0 = nozzle_temperature_ideal(P, T, n)

    # calculate the entropy initial conditions
    s₁ = entropy(model, P, T)

    function objective(Pₜ)
        Tₜ = isentropic_temperature(model, Pₜ, s₁, T0; z=z)
        ρₜ = mass_density(model, Pₜ, Tₜ, z)
        ∫vdP, err = quadgk( (p) -> isentropic_specific_volume(model, p, s₁, T0; z=z),
                            Pₜ, P)
        return -ρₜ*√(2∫vdP)
    end

    return objective
end
```

To better visualize what is going on we can plot the theoretical mass flux, $G$, as a function of nozzle pressures, $P_t$, and can clearly see a maximum around 100 bar, which matches what we expected from the previous method.

<figure>
<img src="/images/relief_valve_sizing_files/figure3.svg" />
<figcaption>The theoretical mass flux as a function of nozzle pressure.</figcaption>
</figure>

Using `Optim.jl` to perform the optimization with the Brent method we can solve the problem easily enough.

```julia
using Optim
```

```julia
function mass_flux_direct_integration(model, P_1, T_1, P_2; z=[1.0])
    obj = objective_factory(model, P_1, T_1; z=z)
    res = optimize(obj, P_2, P_1, Brent())
    return -1*minimum(res)
end
```

```julia
# Optim.jl doesn't play nicely with Unitful.jl
# this wrapper clears this up
function mass_flux_direct_integration(model, P_1::Quantity, T_1::Quantity,
                                      P_2::Quantity; z=[1.0])
    P_1 = ustrip(u"Pa", P_1)
    P_2 = ustrip(u"Pa", P_2)
    T_1 = ustrip(u"K", T_1)
    G = mass_flux_direct_integration(model, P_1, T_1, P_2; z=z)
    return G*1u"kg/m^2/s"
end
```

Direct integration of the VTPR equation of state gives a theoretical mass flux of 54629 kg m^-2 s^-1 ( 54353 kg m^-2 s^-1 from GERG-2008 ). Which is exactly the same as from solving the choked flow energy balance, as expected.

This method could be extended to include liquid and two-phase flows simply through re-writing the `isentropic_temperature` function to be phase-aware and utilize isothermal flash routines. Which shows its flexibility as a method. As it currently stands, this method only handles gases *but*, unlike the energy balance method, the flow does not have to be choked. If the flow is not choked, the maximum will occur at $P_2$ and whatever the isentropic temperature is at that point, the result will simply pop out without any extra effort.

## Comparing the Results

For the sake of completeness, there are two other methods that should be looked at, which are really special cases:
1. the ideal gas case, but using the real compressibility, $Z$, at stagnation conditions, this is the [API 520](#API-2020) standard approach for gases
2. using the isentropic expansion factor, *n* factor, method but calculating *n* at the average of the stagnation and nozzle conditions

These two approaches do better than the basic methods I presented, but I don't think they add enough value on their own. If you have a model of the gas that can give you the compressibility, for example, you are probably better off using either the isentropic expansion factor method or the direct integration method over correcting the ideal gas case (*you* aren't doing the calculations, the computer is, and what's the difference of a few milliseconds running the code?). Once a viable equation of state is in hand, the simplifications are not saving any *actual engineer doing their job* time, they are saving fractions of a second of compute time at the expense of (potentially) significant error.

I think the choice between the first law energy balance and the direct integration technique is more a matter of taste, at least in the case of choked flow. The direct integration method is in the relevant engineering codes/standards, and that is a strong justification for using it.

<figure>
<img src="/images/relief_valve_sizing_files/figure4.svg" />
<figcaption>A comparison of calculated theoretical mass flux for the six methods. The results from the first law energy balance and direct integration are identical.</figcaption>
</figure>

In this case the choice of equation of state did not matter strongly, just for fun I have included a few other common cubic equations of state, they all perform reasonably. However this example is for a single compound that is not strongly associating, it is the type of example where cubic equations of state should work well. The choice of equation of state will be far more important with mixtures and strongly associating substances.

### Opportunities for Optimization

I favoured readability, simplicity, and interoperability with `Unitful`, over performance in several places. For more production ready code there are a few pieces of low-hanging fruit:
1. `Clapeyron.jl` calculates most things, internally, in terms of $v$ and $T$, the convenience functions in terms of $P$, $T$ actually calculate $v$ internally as an intermediary value. This means the code above is, behind the scenes, redundantly calculating $v$ a lot. The functions could be re-written to use the volume forms of the equations, for example `entropy(model, P, T)` replaced with `Clapeyron.VT_entropy(model, V, T)`.
2. The direct integration method is integrating the same part of the volume curve multiple times. Every time the objective function is called it calculates the whole integral as if none of that had ever been done before. There are a lot of opportunities for caching parts of the integral to eliminate many function calls.
3. The optimization method, Brent, is not making use of gradients, I picked it more-or-less because it works and performance is fast enough that I don't notice. Gradient methods are more performant, and there are no doubt many other ways of goosing the performance here. Instead of a bounded method, for example, something like Newton's method could be used using `ForwardDiff` to find the gradient analytically and searching for where it is zero.

## Final Thoughts

I have long been an advocate for engineering to move out of using spreadsheets for everything and to use scripting languages and notebooks like [Jupyter](https://jupyter.org/) and [Pluto](https://plutojl.org/) far more. There are large classes of problems that are easy to solve with code and hard to solve with a spreadsheet. I think almost any calculation using equations of state fit into that category. We end up beholden to commercial software suppliers for calculations that, in my view, engineers should be doing themselves.

Presumably you could do the calculations I laid out above in Excel, at enormous effort, and making liberal use of the solver. Julia, however, has a robust ecosystem for doing all the complicated math, all that really needed to be done was connecting it up. What remains, for the engineer, is assessing the physical system and picking the appropriate methods and thermodynamic models.


For a complete listing of code used to generate data and figures, please see the [corresponding pluto notebook](https://github.com/aefarrell/aefarrell.github.io/blob/main/_notebooks/2024-10-28-relief_valve_sizing.jl)
{: .notice--info}

<h2>References</h2>
<ul>
<li><em><a name="API-2020">API</a> Standard 520, Sizing, Selection, and Installation of Pressure-relieving Devices, Part I - Sizing and Selection</em>. 10th ed. Washington, DC: American Petroleum Institute, 2020.</li>
<li><a name="CCPS-2017">Center for Chemical Process Safety</a>. <em>Guidelines for Pressure Relief and Effluent Handling Systems</em>. 2nd ed. Hoboken, NJ: John Wiley &amp; Sons, 2017.</li>
<li><a name="crowl-2008">Crowl</a>, Daniel A., Lawrence G. Britton, Walter L. Frank, Stanley Grossel, Dennis Hendershot, W. G. High, Robert W. Johnson <em>et al</em>. "Process Safety." in Green, <em>Perry’s Chemical Engineers’ Handbook</em>.</li>
<li><a name="green-2008">Green</a>, Don W., ed. <em>Perry’s Chemical Engineers’ Handbook</em>. 8th ed. New York: McGraw Hill, 2008.</li>
<li><a name="gmehling-2012">Gmehling</a>, J&uuml;rgen, B&auml;rbel Kolbe, Michael Kleiber, and J&uuml;rgen Rary. <em>Chemical Thermodynamics for Process Simulation</em>. Weinheim, DE: Wiley-VCH Verlag &amp; Co., 2012</li>
<li><a name="tilton-2008">Tilton</a>, James N. "Fluid and Particle Dynamics." in Green, <em>Perry’s Chemical Engineers’ Handbook</em>.</li>
</ul>
