### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> title = "Vessel Blowdown"
#> date = "2025-01-24"
#> description = "Pressure vessel blowdown models for ideal gases"
#> 
#>     [[frontmatter.author]]
#>     name = "Allan Farrell"
#>     url = "https://aefarrell.github.io"

using Markdown
using InteractiveUtils

# ╔═╡ 320df4b7-f6a8-4d84-a045-f1efc055e417
using OrdinaryDiffEq, Plots

# ╔═╡ b18f8cd2-d697-11ef-11a8-07514c0e9930
md"""
# Vessel Blowdown - Ideal Gases

A recurring task of mine is to look at some old calculations, done by some previous engineer whose identity is lost to time and organizational flux, and update them to match current reality. Depending on the state of the spreadsheet, and its lack of documentation, this can also mean going down a rabbit hole of research to find where, exactly, a given equation came from and what all the constants in it represent. This post is the result of one of those journeys, trying to track down the source of a model for depressuring a vessel.

![figure 1](figures/pressure_relief.png)

Consider the blow-down of a pressure vessel to a vent stack, wehre the vessel contains a gas. What we want, in this case, is the time it takes to depressure and the pressure in the vessel as a function of time -- the blowdown curve. As a first approximation we can consider the ideal gas case and examine two limiting behaviours for the vessel: when the walls are perfect insulators (the adiabatic case) and when the walls are perfect conductors of heat (the isothermal case). Furthermore we assume the blowdown is through an isentropic nozzle.
"""

# ╔═╡ 4464df47-970c-4fb4-aa46-aad7d2d34cc1
md"""
## The Adiabatic Case

The adiabatic case is often a good approximation for small vessels and early in the blowdown, when the rate of energy lost from the vessel through the bulk transport of the gas is much higher than any heat gained from the environment.

Starting with a mass balance on the vessel:

```math
\frac{dm}{dt} = - w
```

where *m* is the mass inside the vessel and *w* is the mass flow through the valve. Since the volume of the vessel is a constant, *V*, we can write the mass balance as

```math
V \frac{d \rho}{dt} = - w
```

We can perform a change of variables from *ρ* to *P*

```math
V \left( \frac{\partial \rho}{\partial P} \right)_S \frac{dP}{dt} = - w
```

The partial derivative is taken along an isentropic path as the adiabatic expansion within the vessel is isentropic (not because the valve is isentropic).

We can write the mass flow through the nozzle in terms of the theoretical, frictionless, mass velocity *G*, the discharge coefficient $c_D$, and the flow area *A*.

```math
w = c_D A G
```

giving[^1]

```math
\frac{dP}{dt} = - \frac{c_D A}{V} \left( \frac{\partial P}{\partial \rho} \right)_S G
```

[^1]: [Botros, Jungowski, and Weiss](#botros-1989), "Models and Methods of Simulating Gas Pipeline Blowdown."


"""

# ╔═╡ 0cce5b25-f6a0-4df6-8b58-89d68c5b3abc
md"""
### Fully Choked Flow

Assuming the flow through the valve is choked, the velocity in the throat is the sonic velocity $a_t$ and, for an ideal gas, the speed of sound is given by

```math
a = \sqrt{ {k R T} \over M} = \sqrt{ {k P} \over \rho}
```

An ideal gas undergoing an adiabatic expansion from vessel pressure to the pressure in the throat of the valve has the following relationship between density and pressure[^2]

```math
\frac{\rho_t}{\rho} = \left( \frac{P_t}{P} \right)^{\frac{1}{k}}
```

and, for choked flow, the pressure ratio is at maximum at[^3]

```math
\frac{P_t}{P} = \left( { 2 \over {k+1} } \right)^{\frac{k}{k-1}}
```

putting this all together we can write *G* in terms of vessel conditions $\rho$ and $P$

```math
G = \rho_t u_t = \rho_t a_t = \rho_t \sqrt{ {k P_t} \over \rho_t} = \sqrt{k \rho_t P_t}
```

```math
G = \sqrt{k \rho P} \left( 2 \over {k+1} \right)^{\frac{k+1}{2 \left( k - 1 \right)}}
```

From [thermodynamics](https://aefarrell.github.io/2024/10/28/relief_valve_sizing/#fn-8) we know

```math
\left( \frac{\partial P}{\partial \rho} \right)_S = a^2 = \frac{k P}{\rho}
```

and we can put this all together to get

```math
\frac{dP}{dt} = - \frac{c_D A}{V} \left( {k P} \over \rho \right) \sqrt{k \rho P} \left( 2 \over {k+1} \right)^{\frac{k+1}{2 \left( k - 1 \right)}}
```

At this point, it is standard to introduce a time constant $\tau$

```math
\tau = \frac{m_0}{w_0} = \frac{\rho_0 V}{c_D A \sqrt{k \rho_0 P_0}} \left( 2 \over {k+1} \right)^{-\frac{k+1}{2 \left( k - 1 \right)}}
```

or, more clearly,

```math
\frac{1}{\tau} = \frac{c_D A}{V} \sqrt{ {k P_0} \over \rho_0 } \left( 2 \over {k+1} \right)^{\frac{k+1}{2 \left( k - 1 \right)}}
```

Where the subscript 0 indicates the *initial conditions in the vessel*. This simplifies the expression to

```math
\frac{dP}{dt} = -\frac{k}{\tau} P \left( \frac{P}{P_0} \right)^{\frac{k-1}{2k}}
```

Which is seperable and can be integrated to give (after some re-arrangement)

```math
\frac{P}{P_0} = \left( 1 + \left( {k-1} \over 2 \right) \frac{t}{\tau} \right)^{\frac{2k}{1-k}}
```

and the de-pressure time is

```math
t = \frac{2\tau}{1-k} \left( 1 - \left( \frac{P_a}{P_0} \right)^{\frac{1-k}{2k}} \right)
```

Another useful thing to determine is the mass flowrate over time, which can be recovered rather easily recalling

```math
w = -\frac{V}{a^2} \frac{dP}{dt} = -\frac{\rho V}{k P} \frac{dP}{dt}
```

and

```math
\frac{dP}{dt} = -\frac{k}{\tau} P \left( \frac{P}{P_0} \right)^{\frac{k-1}{2k}}
```

we get

```math
w = \frac{\rho V}{\tau} \left( \frac{P}{P_0} \right)^{\frac{k-1}{2k}} = \frac{\rho_0 V}{\tau} \left( \frac{\rho}{\rho_0} \right) \left( \frac{P}{P_0} \right)^{\frac{k-1}{2k}}
```

By recalling the definition of $\tau$ this simplifies to

```math
\frac{w}{w_0} = \left( \frac{P}{P_0} \right)^{{k+1} \over {2k}} = \left( 1 + \left({k-1} \over 2 \right) \frac{t}{\tau} \right)^{\frac{1+k}{1-k}}
```

This final model, for mass flow, is the model most often given in process safety references for blowdown rates. This makes some sense as early in a blowdown the observed pressure curve tend to approximate the adiabatic curve. However (foreshadowing) the isothermal curve leads to higher predicted vessel pressures, and generally higher masss flowrates, which might be more conservative depending on the context.


[^2]: Any physical chemistry textbook, such as [Laidler, Meiser, and Sanctuary](#laidler-2003), *Physical Chemistry*, 79 - 81

[^3]: [Tilton](#tilton-2008), "Fluid and Particle Dynamics," 6-22 - 6-23.

"""

# ╔═╡ 34c796e3-5b2f-4869-b446-614a3f2d4645
md"""
### In the Literature

The adiabatic model is the only simple model given in Lees[^4], with the recommendation to use software such as BLOWDOWN to handle more complex, multiphase, mixtures and heat transfer problems. This is also what my older copy of Perry's gives[^5], albeit with a typo.

!!! note
	Perry's gives the following
	```math
	\frac{w}{w_0} = \left( 1 + \left(\mathbf{k + 1} \over 2 \right) \frac{t}{\tau} \right)^{\frac{1+k}{1-k}}
	```
	Note the sign change, it should be **k-1** not **k+1**, given typical values of *k*~1.4 this actually a huge difference.

Perry's only gives the mass flow, so if you wanted the pressure (and the gas density and temperature) you would need to find some other reference. Or do it yourself, it does sketch out how the equation is derived, if you have the spare time to sit down and integrate.

[^4]: [Lees](#lees-1996), *Loss Prevention in the Process Industries*, 15/44
[^5]: [Crowl *et al.*](#crowl-2008), "Process Safety," 23-57

"""

# ╔═╡ da150d10-3374-43c5-b2c5-b9258c4d0d24
md"""
### The Complete ODE

There are two obvious limitations to this model: it relies on the gas being well approximated by an ideal gas and that the flow out of the vessel is always choked. The first issue I am not going to deal with right now, the second one I think can be easily dealt with by slightly modifying the governing equations.

```math
\frac{dP}{dt} = -\frac{c_D A}{V} a^2 G
```

We can solve this numerically given
```math
\rho = \rho_0 \left(\frac{P}{P_0}\right)^{\frac{1}{k}}
```

"""

# ╔═╡ a28fedec-75d4-4260-88da-3a610ed953d2
function isentropic_mass_flow(P, ρ; k=1.4, Pₐ=101325)
	η = max( Pₐ/P, (2/(k+1))^(k/(k-1)))
	G² = ρ*P*(2k/(k-1))*( η^(2/k) - η^((k+1)/k) )
	G = G² > 0 ? √(G²) : 0
	return G
end

# ╔═╡ 4b997747-ee92-413c-a222-26f84fac6a7f
function speed_of_sound(P, ρ; k=1.4)
	a = √(k*P/ρ)
	return a
end

# ╔═╡ e13e3ccc-8a3c-44fe-bdcd-8e6f82a85d4e
function adiabatic_vessel(P, params, t)
	c, A, V, k, ρ₀, P₀, Pₐ = params
	ρ = ρ₀*(P/P₀)^(1/k)
	a² = speed_of_sound(P, ρ; k=k)^2
	G = isentropic_mass_flow(P, ρ; k=k, Pₐ=Pₐ)
	return-c*A*a²*G/V
end

# ╔═╡ da8cb7ac-cd1b-4aa5-bdfe-883bfa102da0
function depressured_callback(P, t, integrator; tol=0.001)
	c, A, V, k, ρ₀, P₀, Pₐ = integrator.p
	return P - (1+tol)*Pₐ
end

# ╔═╡ 3991aef2-cb3a-4f0b-b199-45e23dbe13d0
md"""
### A Motivating Example

Just to have a real system to think about, I used to SCUBA dive when I was a teenage and had a few mis-haps early on, when I was still figuring things out, accidentally opening the tank valve when the regulator yoke was not fully attached. Blasting air all over the place while I scrambled to shut it off. Typical tanks have capacities ranging from 80 cu. ft. to 100 cu. ft., with working pressures of >3000 psi. That's a pretty high pressure for a relatively small tank. How fast could the tank blowdown if I opened the valve fully and just sat back and watched?
"""

# ╔═╡ 8cb41467-1430-4541-9fec-3301b76e688d
# Ambient conditions
begin
	Pₐ = 101.325e3 # Pa
	Tₐ = 288.15    # K
	ρₐ = 1.21      # kg/m³
end;

# ╔═╡ ee4153d6-cf19-4298-ad69-f1b08bbb5845
md"""

I looked around online and a typical tank with a [80 cu. ft. capacity](https://www.divegearexpress.com/library/articles/calculating-scuba-cylinder-capacities) might have a "water volume" (actual internal volume) of 678 cu. in. (11.11L) and a working pressure of 3000 psi (20.68 MPa). I don't actually know the flow area of a tank valve, I couldn't find it easily, so I'm going to guess it's basically a 5 mm diameter tube when fully open, with a discharge coefficient of 0.85 -- all of this could be firmed up better with some real details of the valve. But this is a start.

"""

# ╔═╡ d6a64773-2e49-46f7-8f76-53df6b5475a8
#Vessel conditions
begin
	c = 0.85
	D = 0.005      # m
	A = 0.25*π*D^2 # m²
	V = 0.01111    # m³
	P₀ = 20.68e6   # Pa
	T₀ = Tₐ
	ρ₀ = ρₐ*(P₀/Pₐ) # ideal gas law
	k = 1.4
end;

# ╔═╡ ad7a013c-2779-49b2-9ca9-096ff1521ba0
md"""
I then set up the differential equation and integrate to get the blowdown curve.
"""

# ╔═╡ 72263187-a1f8-4453-a5da-feb8dcd9f318
begin
    params = (c, A, V, k, ρ₀, P₀, Pₐ)
    t_span = (0.0, 100.0)
    prob = ODEProblem(adiabatic_vessel, P₀, t_span, params)
    sol = solve(prob, Tsit5(),
	            callback=ContinuousCallback(depressured_callback, terminate!))
end;

# ╔═╡ 4b3b77dd-831f-4727-8b3f-2c892d9e8068
begin
	τ = 1/( (c*A/V)*√(k*P₀/ρ₀)*(2/(k+1))^((k+1)/(2*(k-1))) )
	choked_model(t) = ( 1 + 0.5*(k-1)*(t/τ))^((2k)/(1-k))

	ts1 = range(0.0, 1.1*sol.t[end]; length=100)
	p1 = plot(xlabel="Time, s", ylabel="Vessel pressure, Pa")
	plot!(p1, ts1, choked_model.(ts1).*P₀, label="Fully choked model")
	plot!(p1, ts1, sol.(ts1; idxs = 1), label="ODE")
	plot!(p1, [ts1[1], ts1[end]], [Pₐ, Pₐ], label="Ambient Pressure")
	lens!(p1, [18, 22], [0.75,1.25].*1e5, inset = (1, bbox(0.5, 0.3, 0.4, 0.4)))
end

# ╔═╡ 99b62b74-db91-47d3-9be0-bdb45d5b3ae2
md"""

So this is saying it blows down pretty fast, in less than 30s. Probably my guess for the valve area is too large. I did just make it up.

Regarding the models themselves, the adiabatic choked model is a very good approximation to the full ODE until the last few fractions of a second, at which point the models diverge. This likely to be true for any high pressure blowdowns, where the vessel pressure starts well above ~2atm, as in that case the majority of the blowdown will be entirely in the choked flow regime.

To play around with this more, I am first going to detour into creating some helper functions and I think this is a natural point to create a `struct` to contain the vessel parameters.
"""

# ╔═╡ 0bd14ea4-0976-4a35-aa43-150963b3d8e2
begin

struct PressureVessel{F <: Number}
	c::F
	A::F
	V::F
	k::F
	ρ₀::F
	P₀::F
	Pₐ::F
	τ::F
end

PressureVessel(c, A, V, k, ρ₀, P₀, Pₐ, τ) = 
    PressureVessel(promote(c, A, V, k, ρ₀, P₀, Pₐ, τ)...)

function PressureVessel(c, A, V, k, ρ₀, P₀, Pₐ)
	τ = 1/( (c*A/V)*√(k*P₀/ρ₀)*(2/(k+1))^((k+1)/(2*(k-1))) )
	return PressureVessel(c, A, V, k, ρ₀, P₀, Pₐ, τ)
end

end;

# ╔═╡ 528a71da-49de-44c9-aeb7-276dfc1e927a
md"""
Where I have added a helper function to ensure all numbers are of the same type, and calculate the value of τ when the `PressureVessel` type is constructed.

Re-creating the results from above, I start with a definition of the vessel
"""

# ╔═╡ 3738a3a0-19b6-4e69-bd34-40cb332bd00e
vessel = PressureVessel(c, A, V, k, ρ₀, P₀, Pₐ);

# ╔═╡ 2efdeb40-2ab7-4b34-9650-66ca8d167846
md"""

I would like to create some generic functions for the blowdown properties I am interested in: pressure and mass flowrate as functions of time and total blowdown time. To accomodate this I define another type to contain the `VesselBlowdown` solution.

"""

# ╔═╡ abd51f03-9fa7-4611-a6c4-d0b41441fae8
abstract type Blowdown end

# ╔═╡ 4bf39b76-9958-49f6-b436-7e02c27dfa20
struct AdiabaticBlowdown{S} <: Blowdown
	pv::PressureVessel
	sol::S
end

# ╔═╡ ad20eb31-a268-4d27-bf58-30ef690afcee
md"""
Here I add some functions to make a `Blowdown` object act like an iterator with only a single element. This is absolutely pointless except that I just happen to like being able to generate a vector of results by using the "dot" notation, like so

```julia
my_function.(blowdown, time_vector)
```

where I want it to broadcast over the `time_vector`.

"""

# ╔═╡ 755eb895-8434-4e3c-900a-8fdd7bcf5a00
Base.length(::Blowdown) = 1

# ╔═╡ f7271b10-5e05-4ec7-a2e4-84eefd5ff014
Base.iterate(b::Blowdown, state=1) = state > length(b) ? nothing : (b,state+1)

# ╔═╡ 6db7af5e-4696-422a-8ef6-172ff207b71f
md"""

For the simple choked model this is fairly straight forward.

"""

# ╔═╡ 7c6290ff-08e1-42bb-bc2c-7f747a908368
adiabatic_blowdown_choked(vessel::PressureVessel) = AdiabaticBlowdown(vessel,nothing)

# ╔═╡ 5a06c184-209b-4ad7-9e27-9fe81946682c
function blowdown_pressure(bd::AdiabaticBlowdown, t)
	P₀, k, τ = bd.pv.P₀, bd.pv.k, bd.pv.τ
	return P₀*( 1 + 0.5*(k-1)*(t/τ))^((2*k)/(1-k))
end

# ╔═╡ 90bf3581-ff00-44c7-b09b-9a6b78616ae9
function blowdown_time(bd::AdiabaticBlowdown)
	P₀, Pₐ, k, τ = bd.pv.P₀, bd.pv.Pₐ, bd.pv.k, bd.pv.τ
	return (2τ/(1-k))*(1 - (Pₐ/P₀)^((1-k)/2k))
end

# ╔═╡ 3db31507-2ddb-48dd-acfa-c5b58fe996fd
md"""

For the full model the initial step is to integrate the differential equation. As a first guess, I calculate the blowdown time for a fully choked blowdown and set the outer-bound for the integration to 10× this. The integrator will terminate when the pressure reaches ambient and thus the last time stored will be the actual blowdown time.

"""

# ╔═╡ 85c9a592-a910-4e82-8a7f-2231e7a19e27
function adiabatic_blowdown_full(vessel::PressureVessel; solver=Tsit5())
	# unpack the parameters
	c, A, V, k, ρ₀, P₀, Pₐ = vessel.c, vessel.A, vessel.V, 
							 vessel.k, vessel.ρ₀, vessel.P₀, vessel.Pₐ
	params = (c, A, V, k, ρ₀, P₀, Pₐ)

	# estimate the time span needed to fully blowdown
	τ = vessel.τ
	t_bd = (2τ/(1-k))*(1 - (Pₐ/P₀)^((1-k)/2k))
    t_span = (0.0, 10t_bd)

	# set up the ODEProblem and solve
    prob = ODEProblem(adiabatic_vessel, P₀, t_span, params)
    sol = solve(prob, solver,
	            callback=ContinuousCallback(depressured_callback, terminate!))

	return AdiabaticBlowdown(vessel,sol)
end

# ╔═╡ b7625eee-8577-4ecb-ac48-ae06da90ec74
blowdown_time(bd::AdiabaticBlowdown{<:ODESolution}) = bd.sol.t[end]

# ╔═╡ fdb1088b-0e60-4b88-8b4d-ca79bec26120
md"""

At this point I've built up enough machinery that playing around with all sorts of variations to the original case become quite simple. As an example, I look at the same air tank but pressured to 1.5 atm instead.

"""

# ╔═╡ 9578f87c-4286-40f5-833a-d235f7e842eb
test_vessel = PressureVessel(c, A, V, k, ρ₀, 1.5Pₐ, Pₐ);

# ╔═╡ d68d6fc7-df19-42cd-9e35-4715b349455d
md"""
## The Isothermal Case

The other limiting case worth exploring is the isothermal case, which is equivalent to the vessel having perfectly conductive walls and remaining always at thermal equilibrium with the environment. This is often a good approximation for large vessels where the blowdown rate is small relative ot the thermal mass of the gas in the vessel.

Recalling, for the adiabatic case, we had the following

```math
\frac{dP}{dt} = - \frac{c_D A}{V} \left( \frac{\partial P}{\partial \rho} \right)_S G
```

For the isothermal case the vessel is being de-pressured along an isothermal path (not an isentropic path) and so we substitute the appropriate partial derivative[^6]

```math
\frac{dP}{dt} = - \frac{c_D A}{V} \left( \frac{\partial P}{\partial \rho} \right)_T G
```

### Fully Choked Flow

As before, the blowdown is through an isentropic nozzle and we assume that flow is choked

```math
G = \sqrt{k \rho P} \left( \frac{2}{k+1} \right)^{\frac{k+1}{2 \left( k-1 \right)}} = \rho \sqrt{ {k P} \over \rho} \left( \frac{2}{k+1} \right)^{\frac{k+1}{2 \left( k-1 \right)}} 
```

From [thermodynamics](https://aefarrell.github.io/2024/10/28/relief_valve_sizing/#fn-9) we can write the partial derivative as

```math
\left( \frac{\partial P}{\partial \rho} \right)_T = \frac{a^2}{k} = \frac{P}{\rho}
```

Thus

```math
\frac{dP}{dt} = - \frac{c_D A}{V} \frac{P}{\rho} \rho \sqrt{ {k P} \over \rho} \left( \frac{2}{k+1} \right)^{\frac{k+1}{2 \left( k-1 \right)}}
```

where the densities, $\rho$, can be cancelled and, since the vessel is isothermal (i.e. $\frac{P}{\rho}$ is a constant), the various constants can be collected to give

```math
\frac{dP}{dt} = -\frac{P}{\tau}
```

Where $\tau$ is as defined for the adiabatic case. This can easily be integrated to give

```math
\frac{P}{P_0} = \exp \left( \frac{-t}{\tau} \right)
```

it also follows, from the ideal gas law, that

```math
\frac{\rho}{\rho_0} = \exp \left( \frac{-t}{\tau} \right)
```

and

```math
\frac{w}{w_0} = \exp \left( \frac{-t}{\tau} \right)
```

This can also be re-arranged to give the blowdown time[^7]

```math
t = \tau \log \left( \frac{P_0}{P_a} \right)
```

[^6]: [Botros, Jungowski, and Weiss](#botros-1989), "Models and Methods of Simulating Gas Pipeline Blowdown."
[^7]: *N.B.* the $\log \left( \dots \right)$ is the natural log, this matches the convention used in julia
"""

# ╔═╡ e7c35605-6e1e-416f-9474-0ef6fbbb99f2
md"""
### In the Literature

This is the equation seen most often in references for blowdown time. It is also what is going on under the hood with many online calculators for vessel blowdown times. Though, in my experience, this is not always well documented and a modified form is often presented.

The time constant, $\tau$, can be broken up to look like this

```math
\tau = \frac{V}{c_D A} \sqrt{\frac{M}{M_{air} Z_0 T_0}} \sqrt{\frac{M_{air}}{kR}} \left( 2 \over {k+1} \right)^{\frac{-1}{2} \frac{k+1}{k-1}}
```

Where we have made the substitution $Z_0 R T_0$ for $R T_0$ to account for non-ideal behaviour. If the gas has a value of *k* ~ 1.4, we can write

```math
\tau = \mathrm{const} \frac{V}{c_D A} \sqrt{ \frac{SG}{Z_0 T_0}}
```

Where the constant is calculated entirely from the properties of air. Generally, I have found, few references describe where this constant comes from and in particular that it depends implicitly on a particular value for *k*. It often has unit conversions absorbed into it, for example[^8]

```math
t = 5.5 \frac{V}{c_D A} \sqrt{ \frac{SG}{Z_0 T_0} } \log \left( \frac{P_0}{P_a} \right)
```

with the units
+ Blowdown time, *t*, in seconds
+ Vessel volume, *V*, in cubic feet
+ Valve flow area, *A*, in square inches
+ Initial temperature, $T_0$, in Rankine
+ Initial pressure, $P_0$, in psia
+ Ambient pressure, $P_a$, in psia

Though I have found a few sources that leave the value of the constant a mystery for the user to puzzle out[^9]. I was honestly surprised at the quality of the results when I first googled this and looked it up in [Knovel](https://app.knovel.com/kn). The highest ranked results, at the time, were cryptic to the point of uselessness or included obivious mistakes (several referred to *t* as the "interstitial velocity" with units of cm/s, an obvious mis-print being blindly recopied in several places, including some e-books on Knovel where one would hope the quality control would be better). There are a few places with useful derivations[^10] but I think a good starting point is the [Tank Blowdown Math](http://www.et.byu.edu/~wheeler/Tank_Blowdown_Math.pdf) set of notes. It is pretty straight forward and does not require a lot of prior knowledge of the partial derviatives of various thermodynamic state variables.

I personally would not bother with the models that pre-calculate the constant for you. We no longer live in the age of slide-rules. The blowdown time equation for fully choked flow is well within the capabilities of excel or any competent person with a scientific calculator. I think it is easier to justify and explain, will be a better model for gases where *k* is not 1.4, and allows one to incorporate small levels of non-ideality through the [isentropic expansion factor](https://aefarrell.github.io/2024/10/28/relief_valve_sizing/#the-isentropic-expansion-factor) *n*.


[^8]: [Campbell](#campbell-1992), *Gas Conditioning and Processing,* 29; [VANEC](#vanec-2025), "Pressure Volume-Blowdown Time Calculation,"

[^9]: [Temizel *et al.*](#temizel-2019), *Formulas and Calculations for Petroleum Engineering*, 262; [Engineers Edge](#engineers-2025), "Blowdown Time in Unsteady Gas Flow Calculator and Equation."

[^10]: [Wheeler](#wheeler-2019), "Tank Blowdown Math."; [Botros, Jungowski, and Weiss](#botros-1989), "Models and Methods of Simulating Gas Pipeline Blowdown."; [Botros and Van Hardeveld](#botros-2018), *Pipeline Pumping and Compression Systems*, 447; [Saad](#saad-1985), *Compressible Fluid Flow*, 98-103, to list but a few.
"""

# ╔═╡ 8bd85419-3927-4f05-b687-3bf4bb178e09
md"""
The isothermal fully-choked model can be implemented building on the types already created, by first creating an `IsothermalBlowdown` type and associated blowdown functions
"""

# ╔═╡ 3bc45288-754a-4afd-a739-0e245f54f3c3
struct IsothermalBlowdown{S} <: Blowdown
	pv::PressureVessel
	sol::S
end

# ╔═╡ df02a298-cead-4511-813c-aacb89d0e102
isothermal_blowdown_choked(vessel::PressureVessel) = IsothermalBlowdown(vessel,nothing)

# ╔═╡ 6a1d0299-1ffe-49d3-87c5-5b989c6ec4e5
function blowdown_pressure(bd::IsothermalBlowdown, t)
	P₀, τ = bd.pv.P₀, bd.pv.τ
	return P₀*exp(-t/τ)
end

# ╔═╡ 9ebfad46-d713-4240-8da1-eb61942f05b6
function blowdown_mass_rate(bd::IsothermalBlowdown, t)
	ρ₀, V, τ = bd.pv.ρ₀, bd.pv.V, bd.pv.τ
	m₀ = ρ₀*V
	w₀ = m₀/τ
	return w₀*exp(-t/τ)
end

# ╔═╡ c3b703dd-2efa-47c4-90ab-45ee47cad7fd
blowdown_time(bd::IsothermalBlowdown) = bd.pv.τ*log(bd.pv.P₀/bd.pv.Pₐ)

# ╔═╡ 1283a0cc-f235-4a76-b6e8-8cf0fd06c4f0
md"""

In a similar vein as the adiabatic case, the requirement for fully choked flow can be relaxed and the ODE integrated numerically instead, starting with the system

```math
\frac{dP}{dt} = -\frac{c_D A}{V} \frac{a^2}{k} G
```

We can solve this numerically given that, for an isothermal system, the density is given by
```math
\rho = \rho_0 \left(\frac{P}{P_0}\right)
```

"""

# ╔═╡ 75db611d-dc36-4f49-a5b8-1d1daad18589
function isothermal_vessel(P, params, t)
	c, A, V, k, ρ₀, P₀, Pₐ = params
	ρ = ρ₀*(P/P₀)
	a² = speed_of_sound(P, ρ; k=k)^2
	G = isentropic_mass_flow(P, ρ; k=k, Pₐ=Pₐ)
	return-c*A*a²*G/(k*V)
end

# ╔═╡ 0c40d734-a7dc-4efe-955d-cab30e14573b
function isothermal_blowdown_full(vessel::PressureVessel; solver=Tsit5())
	# unpack the parameters
	c, A, V, k, ρ₀, P₀, Pₐ = vessel.c, vessel.A, vessel.V, 
							 vessel.k, vessel.ρ₀, vessel.P₀, vessel.Pₐ
	params = (c, A, V, k, ρ₀, P₀, Pₐ)

	# estimate the time span needed to fully blowdown
	τ = vessel.τ
	t_bd = τ*log(P₀/Pₐ)
    t_span = (0.0, 10t_bd)

	# set up the ODEProblem and solve
    prob = ODEProblem(isothermal_vessel, P₀, t_span, params)
    sol = solve(prob, solver,
	            callback=ContinuousCallback(depressured_callback, terminate!))

	return IsothermalBlowdown(vessel,sol)
end

# ╔═╡ 0e8afb24-8f3e-47fd-babb-584da99d6ebc
blowdown_time(bd::IsothermalBlowdown{<:ODESolution}) = bd.sol.t[end]

# ╔═╡ 8a1976b0-eb42-45c1-9ee1-8fd255bb1d0f
function blowdown_pressure(bd::AdiabaticBlowdown{<:ODESolution}, t)
	if t < blowdown_time(bd)
		return bd.sol(t)
	else
		return bd.sol.u[end]
	end
end

# ╔═╡ a88b827b-6378-4b45-9df1-ded740a9c95c
function blowdown_pressure(bd::IsothermalBlowdown{<:ODESolution}, t)
	if t < blowdown_time(bd)
		return bd.sol(t)
	else
		return bd.sol.u[end]
	end
end

# ╔═╡ f6ce5c05-ccc0-43e2-b71d-53ad74c8df8c
function blowdown_mass_rate(bd::AdiabaticBlowdown, t)
	ρ₀, V, P₀, k, τ = bd.pv.ρ₀, bd.pv.V, bd.pv.P₀, bd.pv.k, bd.pv.τ
	m₀ = ρ₀*V
	w₀ = m₀/τ
	P = blowdown_pressure(bd, t)
	return w₀*(P/P₀)^((k+1)/(2k))
end

# ╔═╡ cc74796c-45c1-49d7-9700-bc01b9c84c38
function blowdown_mass_rate(bd::AdiabaticBlowdown{<:ODESolution}, t)
	if t < blowdown_time(bd)
		# unpack the parameters
		c, A, k, ρ₀, P₀, Pₐ = bd.pv.c, bd.pv.A, bd.pv.k, bd.pv.ρ₀, 
						      bd.pv.P₀, bd.pv.Pₐ
		
		# calculate w = c*A*G
		P = blowdown_pressure(bd, t)
		ρ = ρ₀*(P/P₀)^(1/k)
		G = isentropic_mass_flow(P, ρ; k=k, Pₐ=Pₐ)
	
		return c*A*G
	else
		return 0.0
	end
end

# ╔═╡ 0e2ab020-ad35-4864-a1c3-271d6121a0c9
let	bdc=adiabatic_blowdown_choked(vessel), bdf=adiabatic_blowdown_full(vessel)
	
	ts = range(0, 1.1blowdown_time(bdf); length=100)
	plt = plot(xlabel="Time, s", ylabel="Vessel pressure, Pa")
	plot!(plt, ts, blowdown_pressure.(bdc, ts), label="Fully choked model", color=:blue)
	plot!(plt, [blowdown_time(bdc)], [Pₐ], seriestype=:scatter, color=:blue, label=nothing)
	plot!(plt, ts, blowdown_pressure.(bdf, ts), label="ODE", color=:red)
	plot!(plt, [blowdown_time(bdf)], [1.001*Pₐ], seriestype=:scatter, color=:red, label=nothing)
	plot!(plt, [ts[1], ts[end]], [Pₐ, Pₐ], label="Ambient Pressure", color=:green)
	lens!(plt, [18, 22], [0.75,1.25].*1e5, inset = (1, bbox(0.5, 0.3, 0.4, 0.4)))
end

# ╔═╡ 6a7393ee-3044-4f30-93bb-ed3f3fe05179
begin
	bdc = adiabatic_blowdown_choked(test_vessel)
	bdf = adiabatic_blowdown_full(test_vessel)
	
	ts3 = range(0, 1.1blowdown_time(bdf); length=100)
	p3 = plot(xlabel="Time, s", ylabel="Vessel pressure, Pa")
	plot!(p3, ts3, blowdown_pressure.(bdc, ts3), label="Fully choked model", color=:blue)
	plot!(p3, [blowdown_time(bdc)], [Pₐ], seriestype=:scatter, color=:blue, label=nothing)
	plot!(p3, ts3, blowdown_pressure.(bdf, ts3), label="ODE", color=:red)
	plot!(p3, [blowdown_time(bdf)], [1.001*Pₐ], seriestype=:scatter, color=:red, label=nothing)
	plot!(p3, [ts3[1], ts3[end]], [Pₐ, Pₐ], label="Ambient Pressure", color=:green)
end

# ╔═╡ b5bad24c-4001-4031-8e62-89551f1f6bf0
md"""
Now it is clear that the fully choked model model isn't working well, it predicts a blowdown time of $( round(blowdown_time(bdc); digits=2) )s whereas numerically solving the ODE gives an answer of $( round(blowdown_time(bdf); digits=2) )s, a $( round( (blowdown_time(bdf)/blowdown_time(bdc) -1)*100; digits=0) )% greater predicted blowdown.

That said...I'm being a little coy about something: the full ODE predicts that the vessel will *never* blowdown. The pressure will get closer and closer to ambient but never get there. This is because *G*, for non-choked flow, asymptotically approaches zero as the vessel pressure approaches ambient pressure. How you define *blowdown time* is really a function of *how close to ambient* is close enough. Even if I set the tolerance in the `depressured_callback` function, which terminates the integration once the integrator is within tolerance of the ambient pressure, to zero it would, in reality, simply terminate at the default numerical precision of `DifferentialEquations.jl`. In this case I've said "within 0.1% of ambient is close enough," but that's totally arbitrary.
"""

# ╔═╡ 03f26143-9587-4fa6-a101-4b2aff69c31f
function blowdown_mass_rate(bd::IsothermalBlowdown{<:ODESolution}, t)
	if t < blowdown_time(bd)
		# unpack the parameters
		c, A, k, ρ₀, P₀, Pₐ = bd.pv.c, bd.pv.A, bd.pv.k, bd.pv.ρ₀, 
						      bd.pv.P₀, bd.pv.Pₐ
		
		# calculate w = c*A*G
		P = blowdown_pressure(bd, t)
		ρ = ρ₀*(P/P₀)
		G = isentropic_mass_flow(P, ρ; k=k, Pₐ=Pₐ)
	
		return c*A*G
	else
		return 0.0
	end
end

# ╔═╡ c574aee8-a963-4da5-b69f-29c3081bffb2
let	bdc=adiabatic_blowdown_choked(vessel), bdf=adiabatic_blowdown_full(vessel)
	
	ts = range(0, 1.1blowdown_time(bdf); length=100)
	plt = plot(xlabel="Time, s", ylabel="Mass flow rate, kg/s")
	plot!(plt, ts, blowdown_mass_rate.(bdc, ts), label="Fully choked model", color=:blue)
	plot!(plt, ts, blowdown_mass_rate.(bdf, ts), label="ODE", color=:red)
	plot!(plt, [blowdown_time(bdf)], [0.0], seriestype=:scatter, color=:red, label=nothing)
	lens!(plt, [18, 22], [0,0.01], inset = (1, bbox(0.5, 0.3, 0.4, 0.4)))
end

# ╔═╡ f1f17a1f-39de-4400-81fa-82608a80c503
let	bdc=isothermal_blowdown_choked(vessel), bdf=isothermal_blowdown_full(vessel)
	
	ts = range(0, 1.1blowdown_time(bdf); length=100)
	plt = plot(xlabel="Time, s", ylabel="Vessel pressure, Pa")
	plot!(plt, ts, blowdown_pressure.(bdc, ts), label="Fully choked model", color=:blue)
	plot!(plt, [blowdown_time(bdc)], [Pₐ], seriestype=:scatter, color=:blue, label=nothing)
	plot!(plt, ts, blowdown_pressure.(bdf, ts), label="ODE", color=:red)
	plot!(plt, [blowdown_time(bdf)], [1.001*Pₐ], seriestype=:scatter, color=:red, label=nothing)
	plot!(plt, [ts[1], ts[end]], [Pₐ, Pₐ], label="Ambient Pressure", color=:green)
	lens!(plt, [17, 19], [0.95,1.2].*Pₐ, inset = (1, bbox(0.5, 0.3, 0.4, 0.4)))
end

# ╔═╡ a70d5c15-b4c5-4767-ace7-b1acdc9b515c
let	bdc=isothermal_blowdown_choked(vessel), bdf=isothermal_blowdown_full(vessel)
	
	ts = range(0, 1.1blowdown_time(bdf); length=100)
	plt = plot(xlabel="Time, s", ylabel="Mass flow rate, kg/s")
	plot!(plt, ts, blowdown_mass_rate.(bdc, ts), label="Fully choked model", color=:blue)
	plot!(plt, ts, blowdown_mass_rate.(bdf, ts), label="ODE", color=:red)
	plot!(plt, [blowdown_time(bdf)], [0.0], seriestype=:scatter, color=:red, label=nothing)
	lens!(plt, [17, 19], [0,0.005], inset = (1, bbox(0.5, 0.3, 0.4, 0.4)))
end

# ╔═╡ a6c0405d-6331-4698-90f1-e7b7fa81a9fd
md"""
It is a similar story to the adiabatic case: for systems with a high initial pressure, the flow out of the valve is fully choked for almost the entire blowdown. It is only in the final fraction of a second that the full ODE system deviates from the model that assumes flow is choked all the time.

In most practical situations, the difference would likely be swamped by two much greater problems with these models:
1. the gases are assumed be ideal with constant *k*
2. the vessel is perfectly isothermal (or adiabatic)

Both of these assumptions will have a much greater impact on how well the model fits observed blowdowns than the slight deviation at the end of the blowdown due to non-choked flow.
"""

# ╔═╡ 59f921ad-a620-4497-aa07-9530023aa616
md"""
## Comparing Blowdown Models
"""

# ╔═╡ d09d9cbd-d215-430d-9e5c-5a3b20d711f2
md"""

I think it might be simpler to visualize when the choked flow blowdown models will fall down by looking at the high pressure blowdown, the original example, versus the low pressure blowdown in dimensionless form. In this form, the choked flow blowdown curves (both adiabatic and isothermal) only depend on *k*. They are in fact the exact same curve. All that has changed is where along the curve the blowdown terminates.

"""

# ╔═╡ b18c9907-f0ee-452a-bd51-c6dc3ebf1401
let	abdc=adiabatic_blowdown_choked(vessel), ibdc=isothermal_blowdown_choked(vessel)
	ts = range(0, 1.1blowdown_time(abdc); length=100)
	τs = ts./vessel.τ
	plt = plot(xlabel="t/τ", ylabel="P/P₀", xlims=(0,6), ylims=(0,1), title="Adiabatic and Isothermal Blowdown")
	annotate!(plt, 1, 1, Plots.text("high pressure blowdown", 10, :left))
	plot!(plt, τs, blowdown_pressure.(ibdc, ts)./P₀, label="Isothermal", color=:blue)
	plot!(plt, [blowdown_time(ibdc)]./vessel.τ, [Pₐ]./P₀, seriestype=:scatter, color=:blue, label=nothing)
	plot!(plt, τs, blowdown_pressure.(abdc, ts)./P₀, label="Adiabatic", color=:orange)
	plot!(plt, [blowdown_time(abdc)]./vessel.τ, [Pₐ]./P₀, seriestype=:scatter, color=:orange, label=nothing)
	plot!(plt, [τs[1], τs[end]], [Pₐ, Pₐ]./P₀, label="Ambient Pressure", color=:green)
	lens!(plt, [5, 6], [0,0.01], inset = (1, bbox(0.5, 0.3, 0.4, 0.4)))
end

# ╔═╡ 55f740e7-bbbd-420f-85bc-2cd52b7509f1
md"""
In the high pressure case the blowdown terminates much closer to $\frac{P}{P_0}=0$ and most of the curve is fully choked.
"""

# ╔═╡ f24f515d-4694-476a-b70c-d138ad065123
let	v = test_vessel
	abdc=adiabatic_blowdown_choked(v)
	abdf=adiabatic_blowdown_full(v)
	ibdc=isothermal_blowdown_choked(v)
	ibdf=isothermal_blowdown_full(v)
	τs = range(0.0, 6.0; length=100)
	ts = τs.*v.τ
	plt = plot(xlabel="t/τ", ylabel="P/P₀", xlims=(0,6), ylims=(0,1), title="Adiabatic and Isothermal Blowdown")
	annotate!(plt, 1, 1, Plots.text("low pressure blowdown", 10, :left))
	plot!(plt, τs, blowdown_pressure.(ibdf, ts)./v.P₀, label="Isothermal", color=:blue)
	plot!(plt, τs, blowdown_pressure.(ibdc, ts)./v.P₀, label="... choked flow", color=:blue, linestyle=:dash)
	plot!(plt, τs, blowdown_pressure.(abdf, ts)./v.P₀, label="Adiabatic", color=:orange)
	plot!(plt, τs, blowdown_pressure.(abdc, ts)./v.P₀, label="... choked flow", color=:orange, linestyle=:dash)
	plot!(plt, [τs[1], τs[end]], [Pₐ, Pₐ]./v.P₀, label="Ambient Pressure", color=:green)
	# lens!(plt, [6.5, 7.5], [0.75,1.25], inset = (1, bbox(0.5, 0.3, 0.4, 0.4)))
end

# ╔═╡ 83b4cc12-0c6a-46de-98f1-c3f619662957
md"""
In the low pressure case the blowdown terminates at a much steeper part of the blowdown curve and the departure for non-choking flow is much more apparent.
"""

# ╔═╡ 1a4ed3cd-4609-43e5-bb98-098f99584d54
md"""
It is not immediately clear to me why the adiabatic case is all over the standard references for process safety. If what you care about is the pressure sustained within a vessel, the mass flowrate emitted through a blowdown stack or vent, and the duration of the blowdown, it is almost always more conservative to use the isothermal case. The isothermal (fully choked) model is also just easier to calculate, being just $\exp \left( \frac{-t}{\tau} \right)$.

The adiabatic case will give a better sense of how temperature changes within the vessel. I've largely left it out, but adiabatic blowdown does lead to a significant temperature drop and this cryogenic cooling can be a process hazard on its own. The gas exiting, and the vessel walls themselves, will get quite cold. Anyone who has gone camping in more marginal weather and watched a one-pound propane cylinder develop frost on the outside while cooking has seen this effect in action.[^11] But actually calculating the vessel temperature is almost entirely ignored in blowdown calculations for ideal gases, in my experience.

I can see why the isothermal model is seemingly universally used for estimating blowdown times. The two curves do cross at high $\frac{t}{\tau}$ and so it is not *always* the case that the isothermal model is more conservative, but it likely to be conservative or *close enough*. It is also very easy to calculate.

[^11]: This is also why butane cylinders are often not a good idea for early spring camping (in Canada), the cooling effect is strong enough to cause the butane inside to liquify and the stove won't work very well.
"""

# ╔═╡ 14003d09-7ed4-41df-b7b2-11d1953a637c
md"""
## Final Thoughts

I deliberately set up the ODEs such that there is a clear path to implementing a real gas model through an equation of state. All that really needs to be done is create functions for these three steps:
1. the speed of sound
2. the density as a function of pressure, either along an isentropic path (in the adiabatic case) or along an isothermal path
3. the isentropic mass velocity, *G*

Plugging those into the relevant steps in the `adiabatic_vessel` and `isothermal_vessel` functions changes from the ideal gas case to the real gas case. The rest of the code remains the same and operates unchanged.

In this case I think solving the full ODE for the ideal gas case alone is probably not worth the effort for most cases. The error in assuming an ideal gas, or in assuming one of the limiting heat transfer cases, is probably far larger than the error in assuming fully choked flow for all but the few cases that are near atmospheric pressure. If you are going to be estimating the blowdown for a real gas, then that's different. If you are going to the hassle of setting up and solving the ODE, might as well have as few un-necessary assumptions as you can get away with. It really isn't any more *person* effort, at that point, just more *computer* effort, and when the calculations happen in less than a second, how much less than a second is of little practical importance.
"""

# ╔═╡ 1ab74775-4f00-4552-85a9-b8a679244240
md"""
## References
"""

# ╔═╡ d0d2052a-506e-4272-a05f-ad593dbc127e
html"""
<ul>
<li><a name="botros-1989">Botros</a>, Kamal K., W. M. Jungowski, and M. H. Weiss. "Models and Methods of Simulating Gas Pipeline Blowdown." <em>The Canadian Journal of Chemical Engineering</em> 67, no. 4 (August 1989): 529-539. [<a href="https://doi.org/10.1002/cjce.5450670402">doi: 10.1002/cjce.5450670402</a>]</li>
<li><a name="botros-2018">Botros</a>, Kamal K., and Thomas Van Hardeveld. <em>Pipeline Pumping and Compression Systems - A Practical Approach</em>. 3rd ed. New York: ASME Press. 2018</li>
<li><a name="campbell-1992">Campbell</a>, John M. <em>Gas Conditioning and Processing</em>. vol 2. Tulsa, OK: John M. Campbell & Co. 1992.</li>
<li><a name="crowl-2008">Crowl</a>, Daniel A., Lawrence G. Britton, Walter L. Frank, Stanley Grossel, Dennis Hendershot, W. G. High, Robert W. Johnson <em>et al</em>. "Process Safety." in Green, <em>Perry’s Chemical Engineers’ Handbook</em>.</li>
<li><a name="engineers-2025">Engineers Edge</a>. "Blowdown Time in Unsteady Gas Flow Calculator and Equation." accessed January 13, 2025 [ <a href="https://www.engineersedge.com/calculators/blowdown_time_in_unsteady_gas_16011.htm">url</a> ]</li>
<li><a name="green-2008">Green</a>, Don W., ed. <em>Perry’s Chemical Engineers’ Handbook</em>. 8th ed. New York: McGraw Hill, 2008.</li>
<li><a name="laidler-2003">Laidler</a>, Keith J., John H. Meiser, and Bryan C. Sanctuary. <em>Physical Chemistry</em>. 4th ed. Boston: Houghton Mifflin Co., 2003.</a></li>
<li><a name="lees-1996">Lees</a>, Frank P. <em>Loss Prevention in the Process Industries</em>. 2nd ed. Oxford: Butterworth-Heinemann, 1996.</li>
<li><a name="saad-1985">Saad</a>, Michel A. <em>Compressible Fluid Flow</em>. Englewood Cliffs, NJ: Prentiss-Hall. 1985</li>
<li><a name="temizel-2019">Temizel</a>, Cenk, Tayfun Tuna, Mehmet Melik Oskay, and Luigi Saputelli. <em>Formulas and Calculations for Petroleum Engineering</em>. Cambridge, MA: Gulf Professional Publishing. 2019.</li>
<li><a name="tilton-2008">Tilton</a>, James N. "Fluid and Particle Dynamics." in Green, <em>Perry’s Chemical Engineers’ Handbook</em>.</li>
<li><a name="vanec-2025">VANEC</a>. "Pressure Volume-Blowdown Time Calculation." accessed January 11, 2025 [ <a href="https://www.vanec.com/pressurized-volume-blowdown-time-calculation.html">url</a> ]</li>
<li><a name="wheeler-2019">Wheeler</a>, Dean R. "Tank Blowdown Math." March 13, 2019 [ <a href="http://www.et.byu.edu/~wheeler/Tank_Blowdown_Math.pdf">url</a> ]</li>
</ul>
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
OrdinaryDiffEq = "~6.89.0"
Plots = "~1.40.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "127892c3d8f685ca4b304c074f6873399b8c5d9b"

[[deps.ADTypes]]
git-tree-sha1 = "eea5d80188827b35333801ef97a40c2ed653b081"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.9.0"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown"]
git-tree-sha1 = "96bed9b1b57cf750cca50c311a197e306816a1cc"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.39"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsDatesExt = "Dates"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsTestExt = "Test"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "d80af0733c99ea80575f612813fa6aa71022d33a"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.0"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3640d077b6dafd64ceb8fd5c1ec76f7ca53bcf76"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.16.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "0dd7edaff278e346eb0ca07a7e75c9438408a3ce"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.10.3"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "TruncatedStacktraces"]
git-tree-sha1 = "f8eefbb7e910f59087c4bb09ce670f235758ee4a"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.158.3"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "16611777adf4a818f18c33210895947814180964"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.6.16"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = "Enzyme"
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = "ForwardDiff"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = "PolyesterForwardDiff"
    DifferentiationInterfaceReverseDiffExt = "ReverseDiff"
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
git-tree-sha1 = "9c3a42611e525352e9ad5e4134ddca5c692ff209"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.4"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e51db81749b0777b2147fbe7b783ee79045b8e99"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.4+3"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "8e18940a5ba7f4ddb41fe2b79b6acaac50880a86"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.26.1"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Expronicon]]
deps = ["MLStyle", "Pkg", "TOML"]
git-tree-sha1 = "fc3951d4d398b5515f91d7fe5d45fc31dccb3c9b"
uuid = "6b7a57c9-7cc1-4fdf-b7f5-e857abae3636"
version = "0.8.5"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "cbf5edddb61a43669710cbc2241bc08b36d9e660"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.4"

[[deps.FastPower]]
git-tree-sha1 = "46aee43f62bc2bc06a74e2d668ffeea0a2689c93"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.1.0"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "b10bdafd1647f57ace3885143936749d61638c3b"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.26.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "af49a0851f8113fcfae2ef5027c6d49d0acec39b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.4"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1dc470db8b1131cfc7fb4c115de89fe391b9e780"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.12.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "07649c499349dad9f08dde4243a4c597064663e9"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.6.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "71dd823bf0a5aaa8264e1143b96a0baa631f4ce0"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.7"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "360f6039babd6e4d6364eff0d4fc9120834a2d9a"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.2.1"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "61dfdba58e585066d8bce214c5a51eaa0539f269"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "97d502765cc5cf3a722120f50da03c2474efce04"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.4"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "e4c3be53733db1051cc15ecf573b1042b3a712a1"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.3.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "ddb60aabc5efb0670c7ddde5fedf7ab44520a3e6"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.36.1"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"
    LinearSolveRecursiveArrayToolsExt = "RecursiveArrayTools"

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "8084c25a250e00ae427a379a5b607e7aed96a2dd"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.171"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "LazyArrays", "LineSearch", "LineSearches", "LinearAlgebra", "LinearSolve", "MaybeInplace", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "4d8944f32db2b07a2bdf8477e878bcb9c9ea2308"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "3.15.1"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = "NLsolve"
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+3"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "cd892f12371c287dc50d6ad3af075b088b6f2d48"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.89.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "8e3c5978d0531a961f70d2f2730d1d16ed3bbd12"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.1.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "b4498d40bf35da0b6d22652ff2e9d8820590b3c6"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.1.2"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "TruncatedStacktraces"]
git-tree-sha1 = "1175717a62ab21736a8f5d0d2531d2a6ad3b9e74"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.9.0"
weakdeps = ["EnzymeCore"]

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"

[[deps.OrdinaryDiffEqDefault]]
deps = ["DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport"]
git-tree-sha1 = "c8223e487d58bef28a3535b33ddf8ffdb44f46fb"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.1.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "e63ec633b1efa99e3caa2e26a01faaa88ba6cef9"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.1.0"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "TruncatedStacktraces"]
git-tree-sha1 = "4dbce3f9e6974567082ce5176e21aab0224a69e9"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.1.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "f63938b8e9e5d3a05815defb3ebdbdcf61ec0a74"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.1.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "0260b1657069608f2a2f26f58f65019d8a9070ab"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.2.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLOperators"]
git-tree-sha1 = "795221c662698851328cb7787965ab4a180d9468"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.1.1"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "a7cc74d3433db98e59dc3d58bc28174c6c290adf"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.1.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "925a91583d1ab84f1f0fea121be1abf1179c5926"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.1.1"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "103e017ff186ac39d731904045781c9bacfca2b0"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.1.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport"]
git-tree-sha1 = "9f8f52aad2399d7714b400ff9d203254b0a89c4a"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.1.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "0f81a77ede3da0dc714ea61e81c76b25db4ab87a"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.1.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "d4bb32e09d6b68ce2eb45fb81001eab46f60717a"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.2.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "590561f3af623d5485d070b4d7044f8854535f5a"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.2.1"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "a2a4119f3e35f7982f78e17beea7b12485d179e9"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.2.1"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "ef44754f10e0dfb9bb55ded382afed44cd94ab57"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.1.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "StaticArrays"]
git-tree-sha1 = "a8b7f8107c477e07c6a6c00d1d66cac68b801bbc"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.1.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport"]
git-tree-sha1 = "da525d277962a1b76102c79f30cb0c31e13fe5b9"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.1.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "332f9d17d0229218f66a73492162267359ba85e9"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.1.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "41c09d9c20877546490f907d8dffdd52690dd65f"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.1.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static"]
git-tree-sha1 = "96b47cdd12cb4ce8f70d701b49f855271a462bd4"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.2.0"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "f6683803a58de600ab7a26d2f49411c9923e9721"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.1.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "StaticArrays"]
git-tree-sha1 = "7dbe4ac56f930df5e9abd003cedb54e25cbbea86"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.2.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "348fd6def9a88518715425025eadd58517017325"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.1.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "1b0d894c880e25f7d0b022d7257638cf8ce5b311"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.1.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport"]
git-tree-sha1 = "4e8b8c8b81df3df17e2eb4603115db3b30a88235"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.1.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "96552f7d4619fabab4038a29ed37dd55e9eb513a"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.1.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "81d7841e73e385b9925d5c8e4427f2adcdda55db"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.1.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "45470145863035bb124ca51b320ed35d071cc6c2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.8"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6d38fea02d983051776a856b7df75b30cf9a3c1f"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.16"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "6c62ce45f268f3f958821a1e5192cf91c75ae89c"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.24"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "b034171b93aebc81b3e1890a036d13a9c4a9e3e0"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.27.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "6db1a75507051bc18bfa131fbc7c3f169cc4b2f6"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.23"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "26fea1911818cd480400f1a2b7f6b32c3cc3836a"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.56.4"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "991d2a8900e687e2c693d587daa739c8fda01366"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.0"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "ef388ca9e4921ec5614ce714f8aa59a5cd33d867"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.11"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "25514a6f200219cd1073e4ff23a6324e4a7efe64"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DiffResults", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "44021f3efc023be3871195d8ad98b865001a2fa1"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "1.12.3"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"
    SimpleNonlinearSolveZygoteExt = "Zygote"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "6914df6005bab9940e2a96879a97a43e1fb1ce78"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "0.6.8"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerDataInterpolationsExt = "DataInterpolations"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    DataInterpolations = "82cc6244-b520-54b8-b5a6-8a565e85f1d0"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "b906758c107b049b6b71599b9f928d9b14e5554a"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.23.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsPolyesterExt = "Polyester"
    SparseDiffToolsPolyesterForwardDiffExt = "PolyesterForwardDiff"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DataStructures", "DocStringExtensions", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "f37f046636f8dc353a39279abfefe296db212171"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.8"
weakdeps = ["Colors"]

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsColorsExt = "Colors"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "87d51a3ee9a4b0d2fe054bdd3fc2436258db2603"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.1.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "47091a0340a675c738b1304b58161f3b0839d454"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.10"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f35f6ab602df8413a50c4a25ca14de821e8605fb"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.7"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "4bc96df5d71515b1cb86dd626915f06f4c0d46f5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.33"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3a6f063d690135f5c1ba351412c82bae4d1402bf"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.25"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "be986ad9dac14888ba338c2554dcfec6939e1393"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.2.1"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "e7f5b81c65eb858bed630fe006837b935518aca5"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.70"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+3"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+2"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─b18f8cd2-d697-11ef-11a8-07514c0e9930
# ╟─4464df47-970c-4fb4-aa46-aad7d2d34cc1
# ╟─0cce5b25-f6a0-4df6-8b58-89d68c5b3abc
# ╟─34c796e3-5b2f-4869-b446-614a3f2d4645
# ╟─da150d10-3374-43c5-b2c5-b9258c4d0d24
# ╠═a28fedec-75d4-4260-88da-3a610ed953d2
# ╠═4b997747-ee92-413c-a222-26f84fac6a7f
# ╠═e13e3ccc-8a3c-44fe-bdcd-8e6f82a85d4e
# ╠═da8cb7ac-cd1b-4aa5-bdfe-883bfa102da0
# ╟─3991aef2-cb3a-4f0b-b199-45e23dbe13d0
# ╠═8cb41467-1430-4541-9fec-3301b76e688d
# ╟─ee4153d6-cf19-4298-ad69-f1b08bbb5845
# ╟─d6a64773-2e49-46f7-8f76-53df6b5475a8
# ╟─ad7a013c-2779-49b2-9ca9-096ff1521ba0
# ╠═320df4b7-f6a8-4d84-a045-f1efc055e417
# ╠═72263187-a1f8-4453-a5da-feb8dcd9f318
# ╟─4b3b77dd-831f-4727-8b3f-2c892d9e8068
# ╟─99b62b74-db91-47d3-9be0-bdb45d5b3ae2
# ╠═0bd14ea4-0976-4a35-aa43-150963b3d8e2
# ╟─528a71da-49de-44c9-aeb7-276dfc1e927a
# ╠═3738a3a0-19b6-4e69-bd34-40cb332bd00e
# ╟─2efdeb40-2ab7-4b34-9650-66ca8d167846
# ╠═abd51f03-9fa7-4611-a6c4-d0b41441fae8
# ╠═4bf39b76-9958-49f6-b436-7e02c27dfa20
# ╟─ad20eb31-a268-4d27-bf58-30ef690afcee
# ╠═755eb895-8434-4e3c-900a-8fdd7bcf5a00
# ╠═f7271b10-5e05-4ec7-a2e4-84eefd5ff014
# ╟─6db7af5e-4696-422a-8ef6-172ff207b71f
# ╠═7c6290ff-08e1-42bb-bc2c-7f747a908368
# ╠═5a06c184-209b-4ad7-9e27-9fe81946682c
# ╠═f6ce5c05-ccc0-43e2-b71d-53ad74c8df8c
# ╠═90bf3581-ff00-44c7-b09b-9a6b78616ae9
# ╟─3db31507-2ddb-48dd-acfa-c5b58fe996fd
# ╠═85c9a592-a910-4e82-8a7f-2231e7a19e27
# ╠═8a1976b0-eb42-45c1-9ee1-8fd255bb1d0f
# ╠═cc74796c-45c1-49d7-9700-bc01b9c84c38
# ╠═b7625eee-8577-4ecb-ac48-ae06da90ec74
# ╟─0e2ab020-ad35-4864-a1c3-271d6121a0c9
# ╟─c574aee8-a963-4da5-b69f-29c3081bffb2
# ╟─fdb1088b-0e60-4b88-8b4d-ca79bec26120
# ╠═9578f87c-4286-40f5-833a-d235f7e842eb
# ╟─6a7393ee-3044-4f30-93bb-ed3f3fe05179
# ╟─b5bad24c-4001-4031-8e62-89551f1f6bf0
# ╟─d68d6fc7-df19-42cd-9e35-4715b349455d
# ╟─e7c35605-6e1e-416f-9474-0ef6fbbb99f2
# ╟─8bd85419-3927-4f05-b687-3bf4bb178e09
# ╠═3bc45288-754a-4afd-a739-0e245f54f3c3
# ╠═df02a298-cead-4511-813c-aacb89d0e102
# ╠═6a1d0299-1ffe-49d3-87c5-5b989c6ec4e5
# ╠═9ebfad46-d713-4240-8da1-eb61942f05b6
# ╠═c3b703dd-2efa-47c4-90ab-45ee47cad7fd
# ╟─1283a0cc-f235-4a76-b6e8-8cf0fd06c4f0
# ╠═75db611d-dc36-4f49-a5b8-1d1daad18589
# ╠═0c40d734-a7dc-4efe-955d-cab30e14573b
# ╠═a88b827b-6378-4b45-9df1-ded740a9c95c
# ╠═03f26143-9587-4fa6-a101-4b2aff69c31f
# ╠═0e8afb24-8f3e-47fd-babb-584da99d6ebc
# ╟─f1f17a1f-39de-4400-81fa-82608a80c503
# ╟─a70d5c15-b4c5-4767-ace7-b1acdc9b515c
# ╟─a6c0405d-6331-4698-90f1-e7b7fa81a9fd
# ╟─59f921ad-a620-4497-aa07-9530023aa616
# ╟─d09d9cbd-d215-430d-9e5c-5a3b20d711f2
# ╟─b18c9907-f0ee-452a-bd51-c6dc3ebf1401
# ╟─55f740e7-bbbd-420f-85bc-2cd52b7509f1
# ╟─f24f515d-4694-476a-b70c-d138ad065123
# ╟─83b4cc12-0c6a-46de-98f1-c3f619662957
# ╟─1a4ed3cd-4609-43e5-bb98-098f99584d54
# ╟─14003d09-7ed4-41df-b7b2-11d1953a637c
# ╟─1ab74775-4f00-4552-85a9-b8a679244240
# ╟─d0d2052a-506e-4272-a05f-ad593dbc127e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
