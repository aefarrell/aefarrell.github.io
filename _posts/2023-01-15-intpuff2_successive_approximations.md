---
title: "Integrating a Gaussian puff - mistakes were made"
toc: true
toc_label: "contents"
toc_sticky: true
comments: true
categories:
  - notes
tags:
  - dispersion modeling
tagline: "successive approximations to ... an integrated gaussian puff model"
header:
  overlay_image: /images/gaussian_dispersion_example_files/veeterzy-unsplash-header.jpg
  caption: "Photo by [veeterzy](https://unsplash.com/@veeterzy?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText) on [Unsplash](https://unsplash.com/)"
---


# Integrating a Gaussian puff - mistakes were made

The other day I was working on a project involving Gaussian puff models and I noticed that I had made a significant mistake, a mistake I have made several times without noticing, and one that invalidated [a whole bunch of work I that I had done previously](https://aefarrell.github.io/2022/06/10/integrated_puff/), so I thought this would be a good opportunity to examine my mistake and it's consequences.

## The Gaussian puff model

To re-cap on what a Gaussian puff model even is: for a short duration release (strictly an instantaneous release) of a neutrally buoyant substance at ground-level, the concentration can be modeled as the product of three Gaussian distributions:

$$ c \left(x,y,z,t \right) = \dot{m} \Delta t \cdot g_x(x, t) \cdot g_y(y) \cdot g_z(z) $$

where

$$ g_x(x,t) = {1 \over \sqrt{2\pi} \sigma_x } \exp \left( -\frac{1}{2} \left( x-u t \over \sigma_x \right)^2 \right) \\
g_y(y) = {1 \over \sqrt{2\pi} \sigma_y } \exp \left( -\frac{1}{2} \left( y \over \sigma_y \right)^2 \right) \\
g_z(z) = {2 \over \sqrt{2\pi} \sigma_z } \exp \left( -\frac{1}{2} \left( z \over \sigma_z \right)^2 \right) $$

Where $\dot{m}$ is the mass emission rate, *&Delta;t* the duration of the release, and *u* the ambient windspeed. The coordinates are such that the release point is at the origin, the puff moves in the downwind, *x*, direction while spreading into the crosswind, *y*, and vertical, *z*, directions.

The dispersion parameters, *&sigma;<sub>x</sub>*, *&sigma;<sub>y</sub>*, *&sigma;<sub>z</sub>* are all functions of the downwind distance and the atmospheric stability.


```julia
# class F puff dispersion
# x is in meters
σx(x) = 0.024*x^0.89
σy(x) = σx(x)
σz(x) = 0.05*x^0.61
```




    σz (generic function with 1 method)



## Integrating the puff

What this generates is an instantaneous release of all of the mass in an infinitesimal point that grows as it moves downwind. This isn't terribly realistic for releases of any appreciable duration (all of the mass is released instantly in this model), so a common approach is to break up the release into a sequence of *n* smaller puffs that each capture the mass released over the sub-interval $ { \Delta t \over n } $. Taking the limit as $n \to \infty$ equates to integrating the puff model from *t - &Delta;t* to *t* giving a nice solution in terms of the error function *erf* and ... this is where I made the critical mistake.

The dispersion parameters are functions of the downwind distance, but critically..*to what?* Taken as the downwind distance to the point being calculated, the dispersion parameters are constants (with respect to time) and the problem simplifies to integrating the Gaussian $ g_{x}(x,t) $ with respect to *t*, which is what I had assumed. *However* if the dispersion parameters are actually correlated to the downwind distance of the *cloud center*, which is $x_c = u t$, they are in fact functions of time and this does not work.

This distinction is by no means made obvious in many of the references for chemical hazard analysis. Most are either vague about it or take the dispersion parameters at the downwind distance *of the point being calculated*. My main reference is the CCPS *Guidelines for Consequence Analysis of Chemical Releases* and it does this (page 107-108). As do several workbooks I have seen. However *Lee's Loss Prevention in the Process Industries, 2nd Ed.* (page 15/112) notes that the dispersion parameters for the Pasquill-Gifford puff model (which this is) are given by

$$ \sigma = { C^2 \over 2 } \left( u t \right)^{2-n} $$

where *C* and *n* are some constants from Sutton, and in general the dispersion correlations are functions of travel time with a lot of discussion in the literature of *to what power*. The standard correlations for the dispersion parameters come from Slade [*Meteorology and Atomic Energy*](https://doi.org/10.2172/4492043) which gives some details on how the measurements were actually taken. It certainly seems to me that the downwind distance was to the cloud center, i.e. the experimenters measured the cloud dimensions at the downwind point to which it had traveled. Which makes the travel time and windspeed implicit.

I think it is a reasonable confusion as the dispersion parameters for a continuous release, a Gaussian *plume* model, are indeed functions of the downwind distance to the point being calculated. It is also frequently the case that examples are given for the concentration at the cloud center, in which case the downwind distance at the point being calculated *is* the downwind distance to the cloud center.


## Dispersion nearly-constants

How critical of a mistake is this? For regions far enough from the origin the dispersion parameters do not vary much in the neighborhood of the plume center. This is shown in the plot below where the difference is taken over the interval $ [ x - \sigma_x, x + \sigma_x ] $. At distances further than a few hundred meters the difference is only a few percent. Suggesting that it might not be an unreasonable approximation to assume the dispersion parameters are constants for the purpose of the integral.


    
![svg](/images/intpuff2_successive_approximations_files/output_3_0.svg)
    



## Different approaches to approximation

Another way of approaching this is simply to view it as an *approximation* instead of an error. On the one hand this is a pretty great rhetorical trick: my answer isn't wrong, it's just *differently true*. But it could be the case that this is a useful simplification, just by eye-balling isopleths and looking at limiting behavior [in the previous notebook](https://aefarrell.github.io/2022/06/10/integrated_puff/) it certainly looked reasonable.

To make life easier, going forward, I am going to define a unit-less time
$$ t = { u t^{\prime} \over L } $$

and unit-less distances

$$ x = {x^{\prime} \over L } \\ y = {y^{\prime} \over L } \\ z = {z^{\prime} \over L } $$

where I am abusing notation with the $\prime$ indicates the variable with units, and no $\prime$ indicates it is unitless. A characteristic length, $L$, is introduced to make everything unitless and, due to the dispersion correlations $ L = 1 \mathrm{m} $ is the most convenient.

We can then explore the performance of different approximations to the integrated puff model by only examining the Gaussian distributions -- with no dependence upon $\dot{m}$ or *u*.


```julia
g(ξ,σ) = exp(-0.5*(ξ/σ)^2)/(√(2π)*σ)

gx(x, t) = g((x-t),σx(t))
gy(y, t) = g(y,σy(t))
gz(z, t) = 2*g(z,σz(t))

pf(x,y,z,t; Δt) = gx(x,t)*gy(y,t)*gz(z,t)*Δt
```



### Sum of discrete puffs

The first type of approximation is to divide the release interval into *n* sub-intervals and *n* Gaussian puffs


```julia
function Σpf(x,y,z,t; Δt, n)
    Δt = min(t,Δt)
    δt = Δt/(n-1)
    _sum = 0
    for i in 0:(n-1)
        t′ = t-i*δt
        pf_i = t′>0 ? gx(x,t′)*gy(y,t′)*gz(z,t′)*δt : 0
        _sum += pf_i
    end
    return _sum
end
```



### Integrating assuming constant &sigma;s

The next type of approximation is the one I made in [the previous post](https://aefarrell.github.io/2022/06/10/integrated_puff/#integrated-puffs) wherein $g_x(x,t)$ is integrated with respect to time, treating the &sigma;s as constants.

There is a little sleight of hand as I include the downwind distance dependence of the &sigma;s after the integration (they aren't *actually* constants)


```julia
using SpecialFunctions: erf

function ∫gx(x,t,Δt)
    Δt = min(t,Δt)
    a  = (x-(t-Δt))/(√2*σx(t-Δt))
    b  = (x-t)/(√2*σx(t))
    return erf(b,a)/2
end

∫pf_approx(x,y,z,t; Δt) = ∫gx(x,t,Δt)*gy(y,x)*gz(z,x)
```



### Numerically integrating the full model

Finally, I take advantage of the `QuadGK` package to numerically integrate the Gaussian puff model, including the time dependence of the dispersion parameters.


```julia
using QuadGK: quadgk

function ∫pf(x,y,z,t; Δt)
    Δt = min(t,Δt)
    integral, err = quadgk(τ -> gx(x,τ)*gy(y,τ)*gz(z,τ), t-Δt, t)
    return integral
end
```



## Comparing performance

### Model error

To give a sense of how these successive approximations work, lets examine a series of slices through the cloud. The first is at a constant *x* on the center-line of the release, looking at how the concentration changes with time.

Just by eye-ball the the approximate integral is very close to the numerical exact(ish) integral, as is a large enough number of puffs. Importantly, I think, the approximate integral error is of the same order of magnitude as a large number of puffs -- so this is *at least as good* in a sense as the discrete sum of puffs method, given that we can vary the number of puffs to always make it a better/worse approximation



![svg](/images/intpuff2_successive_approximations_files/output_14_0.svg)


In the crosswind and vertical directions the sum of discrete puffs approximation works decidedly less well, at least at this slice in the cloud, while the approximate integral still works relatively well. I would say it is still *at least as good* as a sum of discrete puffs for a suitably large number of puffs.


![svg](/images/intpuff2_successive_approximations_files/output_16_0.svg)


![svg](/images/intpuff2_successive_approximations_files/output_17_0.svg)



This is, of course, very particular to that point downwind of the release. As we move closer to the origin the integral approximation gets worse, but then so does the sum of discrete puffs model. Especially for a low number of puffs: they become visibly discrete. I think this reinforces that, at least for class F stability, this approximation is in the same ball park as summing over a set discrete Gaussian puffs.


![svg](/images/intpuff2_successive_approximations_files/output_19_0.svg)



### Compute time

Model error is not the only factor in deciding upon an approximation. Since `QuadGK` exists we have to ask ourselves, why would we not always use it? We can answer that by benchmarking the three approaches at a particular point of interest (I don't think the choice of point impacts the calculations at all)


```julia
using BenchmarkTools: @benchmark

# point of interest
x₁ = 100
y₁ = σy(x₁)
z₁ = σz(x₁)
t₁ = x₁
```


Starting with the full numerical integration of the model, this is the *time to beat*. Any approximation that takes longer than ~30&mu;s is literally pointless: it generates worse results and takes longer.


```julia
@benchmark ∫pf(x₁,y₁,z₁,t₁; Δt=10)
```


![png](/images/intpuff2_successive_approximations_files/benchmark1.png)


As we expect, the sequence of discrete puffs is much faster for fewer puffs, and adding an order of magnitude more puffs increases the time by an order of magnitude. At around *n=100* we are no longer gaining anything over the full numerical integration. So, if the near-field matters a lot to you, then this is probably not a great approximation as the number of puffs required to approximate the full numerical integration well takes longer than just doing the integration.


```julia
@benchmark Σpf(x₁,y₁,z₁,t₁; Δt=10, n=10)
```


![png](/images/intpuff2_successive_approximations_files/benchmark2.png)


```julia
@benchmark Σpf(x₁,y₁,z₁,t₁; Δt=10, n=100)
```


![png](/images/intpuff2_successive_approximations_files/benchmark3.png)


Finally we have the integral approximation. This takes ~1/50th the time as the full numerical integration and, by the results above, it potentially performs just as well as the discrete puff approximation. In the examples above it was doing as well as discrete puff approximations that are too large to be worthwhile.


```julia
@benchmark ∫pf_approx(x₁,y₁,z₁,t₁; Δt=10)
```


![png](/images/intpuff2_successive_approximations_files/benchmark4.png)


I also have put no effort into optimizing any of this code, so take this with a grain of salt. Like the examination of the model error this is hardly rigorous, it is more suggestive than anything. It is possible that one could dramatically improve the discrete puff model, or re-write how the models are calculated to be more performant than I have. I prefer to write code that is easy for me to read, and re-uses things, but that does not necessarily translate into fast.

## Conclusions

I think it's worth noting that calculations that take on the order of tens of microseconds, on my crappy old laptop, are *fast*. To make the various plots required calculating the concentration at hundreds of points and my laptop did it all in the blink of an eye. I would say the first choice, all things being equal, would be simply to use the `QuadGK` model and call it a day. In terms of lines of code it is certainly short, all the heavy lifting is being done by the library. It also best captures *what you are trying to achieve*.

If you are doing a huge number of calculations, and can tolerate some model error, then the integral approximation is a good choice. It is the fastest and can perform as well as the discrete puff model. That said, there is an elephant in the room: The two integral approaches strictly require that all of the puffs are moving along the same line, at the same speed. For a great many chemical release scenarios that is entirely the set of assumptions being made, so it works great. However, for more complex atmospheric conditions -- with variable windspeed and direction -- then they don't work at all. Or, at least, it is not obvious to me how to adapt them to work. A slightly tweaked discrete puff model, tracking each puff's individual center location and windspeed, would be quite easy to implement, giving a more *flexible* model overall. This is in fact the how several more complicated atmospheric dispersion modeling tools work.


For a complete listing of code used to generate data and figures, please see the [corresponding julia notebook](https://nbviewer.org/github/aefarrell/aefarrell.github.io/blob/main/_notebooks/2023-01-15-intpuff2_successive_approximations.ipynb)
{: .notice--info}


---
