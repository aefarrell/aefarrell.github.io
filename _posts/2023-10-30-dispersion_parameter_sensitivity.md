---
title: "Messing around with model parameters"
last_modified_at: 2024-01-07
toc: true
toc_label: "Contents"
toc_sticky: true
comments: true
categories:
  - examples
tags:
  - dispersion modeling
tagline: "The importance of choosing the right references"
header:
  overlay_image: /images/gaussian_dispersion_example_files/veeterzy-unsplash-header.jpg
  caption: "Photo by [veeterzy](https://unsplash.com/@veeterzy?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText) on [Unsplash](https://unsplash.com/)"
---

# Messing around with model parameters

Recently I added some alternative correlations to [GasDispersion.jl](https://github.com/aefarrell/GasDispersion.jl), the julia package I put together for basic chemical release modeling, and I thought it would be worthwhile to circle back and look at some of those in more depth.

Typically, when evaluating various release scenarios, key pieces of the model are specified in advance and each scenario uses the same set of assumptions: comparing apples to apples. For a Gaussian plume dispersion model there are really three key correlations used for the model parameters: the wind-speed profile, crosswind dispersion, and vertical dispersion. Correlations for each of these are given in the standard references and there is not, to my mind, any deep reason to prefer one reference over the another. Besides maintaining consistency with other modeling or perhaps with industry practice in a particular area.

This raises the obvious question: how much does it matter which reference you use? Usually one takes the results of a Gaussian plume model with a fair grain of salt, these are "order of magnitude" estimates really. That's what I'm going to look at here.


## Windspeed

The windspeed correlations I am looking at here are the basic power law

$$ u = u_R \left( z \over z_R \right)^p $$

where *u<sub>R</sub>* is the known windspeed at a reference height *z<sub>R</sub>* and *p* is a parameter that depends upon the [Pasquill stability class](https://en.wikipedia.org/wiki/Outline_of_air_pollution_dispersion#The_Pasquill_atmospheric_stability_classes). There are more complex models that incorporate the surface roughness, Monin-Obukhov mixing length, and other measures of stability, they are beyond this analysis.

There are three different standard references used in `GasDispersion.jl` for windspeed: the default which comes from [Spicer and Havens](#spicer-1989), the correlations used by the [EPA Industrial Source Complex (ISC3)](#epa-1995) dispersion models, and the correlations given in the various [CCPS guidance documents](#ccps-1999)

The ISC3 and CCPS correlations are divided into urban and rural terrain and are exactly the same correlations for the unstable classes. They appear to be the correlations given in Hanna et al ([1982](#hanna-1982)). They also bracket the default correlation. Clearly whether or not the terrain is urban is significant, it can lead to a 20-30% difference in estimated windspeed (depending upon elevation).

![svg](/images/dispersion_parameter_sensitivity_files/output_6_0.svg)

For the stable atmospheres the ISC3 and CCPS *rural* correlations are the same. However they are very different for *urban* terrain and they no longer bracket the default correlation. The CCPS *urban* correlations are the same as Hanna et al ([1982](#hanna-1982)), the ISC3 correlations use the parameter *p* = 0.30 and no reference is given in the model specification so I don't know why.

For an urban release scenario, whether or not one choses the default, the ISC3 urban, or the CCPS urban correlation can lead to a 300% difference in windspeed (for class F stability, depending on elevation). Which is a pretty large difference.

![svg](/images/dispersion_parameter_sensitivity_files/output_8_0.svg)


## Plume Dispersion

The more diverse sets of correlations are for the plume dispersion parameters, the crosswind and vertical dispersion. To some extent this is because the early work[<sup id="fnref-1">1</sup>](#fn-1) presented the dispersion parameters graphically and many subsequent authors generated their own curves to fit these plots.

{% capture footnote-1 %}
<a name="fn-1"><strong>1</strong></a>: [Turner](#turner-1970), *Workbook of Atmospheric Dispersion Estimates*, 8-9. [↩](#fnref-1)
{% endcapture %}

<div class="notice">
  {{ footnote-1 | markdownify }}
</div>

### Crosswind Dispersion

Crosswind dispersion can be divided into the various attempts at fitting the curves presented graphically by [Turner](#turner-1970) and those based on Briggs' urban and rural correlations[<sup id="fnref-2">2</sup>](#fn-2)

{% capture footnote-2 %}
<a name="fn-2"><strong>2</strong></a>: [Briggs](#briggs-1973), *Diffusion Estimation for Small Emissions*, 38. Note that the correlations are given with respect to half-width/half-depth. [↩](#fnref-2)
{% endcapture %}

The default correlation is a simple set of correlations of the form

$$ \sigma_y = a x^b $$

which attempts to fit the [Turner](#turner-1970) curves.

The CCPS correlations are from [Briggs](#briggs-1973)) and the ISC3 *urban* correlations are from Briggs as well, the ISC3 *rural* correlations are something else entirely but I suspect are intended to fit the [Turner](#turner-1970) curves. The correlations from the [TNO yellow book](#bakkum-2005) are also a different attempt at fitting the Turner curves. What `GasDispersion,jl` gives as "Turner" is the fit to the Turner curves given in [Lees](#lees-1996).

![svg](/images/dispersion_parameter_sensitivity_files/output_11_0.svg)

Zooming in on the class F curves is illustrative of the lot: most of the lines overlap and hew pretty close to the curve-fit for [Turner](#turner-1970) with the exception of the Briggs' urban/rural correlations. The biggest impact on these model parameters is whether or not a rural/urban terrain is used or not. *Note* these are log-log plots.

![svg](/images/dispersion_parameter_sensitivity_files/output_13_0.svg)

<div class="notice">
  {{ footnote-2 | markdownify }}
</div>

### Vertical Dispersion

The vertical dispersion correlations are decidedly more varied. Varied enough that I'm just going to show them all at full scale[<sup id="fnref-3">3</sup>](#fn-3)

{% capture footnote-3 %}
<a name="fn-3"><strong>3</strong></a>: The correlations given in the [AIChE/CCPS](#ccps-1999) *Guidelines for Consequence Analysis* for urban conditions has typos in the class A, B and D correlations, I have corrected them here to match the Briggs correlations on which they are supposed to be based. [↩](#fnref-3)
{% endcapture %}


![svg](/images/dispersion_parameter_sensitivity_files/output_17_0.svg)

![svg](/images/dispersion_parameter_sensitivity_files/output_18_0.svg)

![svg](/images/dispersion_parameter_sensitivity_files/output_19_0.svg)

![svg](/images/dispersion_parameter_sensitivity_files/output_20_0.svg)

![svg](/images/dispersion_parameter_sensitivity_files/output_21_0.svg)

![svg](/images/dispersion_parameter_sensitivity_files/output_22_0.svg)


For some of these there is *an order of magnitude* spread in vertical dispersion, depending on which model happens to be used. Even when looking only at the correlations that are "universal", i.e. are not for either urban or rural terrains. From this alone one would expect that the concentration profiles would vary by a large amount, depending on which set of correlations one used to model a given scenario.

<div class="notice">
  {{ footnote-3 | markdownify }}
</div>

## An Example

Just to give an example of how this works out, lets look at the emissions from a large stack. I happened to have picked the stack for a large power plant in the Edmonton area: TransAlta's Sundance station. This power plant is on the shores of Lake Wabamun and is pretty rural, it has several stacks but let's consider only Stack 2 and examine the dispersion of SO<sub>2</sub> emissions.

From Alberta's [AEIR Air Emission Rates dataset](https://open.alberta.ca/opendata/aeirairemissionrates) we can pull the mass emission rates for SO<sub>2</sub> as well as the relevant stack dimensions. *Note* this dataset is from 2018 and thus may not represent the current operations at Sundance.


```julia
# TransAlta Sundance - Stack 2
m = 3200/3600 # mass emission rate: 3200kg/h in kg/s
h = 155.5 # stack height, m
d = 7.3   # stack diameter, m
v = 35.6  # stack exit velocity, m/s
T = 439.7 # stack exit temperature, K
```

For the sake of modeling let's assume a class D atmospheric stability with a windspeed at 10m of 2m/s. The atmosphere is otherwise at standard state.


```julia
# assumed weather conditions
uᵣ  = 2  # windspeed, m/s
zᵣ = 10 # windspeed elevation, m
stability = ClassD

# standard state
Pₛ = 101325 # Pa
Tₛ = 273.15 # K
```

We can construct the relevant scenario for `GasDispersion.jl` directly.


```julia
r = VerticalJet(m, Inf, d, v, h, Pₛ, T, 0.0)

a = SimpleAtmosphere(pressure=Pₛ, temperature=Tₛ, windspeed=uᵣ, windspeed_height=zᵣ, stability=stability)

# a dummy substance, since I know a gaussian plume doesn't require any material
# properties I have just left them as NaNs
SO2 = Substance(name=:SulfurDioxide,molar_weight=0.064066,liquid_density=1,boiling_temp=1,
                latent_heat=1,gas_heat_capacity=1,liquid_heat_capacity=1)

scn = Scenario(SO2,r,a)
```




    Substance: SulfurDioxide 
        MW: 0.064066 kg/mol 
        P_v: GasDispersion.Antoine{Float64}(0.007705368698167287, 0.007705368698167287, 0.0) Pa 
        ρ_g: 2.7095140841291006 kg/m^3 
        ρ_l: 1 kg/m^3 
        T_ref: 288.15 K 
        P_ref: 101325.0 Pa 
        k: 1.4  
        T_b: 1.0 K 
        Δh_v: 1 J/kg 
        Cp_g: 1 J/kg/K 
        Cp_l: 1 J/kg/K 
    VerticalJet release:
        ṁ: 0.8888888888888888 kg/s 
        Δt: Inf s 
        d: 7.3 m 
        u: 35.6 m/s 
        h: 155.5 m 
        P: 101325.0 Pa 
        T: 439.7 K 
        f_l: 0.0  
    SimpleAtmosphere atmosphere:
        P: 101325.0 Pa 
        T: 273.15 K 
        u: 2.0 m/s 
        h: 10.0 m 
        rh: 0.0 % 
        stability: ClassD 




The Gaussian plume model is then given by the following, neglecting the effect of plume rise.


```julia
conc = plume(scn, GaussianPlume; plumerise=false);
```

Plotted below are the results for every equation set, at near ground level (at basically "my head" level). Clearly the urban/rural choice is quite important, leading to a ~4&times; greater maximum concentration. The TNO correlations, which uses the default correlation for windspeed and the TNO correlations for the crosswind and vertical dispersion, leads to less dispersion and thus a greater maximum concentration relative to the rest.

![svg](/images/dispersion_parameter_sensitivity_files/output_33_0.svg)

Plotted below is the 172ppbv isopleth, the 1-hr [Ambient Air Quality Objective (AAQO)](https://open.alberta.ca/publications/9781460134856) for SO<sub>2</sub> in Alberta. As we would expect, the correlations that lead to a higher maximum concentration correspond to less overall dispersion and the isopleth is quite a bit smaller for the urban versus rural case and the TNO versus the remaining cases. The scale is in kilometers so this is quite a large difference in area.

![svg](/images/dispersion_parameter_sensitivity_files/output_35_0.svg)

The above was assuming no plume rise, however the relative differences are much more pronounced when plume rise is included.


```julia
conc = plume(scn, GaussianPlume; plumerise=true);
```

Plotted below is the same downwind concentration plot as above, but incorporating the Briggs' plume rise model. Since this leads to a greater overall dispersion, the concentration is much smaller (everything is well below the AAQO at ground level, which is good news). However this adds another dimension along which the models can vary: plume rise is a function of windspeed, and overall dispersion is a function of plume rise. These different sets of correlations lead to the plume rising to a different elevation, and also dispersing to a differing degree, magnifying the differences between them. In this case there is up to a ~30&times; difference between the max concentrations predicted between the urban and rural case.

![svg](/images/dispersion_parameter_sensitivity_files/output_39_0.svg)


## Final Thoughts

I think the above illustrates the necessity of picking a standard set of correlations for use when screening scenarios at a particular plant (e.g. using either the CCPS urban or rural correlations as appropriate for the area around the plant) and being careful to keep these consistent. It also shows how seriously one should take the exact values generated by the models: not very. The dispersion model results are highly sensitive to the choice of correlations, and they are also quite sensitive to the other assumptions that go into a release scenario (e.g. atmospheric stability, wind-speed, mass emission rate). The results are really order of magnitude at best.

It is often the case that chemical plants are situated at the periphery of cities, in areas that blur the line between "urban" and "rural". Also, cities grow and industrial areas fill in. A plant that was essentially rural may, overtime, fill in such that the urban correlations better represent the area. I think it is worth comparing the urban/rural models for a range of plausible results and considering whether assumptions made in the past about the area around the plant are still valid given changes in the area.

There are other correlations, for wind-speed and for dispersion, that take into account the local surface roughness which could be used instead and the sensitivity to the models to assumptions about surface roughness could be evaluated. This would likely lead to a smaller range of values, and give a path for updating the screening model as the area around the plant changes (update the assumed surface roughness and re-run).

For a complete listing of code used to generate data and figures, please see the [corresponding julia notebook](https://github.com/aefarrell/aefarrell.github.io/blob/main/_notebooks/2023-10-30-dispersion_parameter_sensitivity.ipynb)
{: .notice--info}

## References
    
+ <a name="ccps-1999">AIChE/CCPS</a>. 1999. *Guidelines for Consequence Analysis of Chemical Releases.* New York: American Institute of Chemical Engineers
+ <a name="bakkum-2005">Bakkum</a>, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands.
+ <a name="briggs-1973">Briggs</a>, Gary A. 1973. *Diffusion Estimation for Small Emissions. Preliminary Report.* United States. [doi:10.2172/5118833](https://doi.org/10.2172/5118833)
+ <a name="epa-1995">EPA</a>. 1995. *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models, vol 2.* United States Environmental Protection Agency [EPA-454/B-95-003b](https://gaftp.epa.gov/Air/aqmg/SCRAM/models/other/isc3/isc3v2.pdf)
+ <a name="griffiths-1994">Griffiths</a>, R. F. 1994. "Errors in the use of the Briggs parameterization for atmospheric dispersion coefficients." *Atmospheric Environment* 28(17):2861-2865 [doi:10.1016/1352-2310(94)90086-8](https://doi.org/10.1016/1352-2310(94)90086-8)
+ <a name="hanna-1982">Hanna</a>, Steven R., Gary A. Briggs, and Rayford P. Hosker Jr. 1982. *Handbook on atmospheric diffusion*. United States. [doi:10.2172/5591108](https://doi.org/10.2172/5591108).
+ <a name="lees-1996">Lees</a>, Frank P. 1996. *Loss Prevention in the Process Industries, 2nd ed.* Oxford: Butterworth-Heinemann
+ <a name="spicer-1989">Spicer</a>, Thomas O. and Jerry A. Havens. 1989. *User's Guide for the DEGADIS 2.1 Dense Gas dispersion Model* United States Environmental Protection Agency EPA-450/4-89-019
+ <a name="seinfeld-1986">Seinfeld</a>, John H. 1986. *Atmospheric Chemistry and Physics of Air Pollution.* New York: John Wiley and Sons
+ <a name="turner-1970">Turner</a>, D. Bruce. 1970. *Workbook of Atmospheric Dispersion Estimates.* United States Environmental Protection Agency.

