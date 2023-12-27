---
title: "Taking a second look at the Britter-McQuaid model"
last_modified_at: 2023-12-26
toc: true
toc_label: "Contents"
toc_sticky: true
comments: true
categories:
  - notes
tags:
  - dispersion modeling
  - explosion
  - hazard screening
tagline: "Re-evaluating plume extents and determing the explosive mass"
header:
  overlay_image: /images/britter_mcquaid_files/christianbuehner-unsplash-header.jpg
  overlay_filter: 0.25
  caption: "Photo by [christian buehner](https://unsplash.com/es/@christianbuehner?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText) on [Unsplash](https://unsplash.com/)"
---


# Taking a second look at the Britter-McQuaid model

I recently spent some time looking in detail at the Britter-McQuaid workbook model for dense gas dispersion and I thought the plume model deserved some extra attention. Firstly because I believe there is an error in the plume dimensions, and secondly because I think an important feature of top-hat models is often neglected and the Britter-McQuaid workbook model should be used more.


As a re-cap the Britter-McQuaid model[<sup id="fnref-1">1</sup>](#fn-1) is a series of correlations for the dispersion of denser than air gases. These are given as a series of correlation curves and the typical procedure is to interpolate the downwind distance to the concentration of interest, for example to the Lower Flammability Limit (LFL). The model also gives some equations for estimating the plume horizontal and vertical dimensions, where conventional practice is to assume the plume has a rectangular cross-section and a uniform concentration.

{% capture footnote-1 %}
<a name="fn-1"><strong>1</strong></a>: Britter and McQuaid, *[Workbook on the Dispersion of Dense Gases](#britter-1988)*. [↩](#fnref-1)
{% endcapture %}

<div class="notice">
  {{ footnote-1 | markdownify }}
</div>

## A motivating example

Just to have some numbers to look at, I am going to use a scenario adapted from the Burro series of trials of LNG dispersion[<sup id="fnref-2">2</sup>](#fn-2). The release conditions are:
+ release temperature: -162°C
+ release rate: 0.23 m³/s (liquid)
+ release duration: 174 s
+ windspeed at 10m: 10.9 m/s
+ LNG liquid density (at release conditions): 425.6 kg/m³
+ LNG gas density (at release conditions): 1.76 kg/m³

The goal is to find the distance to the Lower Flammability Limit (LFL) which is 5%(v/v) and ultimately work out the extent of the plume and total explosive mass.

{% capture footnote-2 %}
<a name="fn-2"><strong>2</strong></a>: Adapted from AIChE/CCPS, *[Guidelines for Consequence Analysis](#ccps-1999)*, 122. [↩](#fnref-2)
{% endcapture %}


```julia
using Unitful

Tₐ  = 288.15u"K"      # ambient air temperature 
ρₐ  = 1.225u"kg/m^3"  # density of air at 15°C and 1atm
u₁₀ = 10.9u"m/s"      # windspeed at 10m 

ρₗ  = 425.6u"kg/m^3"  # liquid density of LNG, given
ρᵥ  = 1.76u"kg/m^3"   # vapour density of LNG, given
ṁ   = ρₗ*0.23u"m^3/s" # mass release rate
Tᵣ  = (273.15 - 162)u"K" # boiling point of LNG, given
LFL = 0.05 # lower flammability limit, volume fraction

Qₒ = ṁ/ρᵥ # gas volumetric flowrate: mass flowrate divided by gas density
```

    55.618181818181824 m^3 s^-1



First calculate the critical length, *D*, and the dimensionless parameter *&alpha;* for the model


```julia
D = √(Qₒ/u₁₀)

g = 9.806u"m/s^2"
gₒ = g * (ρᵥ - ρₐ )/ ρₐ
α = 0.2*log10(gₒ^2 * Qₒ / u₁₀^5)
```


    -0.4356993075502021



Then, using digitized curves[<sup id="fnref-3">3</sup>](#fn-3), work out the points for the linear interpolation in terms of $\beta = \log_{10}(x/D)$

{% capture footnote-3 %}
<a name="fn-3"><strong>3</strong></a>: Adapted from AIChE/CCPS, *[Guidelines for Consequence Analysis](#ccps-1999)*, 118. [↩](#fnref-3)
{% endcapture %}


```julia
Cs = [ 0.1,         0.05,        0.02,        0.01,        0.005,       0.002]
βs = [ 0.24*α+1.88, 0.36*α+2.16, 0.45*α+2.39, 0.49*α+2.59, 0.59*α+2.80, 0.39*α+2.87]
```


    6-element Vector{Float64}:
     1.7754321661879513
     2.0031482492819275
     2.1939353116024094
     2.376507339300401
     2.5429374085453804
     2.7000772700554214



These points only cover the middle region of the concentration curve, where the concentration ratio, $ { c_m \over c_0 } $, is between 0.1 and 0.002, there is a near-field correlation that needs to be connected for concentration ratios &gt;0.1


```julia
function Cm_nf(x′)
    if x′ > 0
        return 306/(306 + x′^2)
    else
        return 1.0
    end
end

xnf = 30
βnf = log10(xnf)
Cnf = Cm_nf(xnf)
```


    0.2537313432835821



And a far field correlation for when the concentration ratio is &lt;0.002 which is basically just continuing the curve from the last point but such that the concentration decays with *1/x<sup>2</sup>*


```julia
xff = 10^(maximum(βs))
A = minimum(Cs)*xff^2

function Cm_ff(x′; A=A)
    return A/x′^2
end
```


Finally, putting together the pieces: near field correlation, a linear interpolation for the middle of the concentration curve, and a far field correlation, to form the complete concentration function, along with a correction for non-isothermal releases (of which this is an example)


```julia
using Interpolations

itp = interpolate( ([βnf; βs],), [Cnf; Cs], Gridded(Linear()) )
```


    7-element interpolate((::Vector{Float64},), ::Vector{Float64}, Gridded(Linear())) with element type Float64:
     0.2537313432835821
     0.1
     0.05
     0.02
     0.01
     0.005
     0.002


```julia
function Cm(x::Quantity; xnf=xnf, xff=xff, D=D, T′=Tᵣ/Tₐ)
    x′ = x/D
    c′ = if x′ < xnf
        Cm_nf(x′)
    elseif xnf ≤ x′ < xff
        itp(log10(x′))
    else
        Cm_ff(x′)
    end
    
    c = c′ / (c′ + (1 - c′)*T′)
    
    return c
end
```


    
![svg](/images/britter_mcquaid_files/output_14_0.svg)
    



If all one needs is the distance to the LFL there is an easier way of doing this: interpolate the concentrations to find the *&beta;* corresponding to the LFL (after applying the non-isothermal correction). However, if one also requires the plume dimensions the concentration profile is required.

From the concentration profile calculating the downwind distance to the LFL is very straight-forward.


```julia
using Roots

xn = find_zero((x) -> Cm(x) - LFL, (300,400).*1u"m", Roots.Brent())
```


    354.5630187009715 m

<div class="notice">
  {{ footnote-2 | markdownify }}

  {{ footnote-3 | markdownify }}
</div>

## Looking again at plume dimensions

At first glance the workbook seems to be giving the user everything they need to workout the size of the plume, giving the following diagram

![image.png](/images/britter_mcquaid_files/att1.png)

and the following relations for the labeled distances

$$ L_U = {D \over 2} + 2 l_b $$

$$ L_{Ho} = D + 8 l_b $$

$$ L_H = L_{Ho} + 2.5 \sqrt[3]{ l_b x^2 } $$

with the buoyancy scale *l<sub>b</sub>* defined as
$$ l_b = { { g_o Q_o } \over u_{ref}^{3} } $$


```julia
lb = (gₒ*Qₒ)/u₁₀^3
```


    0.18392758812310803 m




```julia
Lᵤ  = D/2 + 2lb
```


    1.4973003373658906 m


```julia
Lₕₒ = D + 8lb
```


    3.7303110272242135 m


```julia
Lₕ(x) = Lₕₒ + 2.5∛(lb*x^2)
```


### Upwind region

The curve given for *L<sub>H</sub>* for *x* &gt; 0 is **not** the curve for *x* &lt; 0, the upwind extent of the plume. This is the blue curve in the figure below. The orange curve is slightly adjusting *L<sub>H</sub>* such that for *x* &lt; 0 the second term is subtracted (so the curve actually converges to zero instead of blowing up to +&infin; as *x* &rarr; -&infin;). The black dots are points taken from the diagram given by Britter and McQuaid, using a graph digitizer and scaling to the actual *L<sub>Ho</sub>* and *L<sub>U</sub>*. Clearly the given curve for *L<sub>H</sub>* is not at all what is shown in the diagram for the upwind region.

A conservative approach to estimating the size of the upwind extent is to assume *L<sub>H</sub>* = *L<sub>Ho</sub>* for *L<sub>U</sub>* &lt; *x* &lt; 0, i.e. making the upwind region a rectangle of width *L<sub>Ho</sub>* and length *L<sub>U</sub>* [<sup id="fnref-4">4</sup>](#fn-4). This is the green curve in the figure below.

{% capture footnote-4 %}
<a name="fn-4"><strong>4</strong></a>: This is the approach taken in the TNO *Yellow Book* (Bakkum and Duijm "[Vapour Cloud Dispersion](#bakkum-2005)"). [↩](#fnref-4)
{% endcapture %}

Alternatively one could "fit" a curve to hit the end points while also having the same power of *x*: $ L_H = L_{Ho} \left( {x + L_U} \over L_U \right)^{2/3} $ where *L<sub>U</sub>* &lt; *x* &lt; 0, this at least retains the same general shape and is the red curve in the figure below. I think this should be taken with the giant caveat that I don't know if insisting on the same power law is truly justified.

    
![svg](/images/britter_mcquaid_files/output_23_0.svg)
    



For most typical cases I would think the upwind region would be a small component of the overall plume and taking the conservative, rectangle, approach would be a small error.

<div class="notice">
  {{ footnote-4 | markdownify }}
</div>

### Vertical extent

The vertical extent is not given on the diagram, but an equation is given in the text, with the note that this comes from continuity, however I think this is incorrect.

$$ L_V = {Q_o \over {u_{ref} L_H} } = { D^2 \over L_H } $$

Suppose a steady state plume with a system boundary such that the plume is sliced along the *y-z* plane at some downwind distance *x*. All of the mass entering the plume, from the source, exits the plume through this plane

![image.png](/images/britter_mcquaid_files/att2.png)

Consider the steady state mass balance

$$ \textrm{mass in} = \textrm{mass out} \\
c_o Q_o = \iint_A c u \,dA = \int_{0}^{\infty} \int_{-\infty}^{\infty} c(x,y,z) u(x,y,z) \,dy \,dz $$

By the nature of a top-hat model the plume cross section is a rectangle with half-width *L<sub>H</sub>* and height *L<sub>V</sub>* and the concentration everywhere *inside* the rectangle is *c<sub>m</sub>*. Assuming a constant advection velocity, *u*, the integral can be simplified to

$$ \iint_A c u \,dA = c_m u \int_{0}^{L_V} \int_{-L_H}^{L_H} \,dy \,dz = 2 c_m u L_H L_V $$

The steady state mass balance is then

$$ c_o Q_o = 2 c_m u L_H L_V $$

and the vertical extent can be solved for with some simple re-arrangement

$$ L_V = { { c_o Q_o } \over { 2 c_m u L_H } } = {1 \over 2}{ c_o \over c_m } { Q_o \over {u L_H} } $$

Setting the advection velocity of the plume to the reference windspeed gives

$$ L_V = {1 \over 2}{ c_o \over c_m } { Q_o \over {u_{ref} L_H} } = {1 \over 2} { c_o \over c_m } { D^2 \over L_H }$$


```julia
Lᵥ(x) = D^2/(2*Cm(x)*Lₕ(x))
```


This is definitely similar to what is given by Britter and McQuaid but with two big differences:
+ it depends upon the concentration
+ it is divided by two

The last point could equally be a mistake in the diagram (I have no real way of checking) as while the diagram shows *L<sub>H</sub>* as the plume *half-width*, the text simply refers to it as the "lateral plume extent", which is ambiguous -- do they mean the entire lateral extent or from the center-line of the plume? 

The TNO Yellow Book gives a different equation[<sup id="fnref-5">5</sup>](#fn-5) for the vertical extent:

{% capture footnote-5 %}
<a name="fn-5"><strong>5</strong></a>: Bakkum and Duijm "[Vapour Cloud Dispersion](#bakkum-2005)," equation 4.104. [↩](#fnref-5)
{% endcapture %}

$$ L_V = {1 \over 2} { Q_o \over {u_{ref} L_H} } = {1 \over 2} { D^2 \over L_H }$$

Which clearly follows from assuming *L<sub>H</sub>* is the *half-width*, and the corresponding figure is labeled as such (using the same equation for *L<sub>H</sub>* as Britter and McQuaid). But it doesn't depend upon concentration.

I think the vertical extent has to depend upon the concentration as otherwise mass will simply disappear from the plume as it extends downwind. There is also the obvious problem that since the plume lateral extent monotonically increases, and the vertical extent is inversely related to it, the vertical extent is monotonically decreasing. In fact it becomes vanishingly small quite quickly. This entirely the opposite of what is observed with actual dense plume dispersion.

This can be seen most clearly in the following figure in which the vertical extent is shown as a function of downwind distance along with the mass flowrate in the plume (i.e. $ c_m u A $ )


    
![svg](/images/britter_mcquaid_files/output_28_0.svg)
    



I think it is fairly obvious that both the Britter-McQuaid and TNO models give silly answers for the vertical extent. Though the corrected curve, the green curve, clearly has problems too: it has an odd bumpiness, as a result of the linear interpolation, and it is also too small due to both assuming the concentration everywhere is equal to the ground level concentration and due to an overly large advection velocity (the windspeed at 10m is quite a bit higher than the windspeed at ~1m).

An alternative approach to using the reference windspeed as the advection velocity is to assume the advection velocity is some constant fraction of the reference velocity, e.g. $u = 0.4 u_{ref} $, which is what Britter and McQuaid use for the instantaneous model.

Another alternative might be to use an *average* windspeed, *&umacr;* over cross-section of the plume as the advection velocity, assuming windspeed is only a function of height.

$$ \bar{u} = { { \iint_A u \,dA } \over A } =  { { \int_{0}^{L_V} u(z) \,dz } \over L_V } $$

Assuming the windspeed follows a powerlaw distribution $ u = u_{ref} \left( z \over z_{ref} \right)^p $ gives

$$ \bar{u} = { { \int_{0}^{L_V} u(z) \,dz } \over L_V } \\
= {1 \over L_V} \int_{0}^{L_V} u_{ref} \left( z \over z_{ref} \right)^p \,dz \\
= { u_{ref} \over {p+1} } \left( L_V \over z_{ref} \right)^p $$

plugging it into the simple mass balance

$$ c_o Q_o = c_m \bar{u} A \\
= c_m {u_{ref} \over {p+1} } \left( L_V \over z_{ref} \right)^p { 2 L_H L_V } $$

re-arranging to solve for *L<sub>V</sub>*

$$ L_V = \left( { {p+1} \over 2 } { c_o \over c_m } z_{ref}^p {Q_o \over {u_{ref} L_H} } \right)^{1 \over {p+1} } \\
= \left( { {p+1} \over 2 } { c_o \over c_m } z_{ref}^p {D^2 \over L_H } \right)^{1 \over {p+1} }$$

The red curve in the figure above is this model, using *p* = 0.15[<sup id="fnref-6">6</sup>](#fn-6)

{% capture footnote-6 %}
<a name="fn-6"><strong>6</strong></a>: AIChE/CCPS, *[Guidelines for Consequence Analysis](#ccps-1999)*, 83. [↩](#fnref-6)
{% endcapture %}

This could also be done using the logarithmic windspeed curve $ u = {u_{\star} \over \kappa} \log \left( z \over z_) \right) $ where $u_{\star}$ is the friction velocity and *z<sub>0</sub>* is the roughness length. Though I don't imagine the expression would work out as nicely.

<div class="notice">
  {{ footnote-5 | markdownify }}

  {{ footnote-6 | markdownify }}
</div>

### Recommendations

For the upwind region, assuming a simple rectangular prism with length *L<sub>U</sub>*, width *2L<sub>Ho</sub>*, height *L<sub>Vo</sub>* and uniform concentration *c<sub>o</sub>* is a conservative approach. Likely the plume downwind of the source will be much larger than the upwind area and so this will be a small overestimate.

The simple mass balance approach to calculating the plume height is a reasonable approach if one simply wants to reference Britter and McQuaid and not have to justify additional assumptions. It is not what is given in the text, but it is what is *described* in the text. The other models for plume height may be more realistic, in the sense that they represent more realistic advection velocities, and will give larger explosive masses for the plume, however they have not been validated against any actual data. That validation may be a worthwhile exercise but is well beyond the scope of this blog post.

## Calculating the explosive mass

The explosive mass in the cloud is the given by the volume integral

$$ m_e = \iiint_V c dV $$

where *V* is defined as the region where *c* &ge; *LFL*[<sup id="fnref-7">7</sup>](#fn-7).

{% capture footnote-7 %}
<a name="fn-7"><strong>7</strong></a>: Some sources recommend *1/2 LFL*. [↩](#fnref-7)
{% endcapture %}


Using the concentration profile and the plume extents, we could work out the function `c(x,y,z)` such that the concentration is returned if we are:
+ within the plume, and
+ the concentration is &ge; *LFL*

To determine the explosive mass in the downwind region this might be done by the following

```julia

cₒ = ustrip(u"kg/m^3", ṁ/Qₒ)

Lₕ(x::Number) = ustrip(u"m", Lₕ(x*1u"m"))
Lᵥ(x::Number) = ustrip(u"m", D)^2/(2*Cm(x)*Lₕ(x))

function c(x,y,z; lim=LFL)
    c_ = Cm(x)
    
    if c_ ≥ lim
        if (abs(y) ≤ Lₕ(x)) && (z ≤ Lᵥ(x))
            return cₒ*c_
        else
            return 0.0
        end
    else
        return 0.0
    end
end

```

```julia

using HCubature: hcubature

x_min, x_max = 0, xn
y_min, y_max = -Lₕ(xn), Lₕ(xn)
z_min, z_max = 0, Lᵥ(xn)

m_e, err = hcubature( c, [x_min, y_min, z_min], [x_max, y_max, z_max])

```

This is a pretty tedious integration, is very inefficient, and doesn't take into account any of the *structure* of the model and it turns out that a top-hat model has some pretty convenient structure.

<div class="notice">
  {{ footnote-7 | markdownify }}
</div>

### A nice property of top hat models

Returning to the integral for the explosive mass, the plume can be divided into an upwind region (*x* &lt; 0) and a downwind region (*x* &ge; 0)

$$ m_e = \iiint_V c \,dV = m_{e,u} + m_{e,d} $$

with the explosive mass of the downwind region being

$$ m_{e,d} = \int_0^{x_n} \iint_A c \,dA \,dx $$

For a top-hat model, since the concentration at a given downwind distance is constant everywhere within the plume cross-section $\iint_A c dA = c_m A$, and, from a mass balance on the plume

$$ c_m A u = c_o Q_o \\
c_m A = { {c_o Q_o} \over u} $$

which is a constant, thus

$$ m_{e,d} = \int_0^{x_n} c_m A \,dx \\
= \int_0^{x_n} { {c_o Q_o} \over u} \,dx \\
= { {c_o Q_o} \over u} x_n $$

For the explosive mass of the upwind region a simple box model gives $ m_{e,u} = 2 c_o L_U L_{Ho} L_{Vo} $. Putting everything together[<sup id="fnref-8">8</sup>](#fn-8)

{% capture footnote-8 %}
<a name="fn-8"><strong>8</strong></a>: this is not specific to the Britter-McQuaid model, it works for *any* top hat model. [↩](#fnref-8)
{% endcapture %}

$$ m_e = 2 c_o L_U L_{Ho} L_{Vo} + { {c_o Q_o} \over u} x_n $$


This can be simplified greatly by setting the advection velocity to *u<sub>ref</sub>*

$$ m_e = 2 c_o L_U L_{Ho} L_{Vo} + { {c_o Q_o} \over u_{ref} } x_n \\
= 2 c_o L_U L_{Ho} {1 \over 2}{D^2 \over L_{Ho} } + c_o D^2 x_n \\
= c_o D^2 \left( L_U + x_n \right)$$


```julia
cₒ = ṁ/Qₒ

mₑ = cₒ*D^2*(Lᵤ+xn)
```


    3197.617661470163 kg


This very simple expression is the obvious strength of a top-hat model: it makes calculating the explosive mass incredibly easy[<sup id="fnref-9">9</sup>](#fn-9). It also retroactively justifies why the Britter McQuaid model is oriented around calculating *x<sub>n</sub>*: that's all you actually need[<sup id="fnref-10">10</sup>](#fn-10).

{% capture footnote-9 %}
<a name="fn-9"><strong>9</strong></a>: Woodward gives a variation on this for top-hat models more generally (Woodward, *[Estimating the Flammable Mass](#woodward-1998)*). [↩](#fnref-9)
{% endcapture %}

{% capture footnote-10 %}
<a name="fn-10"><strong>10</strong></a>: some sources recommend calculating the explosive mass as the region of the plume with the concentration *LFL* &le; *c* &le; *UFL*, in which case $m_e = c_o D^2 \left( x_{n,LFL} - x_{n,UFL} \right)$ [↩](#fnref-10)
{% endcapture %}


If this seems too good to be true, the integration can be performed numerically by taking

$$ \iint_A c \,dA = c_m \cdot 2L_H \cdot L_V $$


```julia
using QuadGK: quadgk

function ∫∫cdA(x)
    if Cm(x) ≥ LFL
        return cₒ*Cm(x)*(2Lₕ(x))*Lᵥ(x)
    else
        return 0.0u"kg/m"
    end
end

m_ed, err = quadgk(∫∫cdA, 0u"m", xn)

m_eu = 2*cₒ*Lᵤ*Lₕₒ*Lᵥ(0u"m")

m_eu + m_ed
```


    3197.617661470163 kg


Which is exactly the same.

Above I claimed the upwind region was "small" relative to the downwind region, this can be shown easily as the mass in each region is directly proportional to the length.


```julia
Lᵤ/(Lᵤ+xn)
```


    0.004205187316042018


Since the mass in the upwind region is &lt;0.5% of the total mass in the cloud, I think the simple box model is justified.

<div class="notice">
  {{ footnote-8 | markdownify }}

  {{ footnote-9 | markdownify }}

  {{ footnote-10 | markdownify }}
</div>

### Added complications

According to Britter and McQuaid the top-hat model generates an overly conservative plume extent and they recommend using given the lateral extent curve up to *2/3 x<sub>n</sub>* and after which connecting to *x<sub>n</sub>* using straight lines, as shown in the plume diagram. This makes the integration for explosive mass a little more complicated.

![image.png](/images/britter_mcquaid_files/att1.png)

For simplicity the plume can be divided into three regions, the upwind region (*x* &lt; 0), the downwind region up to the cutoff (0 &le; *x* &lt; 2/3 *x<sub>n</sub>*), and the downwind cutoff region (2/3 *x<sub>n</sub>* &le; *x* &lt; *x<sub>n</sub>* )

$$ m_{e, \textrm{cut off} } = m_{e,u} + m_{e,d1} + m_{e,d2} $$

The upwind region, *m<sub>e,u</sub>*, and the first downwind region *m<sub>e,d1</sub>* are already known, they are the same as above up to 2/3 *x<sub>n</sub>*. What is left to determine is the explosive mass in the cutoff region.

$$ m_{e,d2} = \int_{2/3 x_n}^{x_n} \iint_A c \,dA \,dx $$

The integral can be re-written to take advantage of *c<sub>m</sub>A* being an invariant for a top-hat model,

$$ m_{e,d2} = \int_{2/3 x_n}^{x_n} \iint_A c \,dA \,dx \\
= c_m A_{\textrm{original} } \int_{2/3 x_n}^{x_n} { A_{\textrm{cut off} } \over A_{\textrm{original} } } \,dx $$

Assuming the vertical extent remains unchanged in this operation, the ratio of areas is the same as the ratio of horizontal extents

$$ { A_{\textrm{cut off} } \over A_{\textrm{original} } } = { L_{H, \textrm{cut off} } \over L_{H, \textrm{original} } } $$

From some simple geometry, the horizontal extent is

$$ L_{H, \textrm{cut off} } = 3 L_{H, 2/3 x_n}  { {x_n - x} \over x_n }$$

Which then leads to

$$ m_{e,d2} = 3 c_o D^2 \int_{2/3 x_n}^{x_n} { L_{H, 2/3 x_n} \over L_H } { {x_n - x} \over x_n } \,dx $$

There is probably a closed form for this integral but it is just as easy to integrate that numerically.

$$ m_{e, \textrm{cut off} } = c_o D^2 \left( L_U + \frac{2}{3}x_n + 3 L_{H, 2/3 x_n} \int_{2/3 x_n}^{x_n} { 1 \over L_H(x) } { {x_n - x} \over x_n } \,dx   \right)$$


```julia
mₑ_cutoff = cₒ*D^2*(Lᵤ + (2/3)*xn 
    + 3*Lₕ((2/3)*xn)*quadgk( (x) -> (xn - x)/(xn*Lₕ(x)), (2/3)*xn, xn)[1] )
```


    2620.489605856347 kg


This works out to be about 20% less than the original explosive mass.


```julia
mₑ_cutoff/mₑ
```


    0.8195131136008078


## Final thoughts

I think the error in the vertical extent may have limited the apparent utility of the Britter-McQuaid model. Most references I have do use the Britter-McQuaid model, noting that it is "reasonably simple to apply, and produces results which appear to be as good as more sophisticated models"[<sup id="fnref-11">11</sup>](#fn-11), however they either claim that it is *only* good for calculating *x<sub>n</sub>* or gloss over how it could be used for anything else. The CCPS references seem consistent in neglecting to mention at all that the model can also estimate the plume extent. So, while I can't imagine I'm the first person to have noticed that the given equation for *L<sub>V</sub>* doesn't work, I have yet to encounter anyone *actually admitting it*.

{% capture footnote-11 %}
<a name="fn-11"><strong>11</strong></a>: AIChE/CCPS, *[Guidelines for Consequence Analysis](#ccps-1999)*, 122. [↩](#fnref-11)
{% endcapture %}

That said, the correction also seems obvious to me: one simply follows what is described in the text which is exactly how Britter and McQuaid calculated the cloud height for the instantaneous model (which is correct) in the same workbook. That the incorrect equation for *L<sub>V</sub>* is repeated in other references[<sup id="fnref-12">12</sup>](#fn-12), with only the *[TNO Yellow Book](#bakkum-2005)* making a correction, while still repeating a critical mistake, strikes me as very odd.

{% capture footnote-12 %}
<a name="fn-12"><strong>12</strong></a>: For example: Lees, *[Loss Prevention](#lees-1996)*, and Casal, *[Evaluation of the Effects of Consequences of Major Accidents in Industrial Plants](#casal-2018)*. [↩](#fnref-12)
{% endcapture %}

The Britter-McQuaid model would seem to be the perfect fit for screening models, which are often only order of magnitude estimates at best anyways. It gives reasonable concentrations, plausible plume extents, and the explosive mass is ridiculously easy to calculate (slightly more tedious if you are using the 2/3 cut-off region but nothing that couldn't be worked out in advance if this was going to be incorporated into a routine calculation tool).

<div class="notice">
  {{ footnote-11 | markdownify }}

  {{ footnote-12 | markdownify }}
</div>


For a complete listing of code used to generate data and figures, please see the [corresponding julia notebook](https://github.com/aefarrell/aefarrell.github.io/blob/main/_notebooks/2023-03-12-Britter-McQuaid.ipynb)
{: .notice--info}


## References

+ <a name="ccps-1999">AIChE/CCPS</a>. 1999. *Guidelines for Consequence Analysis of Chemical Releases.* New York: American Institute of Chemical Engineers
+ <a name="bakkum-2005">Bakkum</a>, E.A. and N.J. Duijm. 2005. "Chapter 4 - Vapour Cloud Dispersion" in *Methods for the Calculation of Physical Effects, CPR 14E* (TNO Yellow Book) Edited by C.J.H. van den Bosch and R.A.P.M. Weterings. The Netherlands. [url](https://repository.tno.nl/islandora/object/uuid:4928209c-5998-4261-9393-3d55073e6e87)
+ <a name="britter-1988">Britter</a>, Rex E. and J. McQuaid. 1988. *Workbook on the Dispersion of Dense Gases. HSE Contract Research Report No. 17/1988*
+ <a name="casal-2018">Casal</a>, Joachim. 2018. *Evaluation of the Effects of Consequences of Major Accidents in Industrial Plants, 2nd Ed.* Amsterdam: Elsevier
+ <a name="lees-1996">Lees</a>, Frank P. 1996. *Loss Prevention in the Process Industries, 2nd ed.* Oxford: Butterworth-Heinemann
+ <a name="woodward-1998">Woodward</a>, John L. 1998. *Estimating the Flammable Mass of a Vapour Cloud*, New York: American Institute of Chemical Engineers