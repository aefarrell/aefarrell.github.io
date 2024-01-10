---
title: "Dynamic Mode Decomposition"
last_modified_at: 2023-12-26
toc: true
toc_label: "Contents"
toc_sticky: true
comments: true
categories:
  - notes
tags:
  - dynamical systems
  - system identification
tagline: "Dynamic mode decomposition of fluid flow problems"
header:
  overlay_image: /images/chuttersnap-unsplash-header.jpg
  overlay_filter: 0.25
  caption: "Photo by [CHUTTERSNAP](https://unsplash.com/@chuttersnap?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText) on [Unsplash](https://unsplash.com/)"
---

# Dynamic Mode Decomposition

Recently I've been playing around with Dynamic Mode Decomposition (DMD) and this notebook compiles my notes and julia code in one place for later reference.

Very generally DMD is an approach to system identification problems that is well suited for high dimensional data and systems with coherent spatio-temporal structures. In particular DMD finds the "best fit" linear approximation to the dynamical system, i.e. it finds the matrix **A** such that

$$ \mathbf{\dot{x} } = \mathbf{A x} $$

Where **x** is the high dimensional state vector for the system. One key strength of DMD is that it allows one to calculate **x**(t) without explicitly calculating **A**. This may not seem like a particularly useful property on its face unless one notes that the matrix **A** is *n*&times;*n* and, for systems with a very large *n* (i.e. very high dimensionality) that can be *huge*. Context for *huge* is also important: a matrix that fits easily in memory on my laptop may be infeasibly *huge* for an embedded system. For control applications, such as MPC, DMD may be a good method for generating approximations that are both *good* and *space efficient*.

## Example: Flow Past a Cylinder

As a motivating example, I am going to use the flow past a cylinder dataset from [Data-Driven Science and Engineering](http://databookuw.com), specifically the [matlab dataset](http://databookuw.com/DATA.zip). This dataset is the simulated vorticity for fluid flow past a cylinder. The vector **x** in this case is the vorticity at every point in the discretized flow field at a particular time; a two dimensional array of 89,351 pixels reshaped into a column vector. The data is a sequence of equally spaced snapshots of the flow field, and ultimately we wish to generate a linear system that best approximates this.

The `MAT` package allows us to import data from matlab data files directly into julia


```julia
using MAT

file = matopen("data/CYLINDER_ALL.mat")

# import the data set
data = read(file, "VORTALL");

# the orinal dimensions of each snapshot
nx = Int(read(file, "nx"))
ny = Int(read(file, "ny"))

# the final dimensions of the data matrix
n, m = size(data)
```

  (89351, 151)


The data set, `data`, has already been processed into the form we need: each column represents a "frame" of the animation. We can walk through the matrix, taking each column and re-shaping it back into a 2D array, and recover the original flow as a movie.

![gif](/images/dynamic_mode_decomposition_files/output_2_0.gif)

The data set has the property that the number of data points at each time step, *n*, is much greater than the number of time steps, *m*. In fact *n* is large enough that the *n*&times;*n* matrix **A** might be unwieldy to store: If we assume it is a dense matrix of 64-bit floats, 8 bytes each, we would need ~64GB of memory just to store it.


```julia
size_A_naive = n*n*8
```

  63868809608


## Exact DMD

DMD provides us a method to both find a *best fit* approximation for **A** while also being more space (and computation) efficient. To get there we first need to define what a *best fit* means.

### Best Fit Matrix

Consider the general linear system **Y** = **AX**, where **Y** is a *n* &times; *m* matrix of outputs, **X** is a *n* &times; *m* matrix of inputs and **A** is an *n* &times; *n* linear transformation matrix. We say that the *best fit* matrix **A** is the matrix that minimizes

$$ \| \mathbf{ A X } - \mathbf{Y} \|_{F} $$

where $ \| \cdots \|_{F}$ is the [Frobenius norm](https://mathworld.wolfram.com/FrobeniusNorm.html).

The solution to which is 

$$ \mathbf{A} = \mathbf{YX}^{\dagger} $$

where **X<sup>&dagger;</sup>** is the [Moore-Penrose pseudoinverse](https://mathworld.wolfram.com/Moore-PenroseMatrixInverse.html) of **X**.<a href="#fn-1" class="sidenote-number"></a><span class="sidenote" id="fn-1">I think this can be shown fairly easily by starting with the definition of the Frobenius norm $ \| \mathbf{ A X } - \mathbf{Y} \|_{F}^{2} = \mathrm{Tr}\left( \left(\mathbf{ A X } - \mathbf{Y}\right)\left(\mathbf{ A X } - \mathbf{Y} \right)^{T} \right) $ and finding the matrix **A** that minimizes that using standard [matrix calculus](https://en.wikipedia.org/wiki/Matrix_calculus), and some properties of the pseudoinverse.</span>



### Singular Value Decomposition

The conventional way of calculating the Moore-Penrose pseudoinverse is to use the [Singular Value Decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition): for a matrix **X** with SVD **X**=**U&Sigma;V<sup>*</sup>**, the pseudoinverse is **X<sup>&dagger;</sup>**=**V&Sigma;<sup>-1</sup>U<sup>*</sup>**. Returning to the best fit matrix **A** we find

$$ \mathbf{A} = \mathbf{Y} \mathbf{V} \mathbf{\Sigma}^{-1} \mathbf{U}^{*} $$

We can calculate a projection of **A** onto the space of the upper singular vectors **U**

$$ \tilde{ \mathbf{A} } = \mathbf{U}^{*} \mathbf{A} \mathbf{U} = \mathbf{U}^{*} \mathbf{Y} \mathbf{V} \mathbf{\Sigma}^{-1} \mathbf{U}^{*} \mathbf{U} = \mathbf{U}^{*} \mathbf{Y} \mathbf{V} \mathbf{\Sigma}^{-1} $$

Which then allows us to *reconstruct* the matrix **A** on demand while only needing to store the matrices **&Atilde;** and **U**, by the following
$$ \mathbf{A} = \mathbf{U} \tilde{ \mathbf{A} } \mathbf{U}^{*} $$

This is useful when *n* &gt; &gt; *m* as **U** is *n*&times;*m* and **&Atilde;** is *m*&times;*m*. For this example this has reduced the memory requirement to ~108MB, a &gt;99.8% reduction


```julia
size_A_exact = (n*m + m*m)*8
```

  108118416


```julia
size_A_exact/size_A_naive
```

  0.0016928202774340957


Returning to the original problem, we have a sequence of discrete snapshots arranged in a matrix such that each column, *k*, is the vector **x**<sub>k</sub>. Our aim, then, is to find the *best fit* matrix **A** for the linear system

$$ \mathbf{x}_{k+1} = \mathbf{A} \mathbf{x}_k $$

for all **x**<sub>k</sub> in our data set. Or in other words, to find the *best fit* matrix **A** for the system

$$ \mathbf{X}_{2} = \mathbf{A} \mathbf{X}_{1} $$

where **X**<sub>1</sub> is the matrix of all of the vectors **x**<sub>k</sub> and **X**<sub>2</sub> is the matrix of the corresponding **x**<sub>k+1</sub>'s.

Though, using DMD, we will instead calculate **&Atilde;** and **U**, leaving us with

$$ \mathbf{x}_{k+1} = \mathbf{U} \mathbf{ \tilde{A} } \mathbf{U}^{*} \mathbf{x}_k $$

To start, we divide the data set into **X**<sub>1</sub> and **X**<sub>2</sub>


```julia
using LinearAlgebra
```


```julia
# dividing into past and future states
X₁ = data[:, 1:end-1];
X₂ = data[:, 2:end];
```

Then compute the SVD of **X**<sub>1</sub>.<a href="#fn-2" class="sidenote-number"></a><span class="sidenote" id="fn-2">The `svd` function in julia returns the singular values in a `Vector`, but for later on it will be more convenient have this as a `Diagonal` matrix.</span>


```julia
# SVD
U, Σ, V = svd(X₁)
Σ = Diagonal(Σ);
```

Then calculate the projection **&Atilde;** (I am pre-computing **YV&Sigma;**<sup>-1</sup> as that will come in handy later)


```julia
# projection
YVΣ⁻¹ = X₂*V*Σ^-1
Ã = U'*YVΣ⁻¹

size(Ã)
```

  (150, 150)


We can then calculate the predicted **x**<sub>k+1</sub>'s, without ever having to actually compute (or store) **A**


```julia
X̂₂_exact = (U*(Ã*(U'*X₁)));
```

As before, we can step through the matrix, extract each frame of the 2D flow field, and animate them, giving us a general sense of how well this worked

![gif](/images/dynamic_mode_decomposition_files/output_11_0.gif)



### Dynamic Modes

Of course this only solves the problem in the discrete case (for control applications that may be all you need). Consider again the system $ \mathbf{\dot{x} } = \mathbf{A x} $, the solution to this differential equation is

$$ \mathbf{x}\left( t \right) = e^{\mathbf{A}t} \mathbf{x}_{0} $$

where **x**<sub>0</sub> is the initial conditions. If the matrix **A** has [eigendecomposition](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix) **&Phi;&Lambda;&Phi;**<sup>-1</sup> then this can be written as

$$ \mathbf{x}\left( t \right) = \mathbf{\Phi} e^{\mathbf{\Lambda}t} \mathbf{\Phi}^{-1} \mathbf{x}_{0} $$

So it would be very convenient if we could get those eigenvalues and eigenvectors, preferably without having to actually compute **A**. 

Recall, by definition, the projection matrix **&Atilde;** is [unitarily similar](https://en.wikipedia.org/wiki/Matrix_similarity) to **A**, which means the eigenvalues are identical. The eigenvectors of **A** can also be recovered from properties of **&Atilde;**: Suppose **&Atilde;** has the eigendecomposition **W&Lambda;W**<sup>-1</sup>

$$ \mathbf{ \tilde{A} } \mathbf{W} = \mathbf{W} \mathbf{\Lambda} \\
\mathbf{U}^{*} \mathbf{A} \mathbf{U} \mathbf{W} = \mathbf{W} \mathbf{\Lambda} \\
\mathbf{U} \mathbf{U}^{*} \mathbf{A} \mathbf{U} \mathbf{W} = \mathbf{U} \mathbf{W} \mathbf{\Lambda} \\
\mathbf{A} \mathbf{\Phi} = \mathbf{\Phi} \mathbf{\Lambda} $$

where

$$ \mathbf{\Phi} = \mathbf{U} \mathbf{W} $$

This is what is given in the original DMD, however more recent work recommends using

$$ \mathbf{\Phi} = \mathbf{Y} \mathbf{V} \mathbf{\Sigma}^{-1} \mathbf{W} $$


```julia
# calculate eigenvectors and eigenvalues
# of projection Ã
Λ, W = eigen(Ã)
    
# reconstruct eigenvectors of A
Φ = YVΣ⁻¹*W;
```

Whether or not the ultimate goal is to generate the continuous system, the eigenvectors and eigenvalues are useful to examine as they represent the *dynamic modes* of the system.


![svg](/images/dynamic_mode_decomposition_files/output_26_0.svg)

I've played somewhat fast and loose with variables: the **A** for the discrete system is *not the same* **A** as the continuous system. Specifically the eigenvalues of the continuous system, &omega; are related to the eigenvalues of the discrete system, &lambda; by the following

$$ \omega_{i} = {\log{ \lambda_{i} } \over \Delta t} $$

where *&Delta;t* is the time step. The eigenvectors are the same, though. So we can generate a function **x**(t) pretty easily:


```julia
# calculate the eigenvalues for 
# the continuous system
Δt = 1
Ω  = Diagonal(log.(Λ)./Δt)

# precomputing this
Φ⁻¹x₀ = Φ\X₁[:,1]

# continuous system
x̂(t) = real( Φ*exp(Ω .* t)*Φ⁻¹x₀ )
```

![gif](/images/dynamic_mode_decomposition_files/output_15_0.gif)


## Refactoring

Through taking the SVD, the eigenvalue decomposition, and projections, DMD involves generating a whole bunch of matrices, which can be really unwieldy to manage without some structure. The low hanging fruit for refactoring is to introduce a `struct` to store those matrices.


```julia
struct DMD
    r::Integer  # Dimension
    U::Matrix   # Upper Singular Vectors
    Ã::Matrix   # Projection of A
    Λ::Diagonal # Eigenvalues of A
    Φ::Matrix   # Eigenvectors of A
end
```

Then we can introduce a method that takes an input matrix **X** and output matrix **Y** and returns the corresponding `DMD` object. We can take advantage of multiple dispatch to to add further methods, such as for the case where we have a single data matrix **X** and wish to calculate the DMD on the "future" and "past" matrices.


```julia
function DMD(Y::Matrix, X::Matrix)
    # dimension
    r = rank(X)
    
    # Full SVD
    U, Σ, V = svd(X)
    Σ = Diagonal(Σ)
    
    # projection
    YVΣ⁻¹ = Y*V*Σ^-1
    Ã = U'*YVΣ⁻¹
    
    # calculate eigenvectors and eigenvalues
    # of projection Ã
    Λ, W = eigen(Ã)
    Λ = Diagonal(Λ)
    
    # reconstruct eigenvectors of A
    Φ = YVΣ⁻¹*W
    
    return DMD(r,U,Ã,Λ,Φ)
end

function DMD(X::Matrix)
    X₁ = X[:, 1:end-1]
    X₂ = X[:, 2:end]
    return DMD(X₂, X₁)
end
```

We can check that this is doing what it is supposed to be doing by comparing with what we have already done


```julia
d = DMD(data)

# This produces the same result as before
d.Φ == Φ && d.Λ == Diagonal(Λ)
```

  true


If you were to build this into a larger project, it would be worthwhile to define some actual unit tests to validate that the DMD is working properly.

### Discrete System

Since we have a DMD type to work with, we can also refactor how discrete systems are generated. In this case I have defined a `struct` for the discrete system, and then added a method such that any discrete system acts as a callable **x**<sub>k+1</sub>=*f*(**x**<sub>k</sub>)


```julia
struct DiscreteSys
    Ã::Matrix
    U::Matrix
end

function DiscreteSys(d::DMD)
    return DiscreteSys(d.Ã,d.U)
end

function (ds::DiscreteSys)(xₖ)
    return (ds.U*(ds.Ã*(ds.U'*xₖ)))
end
```


```julia
ds = DiscreteSys(d)

# This produces the same result as before
X̂₂_exact == ds(X₁)
```

  true


### Continuous System

Similarly we can refactor the generation of continuous systems, first by defining a `struct` for the continuous system, then by adding a method **x**<sub>t</sub>=*f*(*t*). This requires a little more information: we need to keep track of the initial state of the system **x**<sub>0</sub> as well as the step size &Delta;*t*


```julia
struct ContinuousSys
    Φ⁻¹x₀::Vector
    Ω::Diagonal
    Φ::Matrix
end

function ContinuousSys(d::DMD, x₀, Δt=1)
    Φ⁻¹x₀ = d.Φ\x₀
    Ω = Diagonal(log.(d.Λ.diag)./Δt)
    return ContinuousSys(Φ⁻¹x₀, Ω, d.Φ)
end

function (cs::ContinuousSys)(t)
    return real( cs.Φ*exp(cs.Ω .* t)*cs.Φ⁻¹x₀ )
end
```


```julia
cs = ContinuousSys(d, X₁[:,1]);

# This produces the same result as before
x̂(150) == cs(150)
```

  true


### Large Systems

I have been using the default tools in julia, which work well for small matrices. If you are planning on doing DMD on enormous matrices then it is worth investigating packages such as [IterativeSolvers.jl](https://iterativesolvers.julialinearalgebra.org/stable/), [Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl), [KrylovKit.jl](https://jutho.github.io/KrylovKit.jl/stable/) and others to find better ways than vanilla `svd` and `eigen`. It also may be worth thinking about refactoring the problem to be matrix-free, though that is way beyond the scope of these notes.

## Reduced DMD

Whenever a problem involves computing the SVD of a matrix, dimensionality reduction lurks about in the shadows, winking suggestively. By the [Eckart-Young theorem](https://en.wikipedia.org/wiki/Low-rank_approximation#Proof_of_Eckart–Young–Mirsky_theorem_(for_Frobenius_norm)) we know that the best rank *r* approximation to a matrix **X**=**U&Sigma;V**<sup>T</sup> is the truncated SVD **X**<sub>r</sub>=**U**<sub>r</sub>**&Sigma;**<sub>r</sub>**V**<sub>r</sub><sup>T</sup>, i.e. the SVD truncated to the *r* largest singular values (and corresponding singular vectors). So an obvious step for dimensionality reduction in DMD is substitute a truncated SVD for the full SVD.


```julia
function DMD(Y::Matrix, X::Matrix, r::Integer)   
    # full SVD
    U, Σ, V = svd(X)
    
    # truncating to rank r
    @assert r ≤ rank(X)
    U = U[:, 1:r]
    Σ = Diagonal(Σ[1:r])
    V = V[:, 1:r]
    
    # projection
    YVΣ⁻¹ = Y*V*Σ^-1
    Ã = U'*YVΣ⁻¹
    
    # calculate eigenvectors and eigenvalues
    # of projection Ã
    Λ, W = eigen(Ã)
    Λ = Diagonal(Λ)
    
    # reconstruct eigenvectors of A
    Φ = YVΣ⁻¹*W
    
    return DMD(r,U,Ã,Λ,Φ)
end

function DMD(X::Matrix, r::Integer)
    X₁ = X[:, 1:end-1]
    X₂ = X[:, 2:end]
    return DMD(X₂, X₁, r)
end
```

One consequence of truncation, however, is that the resulting matrix **U**<sub>r</sub> is only semi-unitary, in particular

$$ \mathbf{U}_{r}^{*} \mathbf{U}_{r} = \mathbf{I}_{r \times r} $$

but

$$ \mathbf{U}_{r} \mathbf{U}_{r}^{*} \ne \mathbf{I}_{n \times n} $$

This leads to a complication as the matrix **U** is required to be unitary, in particular when recovering **A** from the projection matrix **&Atilde;**, and also when recovering the eigenvalues and eigenvectors of **A** from **&Atilde;**.

But, supposing that this at least approximately works, we are still left with the problem of picking an appropriate value for *r*. One could look at the singular values and pick one based on structure. For this problem it looks like an elbow happens at *r=45*.


![svg](/images/dynamic_mode_decomposition_files/output_47_0.svg)

We can then generate a set of predictions for the reduced DMD, with *r=45*, and compare with the exact DMD


```julia
ds_45 = DiscreteSys(DMD(data, 45))
X̂₂_45 = ds_45(X₁)

norm(X₂ - X̂₂_45) # Frobenius norm
```

  0.005459307491383062



```julia
norm(X₂ - X̂₂_exact)
```

  0.0005597047465277092


An alternative is to specify how much of the variance in the original data set needs to be captured. The singular values are a measure of the variance in the data, and so keeping the top *p* percent of the total variance equates to keeping the top *p* percent of the sum of all of the singular values.

That is to say we calculate the *r* such that

$$ { {\sum_{i}^{r} \sigma_i} \over {\sum_{i}^{m} \sigma_i} } \le p  $$

where &sigma;<sub>i</sub> is the *i*th singular value (in order of largest to smallest).


```julia
function DMD(Y::Matrix, X::Matrix, p::AbstractFloat)
    @assert p>0 && p≤1
    
    # full SVD
    U, Σ, V = svd(X)
    
    # determine required rank
    r = minimum( findall( >(p), cumsum(Σ)./sum(Σ)) )
    
    # truncate
    @assert r ≤ rank(X)
    U = U[:, 1:r]
    Σ = Diagonal(Σ[1:r])
    V = V[:, 1:r]
    
    # projection
    YVΣ⁻¹ = Y*V*Σ^-1
    Ã = U'*YVΣ⁻¹
    
    # calculate eigenvectors and eigenvalues
    # of projection Ã
    Λ, W = eigen(Ã)
    Λ = Diagonal(Λ)
    
    # reconstruct eigenvectors of A
    Φ = YVΣ⁻¹*W
    
    return DMD(r,U,Ã,Λ,Φ)
end

function DMD(X::Matrix, p::AbstractFloat)
    X₁ = X[:, 1:end-1]
    X₂ = X[:, 2:end]
    return DMD(X₂, X₁, p)
end
```

Capturing 99% of the variance, in this case, requires only keeping the first 14 singular values.

![svg](/images/dynamic_mode_decomposition_files/output_54_0.svg)


![gif](/images/dynamic_mode_decomposition_files/output_29_0.gif)

There are also methods for finding the optimal rank for truncated SVD for a data set that involves gaussian noise which I am not going to go into here.

So, supposing that *p*=0.99 works for us, how much further have we reduced the size of our matrices?


```julia
# for p=0.99, r=14
r = 14
size_A_reduced = (n*r + r*r)*8
```

  10008880


To recover the (approximate) **A** matrix we only need to store 10MB, a ~91% reduction over the exact DMD


```julia
size_A_reduced/size_A_exact
```

  0.09257331331972159


and a >99.98% reduction of the naive case (recall the naive approach of storing the entire **A** matrix would take ~64GB)


```julia
size_A_reduced/size_A_naive
```

  0.00015670998193688458

### Truncated SVD and Large Systems

In the above code I simply calculated the full SVD and then truncated it after the fact. If *m* (the rank of **X**) is particularly large, then this can be hilariously inefficient. In those cases it may be worth writing a method that uses [TSVD.jl](https://tsvd.julialinearalgebra.org/stable/) to efficiently calculate *only* the first *r* singular values -- as opposed to calculating all *m* singular values and then chucking out most of them.

## Compressed DMD

Compressed DMD attempts to tackle the *slowest step* in the DMD algorithm: calculating the SVD. An SVD on full data is $\mathcal{O}\left( n m^2 \right)$ if we instead *compress* the data from *n* dimensions to *k* dimensions then the cost of the SVD is reduced to either $\mathcal{O}\left( k m^2 \right)$ (when *k*&gt;*m*) or $\mathcal{O}\left( m k^2 \right)$ (when *k* &lt; *m*), which for large *n* can be a dramatic speed-up.

Suppose we have some *k*&times;*n* unitary matrix **C** which compresses our input matrix **X** into the *compressed* input matrix **X**<sub>c</sub> and our output matrix **Y** into the *compressed* output matrix **Y**<sub>c</sub>

$$ \mathbf{X}_c = \mathbf{C} \mathbf{X} \\ \mathbf{Y}_c = \mathbf{C} \mathbf{Y} $$

We suppose again that **X** has the SVD **X**=**U&Sigma;V<sup>*</sup>**, then

$$ \mathbf{X}_c = \mathbf{C} \mathbf{X} = \mathbf{C} \mathbf{U} \mathbf{\Sigma} \mathbf{V}^{*} $$

and, since **C** is unitary, the SVD of **X**<sub>c</sub> is

$$ \mathbf{X}_c = \mathbf{U}_c \mathbf{\Sigma} \mathbf{V}^{*} $$

where **U**<sub>c</sub>=**CU** is the upper singular values of the compressed input matrix.

The projection matrix **&Atilde;**<sub>c</sub> of the compressed input matrix is

$$ \mathbf{ \tilde{A} }_c = \mathbf{U}^{*}_c \mathbf{Y}_c \mathbf{V} \mathbf{\Sigma}^{-1} \\
= \left( \mathbf{CU} \right)^{*} \mathbf{C} \mathbf{Y} \mathbf{V} \mathbf{\Sigma}^{-1} \\
= \mathbf{U}^{*} \mathbf{C}^{*} \mathbf{C} \mathbf{Y} \mathbf{V} \mathbf{\Sigma}^{-1} \\
= \mathbf{U}^{*} \mathbf{C}^{*} \mathbf{C} \mathbf{Y} \mathbf{V} \mathbf{\Sigma}^{-1} \\
= \mathbf{U}^{*} \mathbf{Y} \mathbf{V} \mathbf{\Sigma}^{-1} \\
= \mathbf{ \tilde{A} } $$

and so we should recover the same eigenvalues and eigenvectors as from the *uncompressed* data.


```julia
using SparseArrays

function cDMD(Y::Matrix, X::Matrix, C::AbstractSparseMatrix)   
    # determining dimensionality
    r = rank(X)
       
    # compress the X and Y
    Xc = C*X
    Yc = C*Y
    
    # singular value decomposition
    Uc, Σc, Vc = svd(Xc)
    Σc = Diagonal(Σc)
       
    # projection
    Ã = Uc'*Yc*Vc*inv(Σc)
    U = C'*Uc
    
    # calculate eigenvectors and eigenvalues
    # of projection Ã
    Λ, W = eigen(Ã)
    Λ = Diagonal(Λ)
    
    # reconstruct eigenvectors of A
    Φ = Y*Vc*inv(Σc)*W
    
    return DMD(r,U,Ã,Λ,Φ)
end

function cDMD(X::Matrix, C::AbstractSparseMatrix)
    X₁ = X[:, 1:end-1]
    X₂ = X[:, 2:end]
    return cDMD(X₂, X₁, C)
end
```

The giant caveat is: how do we generate a *unitary* compression matrix? In fact we can relax this condition if we simply want to recover the eigenvalues and eigenvectors of **A**. It is enough that the data is sparse in some basis and that the compression matrix is incoherent with respect to that basis.

We can think of **C** as a set of *k* (1&times;*n*)-row vectors that project an *n* dimensional vector **x** onto a *k* dimensional space. There are several ways of finding the basis for this projection -- e.g. a uniform random projection or a gaussian projection -- but by far the simplest is to pick a random subset of *k* single pixels and only take the measurements from those pixels.


```julia
function cDMD(Y::Matrix, X::Matrix, k::Integer)
    n, m = size(X)
    @assert k≤n
    
    # build (sparse) compression matrix
    C = spzeros(k, n)
    for i in 1:k
        C[i,rand(1:n)] = 1
    end

    return cDMD(Y, X, C)
end
   
function cDMD(X::Matrix, k::Integer)
    X₁ = X[:, 1:end-1]
    X₂ = X[:, 2:end]
    return cDMD(X₂, X₁, k)
end
```

Suppose we sample at 300 randomly chosen points in the flow field to form the compression matrix


```julia
k = 300

C = spzeros(k, n)
for i in 1:k
    C[i,rand(1:n)] = 1
end
```

That is to say we are only sampling the vorticity at the green dots. This reduces the dimensionality of the data *going in* to the DMD algorithm from 89351 to 300.

![gif](/images/dynamic_mode_decomposition_files/output_36_0.gif)

We can generate a few different compressed DMDs to get a sense of how this impacts the overall performance (in terms of the Frobenius norm) and, much like we saw with reduced DMD, there are diminishing returns,


![svg](/images/dynamic_mode_decomposition_files/output_73_0.svg)

Using the compression matrix from above, we can generate a compressed DMD<a href="#fn-3" class="sidenote-number"></a><span class="sidenote" id="fn-3">While we can reconstruct the eigenvalues and eigenvectors quite successfully, I don't believe we adequately reconstruct **U**, and so this really only works for the *continuous* system. The reconstruction of **U** strongly depends on **C** being unitary and I don't think that condition can be relaxed.</span>


![gif](/images/dynamic_mode_decomposition_files/output_38_0.gif)


The compressed DMD does not actually reduce the storage size of any of the matrices, it is more a technique to speed up the calculation of the SVD. Compressed DMD and reduced DMD can be combined: first by compressing the *n*&times;*m* matrix **X** to a *k*&times;*m* matrix **X**<sub>c</sub> and then finding the best rank *r* approximation to the compressed matrix by truncating the SVD to the *r* largest singular values. The reduction step reduces the memory requirements and, if truncated SVD is used as well, this could significantly improve performance for enormous systems.

There is a related approach called *compressed sensing* DMD, in which the full state vector is not available in the first place. A much smaller dimension set of measurements is sampled and the full state DMD generated using the same general idea as compressed DMD. It isn't that much of a leap from what is above, just with a convex optimization step added to reconstruct the actual state matrix for a given set of measurements.


## Physics Informed DMD

The idea behind physics informed DMD is that the physics of the system imposes *structure* upon the solution, which we can build into the DMD algorithm. This way we generate results that are consistent with physical reality. Which is to say that we are not merely finding the best fit matrix **A**, we are finding the best fit matrix **A** *subject to* some constraints on its structure. [The paper](https://arxiv.org/pdf/2112.04307.pdf) I am using as a reference gives a nice table of different types of flow problems and the sort of structure one might want to impose upon the solution

![image.png](/images/dynamic_mode_decomposition_files/att1.png)

Conveniently the flow past a cylinder example is on that table (that definitely wasn't a motivating factor for choosing it as the example in the first place, nope, not at all) and what we want to impose on the solution is conservation of energy. Conservation of energy in this case equates to requiring that **A** be unitary, which is the standard [procrustes problem](https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem)

We modify the best fit such that we are looking for the **A** matrix that minimizes

$$ \| \mathbf{ A X } - \mathbf{Y} \|_{F} \\
\textrm{ subject to } \mathbf{A}^{*} \mathbf{A} = \mathbf{I} $$

For which the standard solution is to define a matrix **M**

$$ \mathbf{M} = \mathbf{Y} \mathbf{X}^{*} $$

supposing **M** has SVD

$$ \mathbf{M} = \mathbf{U}_{M} \mathbf{\Sigma}_{M} \mathbf{V}_{M}^{*} $$

then the solution is

$$ \mathbf{A} = \mathbf{U}_{M} \mathbf{V}_{M}^{*} $$

Of course we can't directly compute **M** in many cases for the same reason that we can't directly compute **A** : it would be a *n*&times;*n* matrix and for large *n* that would be enormous. So instead we project **X** and **Y** onto the upper singular values of **X** and solve the procrustes problem in that smaller space:

$$ \mathbf{X} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^{*} $$

$$ \mathbf{ \tilde{X} } = \mathbf{U}^{*} \mathbf{X} $$

$$ \mathbf{ \tilde{Y} } = \mathbf{U}^{*} \mathbf{Y} $$

$$ \mathbf{ \tilde{M} } = \mathbf{ \tilde{Y} } \mathbf{ \tilde{X} }^{*} = \mathbf{U}^{*} \mathbf{Y} \mathbf{X}^{*} \mathbf{U} = \mathbf{U}^{*} \mathbf{M} \mathbf{U}$$

since SVD is invariant to left and right unitary transformations, the SVD of the projected $\mathbf{ \tilde{M} }$ is

$$ \mathbf{ \tilde{M} } = \mathbf{U}_{ \tilde{M} } \mathbf{\Sigma}_{M} \mathbf{V}_{ \tilde{M} }^{*} $$

where 

$$ \mathbf{U}_{ \tilde{M} } = \mathbf{U}^{*} \mathbf{U}_{ M } \textrm{ and } \mathbf{V}_{ \tilde{M} } = \mathbf{U}^{*} \mathbf{V}_{M} $$

and the A matrix which solves the projected procrustes problem is

$$ \mathbf{ \tilde{A} } = \mathbf{U}_{ \tilde{M} } \mathbf{V}_{ \tilde{M} }^{*} = \mathbf{U}^{*} \mathbf{U}_{ M } \mathbf{V}_{M}^{*} \mathbf{U} = \mathbf{U}^{*} \mathbf{A} \mathbf{U} $$

which is *exactly* the projected **A** matrix we need to proceed with reconstructing the eigenvalues and eigenvectors as per the standard DMD algorithm.


```julia
# this is piDMD *only* for the case where A must be unitary
# see arXiv:2112.04307 for details on the alternative cases
function piDMD(Y::Matrix, X::Matrix)
    # dimension
    r = rank(X)
    
    # Full SVD
    U, _, _ = svd(X)
    
    # projection
    Ỹ = U'*Y
    X̃ = U'*X
    M̃ = Ỹ*X̃'
    
    # solve procrustes problem
    Uₘ, _, Vₘ = svd(M̃)
    Ã = Uₘ*Vₘ'
    
    # calculate eigenvectors and eigenvalues
    # of projection Ã
    Λ, W = eigen(Ã)
    Λ = Diagonal(Λ)
    
    # reconstruct eigenvectors of A
    Φ = U*W
    
    return DMD(r,U,Ã,Λ,Φ)
end

function piDMD(X::Matrix)
    X₁ = X[:, 1:end-1]
    X₂ = X[:, 2:end]
    return piDMD(X₂, X₁)
end
```


![gif](/images/dynamic_mode_decomposition_files/output_40_0.gif)

We can compare the Frobenius norm of the actual data versus the predicted, and it's clear the physics informed DMD does not generate as good of a fit as exact DMD. Though it could equally be the case that the exact DMD is over-fitting.


```julia
norm(X₂ - X̂₂_pi, 2)
```

  18.35684111920036


```julia
norm(X₂ - X̂₂_exact, 2)
```

  0.0005597047465277092


The main reason why you would pursue physics informed DMD, though, is not necessarily to generate a better fit as much as to generate better (or more physically realistic) dynamic modes.

Similarly to compressed DMD, physics informed DMD can also be combined with reduced DMD. In this case there are two SVD steps but only the upper singular values of **X**, the **U** matrix, needs to be truncated. The second SVD proceeds without truncation.

For a complete listing of code used to generate data and figures, please see the [corresponding julia notebook](https://github.com/aefarrell/aefarrell.github.io/blob/main/_notebooks/2022-12-18-dynamic_mode_decomposition.ipynb)
{: .notice--info}

## References

+ Baddoo, Peter J., Benjamin Herrmann, Beverley J. McKeon, J. Nathan Kutz, and Steven L. Brunton. "Physics-informed dynamic mode decomposition (piDMD)." (2021) [arXiv:2112.04307](https://arxiv.org/abs/2112.04307) with code available on [github](https://github.com/baddoo/piDMD)
+ Bai, Zhe, Eurika Kaiser, Joshua L. Proctor, J. Nathan Kutz, and Steven L. Brunton. "Dynamic Mode Decomposition for Compressive System Identification." *AIAA Journal*, 58 (2020):561-574 doi:[10.2514/1.J057870](https://doi.org/10.2514/1.J057870)
+ Brunton, Steven L. and J. Nathan Kutz. *[Data Driven Science and Engineering](http://databookuw.com)*. Cambridge: Cambridge University Press, 2019.
>this an excellent resource for more than just the details of DMD (chapter 7). It is more than just a book as well: there are several videos of lectures going through all of the details.

+ Brunton, Steven L., Joshua L. Proctor, and J. Nathan Kutz. "Compressive sampling and dynamic mode decomposition." (2013) [arXiv:1312.5186](https://arxiv.org/abs/1312.5186)
+ Brunton, Steven L., Joshua L. Proctor, Jonathan H. Tu, and J. Nathan Kutz. "Compressed sensing and dynamic mode decomposition." *Journal of Computational Dynamics*, 2 (2015): 165-191. doi: [10.3934/jcd.2015002](https://www.aimsciences.org/article/doi/10.3934/jcd.2015002)
+ Schmid, Peter J. "Dynamic mode decomposition of numerical and experimental data." *Journal of Fluid Mechanics*, 656 (2010):5-28 doi:[10.1017/S0022112010001217](https://doi.org/10.1017/S0022112010001217)
+ Tu, Jonathan H., Clarence W. Rowley, Dirk Martin Luchtenburg, Steven L. Brunton, and J. Nathan Kutz. "On dynamic mode decomposition: Theory and applications." *Journal of Computational Dynamics*, 1 (2014): 391-421. doi:[10.3934/jcd.2014.1.391](https://doi.org/10.3934/jcd.2014.1.391)

