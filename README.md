Fast Normal CDF
================

### Introduction

The standard normal CDF
![\\Phi(x)](https://latex.codecogs.com/svg.latex?%5CPhi%28x%29
"\\Phi(x)") is an important function in a broad range of statistical
problems. When we need to evaluate the function many times (for example
in numerical integration), the computation performance may become an
issue.

One way to fast evaluate the function is to use a look-up table, that
is, we pre-compute a set of pairs ![(x\_i,
\\Phi(x\_i))](https://latex.codecogs.com/svg.latex?%28x_i%2C%20%5CPhi%28x_i%29%29
"(x_i, \\Phi(x_i))")\` and then use interpolation to approximate the
function value of a given ![x](https://latex.codecogs.com/svg.latex?x
"x").

This simple library calculates the
![\\Phi(x)](https://latex.codecogs.com/svg.latex?%5CPhi%28x%29
"\\Phi(x)") function using piecewise linear interpolation. The
approximation error is guaranteed to be no greater than ![\\epsilon
= 10^{-7}](https://latex.codecogs.com/svg.latex?%5Cepsilon%20%3D%2010%5E%7B-7%7D
"\\epsilon = 10^{-7}").

### Installation

The package can be installed from Github by calling:

``` r
remotes::install_github("boennecd/fastncdf")
```

### Algorithm

We need to first determine the knots
![x\_i](https://latex.codecogs.com/svg.latex?x_i "x_i") that we want to
pre-compute. Since ![\\Phi(-x) = 1 -
\\Phi(x)](https://latex.codecogs.com/svg.latex?%5CPhi%28-x%29%20%3D%201%20-%20%5CPhi%28x%29
"\\Phi(-x) = 1 - \\Phi(x)"), we only need to consider non-negative
![x\_i](https://latex.codecogs.com/svg.latex?x_i "x_i")’s.

For ![x \> \\Phi^{-1}(1 - \\epsilon)
= 5.1993376](https://latex.codecogs.com/svg.latex?x%20%3E%20%5CPhi%5E%7B-1%7D%281%20-%20%5Cepsilon%29%20%3D%205.1993376
"x \> \\Phi^{-1}(1 - \\epsilon) = 5.1993376"), we set ![\\Phi(x)
= 1](https://latex.codecogs.com/svg.latex?%5CPhi%28x%29%20%3D%201
"\\Phi(x) = 1") and hence the error is bounded by
![\\epsilon](https://latex.codecogs.com/svg.latex?%5Cepsilon
"\\epsilon"). Let ![x\_0
= 0](https://latex.codecogs.com/svg.latex?x_0%20%3D%200 "x_0 = 0"),
![x\_i = ih](https://latex.codecogs.com/svg.latex?x_i%20%3D%20ih
"x_i = ih"), ![i = 0, 1, ...,
N](https://latex.codecogs.com/svg.latex?i%20%3D%200%2C%201%2C%20...%2C%20N
"i = 0, 1, ..., N"), where ![N](https://latex.codecogs.com/svg.latex?N
"N") is the smallest integer such that ![N h
\> 5.1993376](https://latex.codecogs.com/svg.latex?N%20h%20%3E%205.1993376
"N h \> 5.1993376"). Then we need to determine the interval width
![h](https://latex.codecogs.com/svg.latex?h "h") to satisfy the error
bound.

For piecewise linear interpolation, the error is bounded by

  
![E(t) \\leq 1/8 \\cdot \\lVert
f''\\rVert\_{\\infty}h^2](https://latex.codecogs.com/svg.latex?E%28t%29%20%5Cleq%201%2F8%20%5Ccdot%20%5ClVert%20f%27%27%5CrVert_%7B%5Cinfty%7Dh%5E2
"E(t) \\leq 1/8 \\cdot \\lVert f''\\rVert_{\\infty}h^2")  

(Source
<http://pages.cs.wisc.edu/~amos/412/lecture-notes/lecture09.pdf>)

Since ![\\Phi''(x) = \\phi'(x) = -x
\\phi(x)](https://latex.codecogs.com/svg.latex?%5CPhi%27%27%28x%29%20%3D%20%5Cphi%27%28x%29%20%3D%20-x%20%5Cphi%28x%29
"\\Phi''(x) = \\phi'(x) = -x \\phi(x)"), it can be shown that ![\\lVert
\\Phi''\\rVert\_\\infty = \\psi(1)
= 0.2419707](https://latex.codecogs.com/svg.latex?%5ClVert%20%5CPhi%27%27%5CrVert_%5Cinfty%20%3D%20%5Cpsi%281%29%20%3D%200.2419707
"\\lVert \\Phi''\\rVert_\\infty = \\psi(1) = 0.2419707").

Therefore ![h](https://latex.codecogs.com/svg.latex?h "h") can be
calculated as:

``` r
(h <- sqrt(8 / dnorm(1) * 1e-7))
```

    ## [1] 0.001818292

So the ![x\_i](https://latex.codecogs.com/svg.latex?x_i "x_i")s and
![y\_i](https://latex.codecogs.com/svg.latex?y_i "y_i")s values are:

``` r
x <- seq(0, qnorm(1 - 1e-7) + h, by = h)
length(x)
```

    ## [1] 2861

``` r
y <- pnorm(x)
```

We can then call `dput(x)` and `dput(y)` to get the data we need.

### Performance

We compare the speed of `fastpnorm()` and `fastpnorm_preallocated()`
with `pnorm()` in R:

``` r
library(fastncdf)
x <- seq(-6, 6, by = 1e-6)
system.time(y <- pnorm(x))
```

    ##    user  system elapsed 
    ##   0.502   0.012   0.515

``` r
system.time(fasty <- fastpnorm(x))
```

    ##    user  system elapsed 
    ##   0.034   0.023   0.057

``` r
system.time(fasty_prec <- fastpnorm(x, TRUE))
```

    ##    user  system elapsed 
    ##   0.124   0.016   0.139

``` r
range(y - fasty)
```

    ## [1] -9.99999e-08  9.99999e-08

``` r
range(y - fasty_prec)
```

    ## [1] -9.99999e-08  9.99999e-08

``` r
# if we already had a vector with values then we can use a faster version
res <- rep(0., length(x))
system.time(fastpnorm_preallocated(x, res))
```

    ##    user  system elapsed 
    ##   0.028   0.000   0.028

``` r
all.equal(res, fasty)
```

    ## [1] TRUE

We plot the error versus the quantile below:

``` r
par(mar = c(5, 5, 1, 1))
xs <- seq(-6, 6, length.out = 2000)
plot(xs, fastpnorm(xs) - pnorm(xs), type = "h",
     bty = "l", xlab = expression(x), ylab = "Error")
abline(h = 0, lty = 2)
```

![](man/README_files/err_plt-1.png)<!-- -->

``` r
plot(xs, fastpnorm(xs, TRUE) - pnorm(xs), type = "h",
     bty = "l", xlab = expression(x), ylab = "Error")
abline(h = 0, lty = 2)
```

![](man/README_files/err_plt-2.png)<!-- -->
