## Fast Normal CDF

### Introduction

The standard normal CDF `Φ(x)` is an important function in a broad range of
statistical problems. When we need to evaluate the function many times
(for example in numerical integration), the computation performance may become
an issue.

One way to fast evaluate the function is to use a look-up table, that is,
we pre-compute a set of pairs `(x[i], Φ(x[i]))` and then use interpolation
to approximate the function value of a given `x`.

This simple library calculates the `Φ(x)` function using piecewise linear
interpolation. The approximation error is guaranteed to be no greater than
`ε = 1e-7`.

### Algorithm

We need to first determine the knots `x[i]` that we want to pre-compute.
Since `Φ(-x) = 1 - Φ(x)`, we only need to consider non-negative `x[i]`'s.

For `x > Φ^(-1)(1 - 1e-7) = 5.199338`, we set `Φ(x) = 1` and hence the error is
bounded by `ε`. Let `x[0] = 0`, `x[i] = i * h`, `i = 0, 1, ..., N`, where
`N` is the smallest integer such that `N * h > 5.199338`.
Then we need to determine the interval width `h` to satisfy the error bound.

For piecewise linear interpolation, the error is bounded by

`E(t) ≤ 1/8 * ||f''||∞ * h^2`

(Source [http://pages.cs.wisc.edu/~amos/412/lecture-notes/lecture09.pdf](http://pages.cs.wisc.edu/~amos/412/lecture-notes/lecture09.pdf))


Since `Φ''(x) = φ'(x) = -x * φ(x)`, it can be shown that `||Φ''||∞ = φ(1) = 0.2419707`.

Therefore `h` can be calculated as


```r
h = sqrt(8 / dnorm(1) * 1e-7)
h
```

```
## [1] 0.001818292
```

So the `x` and `y` values are


```r
x = seq(0, qnorm(1 - 1e-7) + h, by = h)
length(x)
```

```
## [1] 2861
```

```r
y = pnorm(x)
```

We write the data to a header file `fastncdf_data.h`:


```r
op = options(digits = 15)
f = "src/fastncdf_data.h"
wrt = function(...) cat(..., "\n", sep = "", file = f, append = TRUE)

cat("static const double fastncdf_max = ", x[length(x)], ";\n", sep = "", file = f)
wrt("static const double fastncdf_h = ", h, ";")
wrt("static const double fastncdf_x [] = {")
con = textConnection("xdata", "w")
write(x, file = con, ncolumns = 5, sep = ", ")
close(con)
wrt(paste(xdata, collapse = ",\n"))
wrt("};")
wrt("static const double fastncdf_y [] = {")
con = textConnection("ydata", "w")
write(y, file = con, ncolumns = 5, sep = ", ")
close(con)
wrt(paste(ydata, collapse = ",\n"))
wrt("};")

options(op)
```

### Performance

We compare the speed of `fastncdf()` and `pnorm()` in R.


```r
library(Rcpp)
sourceCpp("test.cpp")

x = seq(-6, 6, by = 1e-6)
system.time(y <- pnorm(x))
```

```
##    user  system elapsed
##   1.032   0.028   1.058
```

```r
system.time(fasty <- fastncdf(x))
```

```
##    user  system elapsed
##   0.124   0.022   0.145
```

```r
max(abs(y - fasty))
```

```
## [1] 9.99999e-08
```
