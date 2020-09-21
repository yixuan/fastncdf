#ifndef FASTNCDF_H
#define FASTNCDF_H
#include <cstddef>

/** computes a pnorm approximation for one point. */
double fastncdf(double const x, bool const precise_tail) noexcept;

/** computes a pnorm approximation for more points. */
void fastncdf(double const *x, double * res, size_t const n,
              bool const precise_tail)
  noexcept;

#endif /* FASTNCDF_H */
