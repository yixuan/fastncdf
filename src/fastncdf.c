#include "fastncdf.h"

double fastncdf_pos(double x)
{
    #include "fastncdf_data.h"

    if(x >= fastncdf_max)  return 1.0;

    const int i = (int)(x / fastncdf_h);
    const double w = (x - fastncdf_x[i]) / fastncdf_h;
    return w * fastncdf_y[i + 1] + (1.0 - w) * fastncdf_y[i];
}

double fastncdf(double x)
{
    if(x < 0)
        return 1.0 - fastncdf_pos(-x);

    return fastncdf_pos(x);
}
