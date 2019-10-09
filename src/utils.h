#include <RcppArmadillo.h>
#include "stephens.h"

double update_alpha(double, double, double, int, int);
void print_clusters(std::vector < std::vector < int > >);