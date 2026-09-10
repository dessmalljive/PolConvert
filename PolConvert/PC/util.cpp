#include <math.h>
#include "util.h"


bool areSame(double a, double b)
{
	return fabs(a-b) < EPSILON;
}
