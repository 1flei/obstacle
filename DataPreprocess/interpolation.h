#pragma once

#include <vector>
#include <cmath>
#include "../util.h"
#include <algorithm>
#include <iostream>
#include <cassert>

std::vector<Pnt3> interpolate(const std::vector<Pnt3>& input, double time_interval, double dist_constraint2=1e9);