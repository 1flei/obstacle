#pragma once
#include <string> 
#include <memory>
#include <unordered_map>
#include "lsh_density_bks.h"
#include "lsh_density.h"
#include "../util.h"

void build_index(std::unordered_map<int, Traj>& traj_dict, 
    int len, 
    int step, 
    double interval, 
    LSH_KDE_BKS& index);

void build_index(std::unordered_map<int, Traj>& traj_dict, 
    int len, 
    int step, 
    double interval, 
    LSH_KDE& index);