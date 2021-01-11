#include "interpolation.h"


std::vector<Pnt3> interpolate(const std::vector<Pnt3>& input, 
    double time_interval, 
    double dist_constraint2)
{
    std::vector<Pnt3> ret;
    if(input.size() <= 1){
        return ret;
    }

    std::vector<Pnt3> input_copy(input);
    std::sort(input_copy.begin(), input_copy.end());
    double mint = input_copy[0].t;
    double maxt = input_copy.back().t;

    double startt = ceil(mint/time_interval)*time_interval;
    int idx_before=0;
    int idx_after=1;
    // std::cout << startt << ", " << mint << ", " << maxt << std::endl;
    for(double curt = startt; curt < maxt; curt+= time_interval) {
        while(input_copy[idx_after].t < curt){
            idx_before++;
            idx_after++;
        }
        double t0 = input_copy[idx_before].t;
        double x0 = input_copy[idx_before].x;
        double y0 = input_copy[idx_before].y;
        double t1 = input_copy[idx_after].t;
        double x1 = input_copy[idx_after].x;
        double y1 = input_copy[idx_after].y;
        if(t0 < curt-time_interval || t1 >= curt+time_interval){
            //do not interpolate in this case
            continue;
        }
        double dist2 = (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1);
        if(dist2 >= dist_constraint2){
            continue;
        }
        double curx = ((curt-t0)*x0 + (t1-curt)*x1)/(t1-t0);
        double cury = ((curt-t0)*y0 + (t1-curt)*y1)/(t1-t0); 
        ret.emplace_back(curt, curx, cury);
    }
    return ret;
}
