#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "geokm_utils.hpp"
#pragma once

using namespace std;

class annulus_utils{

    public:

    template<typename Tfloat>
    inline Tfloat point_point_inner_product(vector<Tfloat> &point);

    template<typename Tfloat>
    inline Tfloat point_center_inner_product(vector<Tfloat> &point, vector<Tfloat> &center);  

};


template<typename Tfloat>
inline Tfloat annulus_utils::point_point_inner_product(vector<Tfloat> &point){

    Tfloat sum = 0;

    for(int i=0; i< point.size(); i++){
        sum += (point[i] * point[i]);
    }
    return sum;
}

template<typename Tfloat>
inline Tfloat point_center_inner_product(vector<Tfloat> &point, vector<Tfloat> &center){

    float sum = 0;

    for(int i=0; i< point.size(); i++){
        sum += (point[i] * center[i]);
    }
    return sum;
}
