#include <iostream>
#include <vector>
#include "math.h"
#include <map>
#include <algorithm>
#include "algo_utils.hpp"
#include "misc_utils.hpp"
#include <numeric>
#pragma once

using namespace std;

class elkan_utils{
    template <typename Tfloat>
    inline void find_closest_center(const vector<vector<Tfloat> > &centroids, 
    vector<vector<Tfloat> > &centroid_dists,
    unsigned long long int &dist_counter);

    template <typename Tfloat, typename Tint>
    inline void update_elk_bounds(vector<Tfloat> &upper_bounds, vector<vector<Tfloat> > &lower_bounds,
    vector<Tfloat> &centroid_motion, vector<Tint> &assigned_clusters);
};


template <typename Tfloat>
inline void find_centroid_distances(const vector<vector<Tfloat> > &centroids, 
vector<vector<Tfloat> > &centroid_dists, vector<Tfloat> &closest_center_dist, 
unsigned long long int &dist_counter){

    algorithm_utils alg_utils;
    Tfloat dist = 0;

    for (int i=0; i<centroids.size(); i++){

        centroid_dists[i][i] = std::numeric_limits<float>::max();

        for (int j=i+1; j < centroids.size(); j++){

            dist = alg_utils.calc_euclidean(centroids[i], centroids[j], dist_counter);
            dist /= 2;
            centroid_dists[i][j] = dist;
            centroid_dists[j][i] = dist;

            if (dist < closest_center_dist[i])
                closest_center_dist[i] = dist;
                
            if (dist < closest_center_dist[j])
                closest_center_dist[j] = dist;
        }
    }
}


template <typename Tfloat, typename Tint>
inline void update_elk_bounds(vector<Tfloat> &upper_bounds, vector<vector<Tfloat> > &lower_bounds,
vector<Tfloat> &centroid_motion, vector<Tint> &assigned_clusters){

    // Updating bounds
    for (int i = 0; i < assigned_clusters.size(); i++){
        upper_bounds[i] += centroid_motion[assigned_clusters[i]];

        for (int j = 0; j<centroid_motion.size(); j++){
            lower_bounds[i][j] -= centroid_motion[j];

        }
    }

}
       