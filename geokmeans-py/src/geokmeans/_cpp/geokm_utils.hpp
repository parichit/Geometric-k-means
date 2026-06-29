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

class dckm_utils{

    template <typename Tfloat>
    inline void find_midpoints(const vector<Tfloat> &center1, const vector<Tfloat> &center2, 
    vector<vector<vector<Tfloat> > > &midpoint, vector<vector<vector<Tfloat> > > &affine, 
    int &curr_center, int &ot_cen);

    template <typename Tfloat>
    inline bool find_context_direction(const vector<Tfloat> &actual_point, 
    const vector<Tfloat> &centroid_vector, const vector<Tfloat> &midpoint,
    float &vec_sum);

    template <typename TD, typename TI>
    void restore_radius(vector<vector <TD> > &dist_matrix,
    vector<TI> &assigned_clusters, 
    vector<vector <TD> > &cluster_size);

    template <typename TD, typename TI>
    void find_neighbors(vector<vector <TD> > &centroids, 
    vector<vector <TD> > &center_dist_mat, vector<TD> &closest_center,
    vector<vector <TD> > &cluster_size, 
    vector<vector<TI> > &neighbors, vector<vector<vector <TD> > > &mid_points, 
    vector<vector<vector <TD> > > &affine_vectors, 
    unsigned long long int &dist_counter);

    template <typename Tfloat, typename Tint>
    inline void determine_data_expression(vector<vector<Tfloat> > &dataset, 
    vector<vector <Tfloat> > &centroids, vector<vector<Tfloat> > &cluster_size, 
    vector<vector <Tfloat> > &center_dist_mat, vector<Tfloat> &closest_center,
    vector<Tint> &assigned_clusters, 
    vector<vector<Tint> > &neighbors,
    vector<vector<vector <Tfloat> > > &affine_vectors, 
    vector<vector<vector <Tfloat> > > &mid_points, 
    unsigned long long int &dist_counter, vector<Tint> &temp);

    template <typename Tfloat, typename Tint>
    inline void determine_data_expression_record(vector<vector<Tfloat> > &dataset, 
    vector<vector <Tfloat> > &centroids, vector<vector<Tfloat> > &cluster_size, 
    vector<vector <Tfloat> > &center_dist_mat, vector<Tfloat> &closest_center,
    vector<Tint> &assigned_clusters, 
    vector<vector<Tint> > &neighbors,
    vector<vector<vector <Tfloat> > > &affine_vectors, 
    vector<vector<vector <Tfloat> > > &mid_points, 
    unsigned long long int &dist_counter, double &le_exec, 
    double &lhe_exec, double &he_exec, vector<Tint> &temp);

    template <typename Tfloat, typename Tint>
    inline void update_bounds(vector<Tfloat> &upper_bound, vector<Tfloat> &lower_bound,
    vector<Tfloat> &centroid_motion, vector<Tint> &assigned_clusters);

    template <typename Tfloat>
    inline bool find_context_direction_new(const vector<Tfloat> &actual_point, 
    const vector<Tfloat> &centroid_vector, const vector<Tfloat> &midpoint);

    template <typename Tfloat>
    inline void find_closest_center(vector<vector<Tfloat> > &centroids, vector<Tfloat> &closest_center_dist,
    unsigned long long int &dist_counter);

    template <typename TD, typename Tint>
    inline void find_exp_neighbors(vector<vector <TD> > &centroids, 
    vector<vector <TD> > &center_dist_mat, vector<TD> &closest_center,
    vector<vector<pair<TD, Tint>> > &center_center_dist,
    unsigned long long int &dist_counter);

    template <typename TD, typename TI>
    inline void find_neighbors_bounded(vector<vector <TD> > &centroids, 
    vector<vector <TD> > &center_dist_mat, vector<TD> &closest_center,
    vector<vector <TD> > &cluster_size, 
    vector<vector<pair<TD, TI>> > &neighbors, vector<vector<vector <TD> > > &mid_points, 
    vector<vector<vector<TD>> > &affine_vectors,
    vector<vector<TD> > &temp_master, vector<TD> &temp_midpoint, vector<TD> &temp_affine, 
    vector<vector<TD> > &midpoint_holder, vector<vector<TD> > &affine_holder, 
    unsigned long long int &dist_counter);

    template <typename TD, typename TI>
    inline void find_neighbors_p(vector<vector <TD> > &centroids, 
    vector<vector <TD> > &center_dist_mat, vector<TD> &closest_center,
    vector<vector <TD> > &cluster_size, 
    vector<vector<TI> > &neighbors, vector<vector<vector <TD> > > &mid_points, 
    vector<vector<vector <TD> > > &affine_vectors,  
    unsigned long long int &dist_counter, vector<TD> &center_mean);
};


template <typename Tfloat>
inline void find_midpoints(const vector<Tfloat> &center1, const vector<Tfloat> &center2, 
vector<vector<vector<Tfloat> > > &midpoint, vector<vector<vector<Tfloat> > > &affine, 
int &curr_center, int &ot_cen){

    for (int i=0; i<center1.size(); i++){
        midpoint[curr_center][ot_cen][i] = (center1[i] + center2[i])/2;
        // midpoint[ot_cen][curr_center][i] = midpoint[curr_center][ot_cen][i];

        affine[curr_center][ot_cen][i] = center2[i] - midpoint[curr_center][ot_cen][i];
        // affine[ot_cen][curr_center][i] = -affine[curr_center][ot_cen][i];
    }
}


template <typename Tfloat>
inline bool find_context_direction(const vector<Tfloat> &actual_point, 
const vector<Tfloat> &centroid_vector, const vector<Tfloat> &midpoint, float &vec_sum){

    vec_sum = 0.0;
    // float temp = 0;
    
    for (int i=0; i<midpoint.size(); i++){
        // temp = (actual_point[i] - midpoint[i]);
        // temp = midpoint[i] - actual_point[i];
        // if ((temp >0 & centroid_vector[i] <0)|| (temp <0 & centroid_vector[i] > 0))
        //     return false;
            vec_sum += ((actual_point[i] - midpoint[i]) * centroid_vector[i]);
    }

    if (vec_sum>0)
        return true;

    return false;
}




template <typename TD, typename TI>
inline void restore_radius(vector<vector <TD> > &dist_matrix,
vector<TI> &assigned_clusters, 
vector<vector <TD> > &cluster_size){

    for (int i=0; i<cluster_size.size(); i++){
        for (int j=0; j< assigned_clusters.size(); j++){
                if (assigned_clusters[j] == i){
                    if(dist_matrix[j][i] > cluster_size[i][1]){
                            cluster_size[i][1] = dist_matrix[j][i];
                    }
                }
            }
        }
    }



template <typename Tfloat>
inline void find_closest_center(vector<vector<Tfloat> > &centroids, vector<Tfloat> &closest_center_dist,
unsigned long long int &dist_counter){

algorithm_utils alg_utils;
float dist = 0;

for (int i=0; i<closest_center_dist.size(); i++)
    closest_center_dist[i] = std::numeric_limits<float>::max();

for (int i=0; i<centroids.size(); i++){
    for (int j=i+1; j<centroids.size(); j++){
            dist = alg_utils.calc_euclidean(centroids[i], centroids[j], dist_counter);
            dist = dist/2;
            if (dist < closest_center_dist[i])
                closest_center_dist[i] = dist;
            if (dist < closest_center_dist[j])
                closest_center_dist[j] = dist;
        }
    }
}


template <typename TD, typename TI>
inline void find_neighbors(vector<vector <TD> > &centroids, 
vector<vector <TD> > &center_dist_mat, vector<TD> &closest_center,
vector<vector <TD> > &cluster_size, 
vector<vector<TI> > &neighbors, vector<vector<vector <TD> > > &mid_points, 
vector<vector<vector <TD> > > &affine_vectors,  
unsigned long long int &dist_counter){

    TD dist = 0;
    TD radius = 0;
    algorithm_utils alg_utils;
    vector<TI> temp1;

    int curr_center = 0, ot_center = 0, cnt = 0;
    float shortestDist2 = 0;
    
    for (int i=0; i<closest_center.size(); i++)
        closest_center[i] = std::numeric_limits<float>::max();

    
    // Calculate inter-centroid distances
    for(curr_center=0; curr_center<centroids.size(); curr_center++){

        center_dist_mat[curr_center][curr_center] = std::numeric_limits<float>::max();
        
        for (ot_center=curr_center+1; ot_center<centroids.size(); 
        ot_center++){    
                
                dist = alg_utils.calc_euclidean(centroids[curr_center], centroids[ot_center], dist_counter);
                dist /= 2;
                center_dist_mat[curr_center][ot_center] = dist;
                center_dist_mat[ot_center][curr_center] = center_dist_mat[curr_center][ot_center];

                if (dist < closest_center[curr_center])
                    closest_center[curr_center] = dist;
                
                if (dist < closest_center[ot_center])
                    closest_center[ot_center] = dist;

        }

    }

    // Calculate inter-centroid distances
    for(curr_center=0; curr_center<centroids.size(); curr_center++){

        radius = cluster_size[curr_center][1] + closest_center[curr_center];
        cnt = 0;
        
        for (ot_center=0; ot_center<centroids.size(); 
        ot_center++){  

            // Start neighbor finding
            if ((curr_center != ot_center) && 
            (center_dist_mat[curr_center][ot_center] <= radius)){
                                        
                temp1.push_back(ot_center);
           
                // Get the mid-point coordinates for this pair of centroids
                find_midpoints(centroids[curr_center], centroids[ot_center], mid_points, affine_vectors, curr_center, ot_center);
                cnt++;
            }
        
        }   

            if (cnt>=1){
                neighbors[curr_center] = temp1;
            }
 
            else if(cnt == 0){
                temp1.push_back(-100);
                neighbors[curr_center] = temp1;
            }

            cluster_size[curr_center][2] = cnt;
            temp1.clear();
    }
}


template <typename TD, typename TI>
inline void find_neighbors_p(vector<vector <TD> > &centroids, 
vector<vector <TD> > &center_dist_mat, vector<TD> &closest_center,
vector<vector <TD> > &cluster_size, 
vector<vector<TI> > &neighbors, vector<vector<vector <TD> > > &mid_points, 
vector<vector<vector <TD> > > &affine_vectors,  
unsigned long long int &dist_counter, vector<TD> &center_mean){

    TD dist = 0;
    TD radius = 0;
    algorithm_utils alg_utils;
    vector<TI> temp1;

    int curr_center = 0, ot_center = 0, cnt = 0;
    float shortestDist2 = 0;
    
    for (int i=0; i<closest_center.size(); i++)
        closest_center[i] = std::numeric_limits<float>::max();
    
    // center mean for Pearson
    calc_center_mean(centroids, center_mean);


    // Calculate inter-centroid distances
    for(curr_center=0; curr_center<centroids.size(); curr_center++){

        center_dist_mat[curr_center][curr_center] = std::numeric_limits<float>::max();
        
        for (ot_center=curr_center+1; ot_center<centroids.size(); 
        ot_center++){    
                
            dist = alg_utils.calc_pearson(centroids[curr_center], centroids[ot_center], center_mean, dist_counter, curr_center,
            ot_center);

            dist /= 2;
            
            center_dist_mat[curr_center][ot_center] = dist;
            center_dist_mat[ot_center][curr_center] = center_dist_mat[curr_center][ot_center];

            if (dist < closest_center[curr_center])
                closest_center[curr_center] = dist;
            
            if (dist < closest_center[ot_center])
                closest_center[ot_center] = dist;

        }

    }

    // Calculate inter-centroid distances
    for(curr_center=0; curr_center<centroids.size(); curr_center++){

        radius = cluster_size[curr_center][1] + closest_center[curr_center];
        cnt = 0;
        
        for (ot_center=0; ot_center<centroids.size(); 
        ot_center++){  

            // Start neighbor finding
            if ((curr_center != ot_center) && 
            (center_dist_mat[curr_center][ot_center] <= radius)){
                                        
                temp1.push_back(ot_center);
           
                // Get the mid-point coordinates for this pair of centroids
                find_midpoints(centroids[curr_center], centroids[ot_center], mid_points, affine_vectors, curr_center, ot_center);
                cnt++;
            }
        
        }   

            if (cnt>=1){
                neighbors[curr_center] = temp1;
            }
 
            else if(cnt == 0){
                temp1.push_back(-100);
                neighbors[curr_center] = temp1;
            }

            cluster_size[curr_center][2] = cnt;
            temp1.clear();
    }
}




template <typename TD, typename TI>
inline void find_neighbors_bounded(vector<vector <TD> > &centroids, 
vector<vector <TD> > &center_dist_mat, vector<TD> &closest_center,
vector<vector <TD> > &cluster_size, 
vector<vector<pair<TD, TI>> > &neighbors, vector<vector<vector <TD> > > &mid_points, 
vector<vector<vector<TD>> > &affine_vectors,
vector<vector<TD> > &temp_master, vector<TD> &temp_midpoint, vector<TD> &temp_affine, 
vector<vector<TD> > &midpoint_holder, vector<vector<TD> > &affine_holder, 
unsigned long long int &dist_counter){

    TD dist = 0;
    TD radius = 0;
    algorithm_utils alg_utils;
    vector<pair<float, int>> temp_pair;

    int curr_center = 0, ot_center = 0, cnt = 0, closest_neighbor = 0, temp_val = 0, closest_ind = 0;
    float shortestDist2 = 0;
    
    for (int i=0; i<closest_center.size(); i++)
        closest_center[i] = std::numeric_limits<float>::max();

    
    // Calculate inter-centroid distances
    for(curr_center=0; curr_center<centroids.size(); curr_center++){

        center_dist_mat[curr_center][curr_center] = std::numeric_limits<float>::max();
        
        for (ot_center=curr_center+1; ot_center<centroids.size(); 
        ot_center++){    
                
                dist = alg_utils.calc_euclidean(centroids[curr_center], centroids[ot_center], dist_counter);
                dist /= 2;
                center_dist_mat[curr_center][ot_center] = dist;
                center_dist_mat[ot_center][curr_center] = center_dist_mat[curr_center][ot_center];

                if (dist < closest_center[curr_center])
                    closest_center[curr_center] = dist;
                
                if (dist < closest_center[ot_center])
                    closest_center[ot_center] = dist;

            }
    }


    // Calculate inter-centroid distances
    for(curr_center=0; curr_center<centroids.size(); curr_center++){

        radius = cluster_size[curr_center][1] + closest_center[curr_center];
        cnt = 0;
        
        for (ot_center=0; ot_center<centroids.size(); 
        ot_center++){  

            // Start neighbor finding
            if ((curr_center != ot_center) && 
            (center_dist_mat[curr_center][ot_center] <= radius)){
                                        
                temp_pair.push_back(std::make_pair(center_dist_mat[curr_center][ot_center], ot_center));
           
                // Get the mid-point coordinates for this pair of centroids
                find_midpoints(centroids[curr_center], centroids[ot_center], mid_points, affine_vectors, curr_center, ot_center);
                cnt++;
            }
        
        }   

            if (cnt>=1){
                std::sort(temp_pair.begin(), temp_pair.end());
                neighbors[curr_center] = temp_pair;
            }
 
            else if(cnt == 0){
                temp_pair.push_back(make_pair(-100, -100));
                neighbors[curr_center] = temp_pair;
            }

            cluster_size[curr_center][2] = cnt;
            temp_pair.clear();
    }
}




template <typename TD, typename Tint>
inline void find_exp_neighbors(vector<vector <TD> > &centroids, 
vector<vector <TD> > &center_dist_mat, vector<TD> &closest_center,
vector<vector<pair<TD, Tint> > > &center_center_dist,
unsigned long long int &dist_counter){

    TD dist = 0;
    TD radius = 0;
    algorithm_utils alg_utils;


    int curr_center = 0, ot_center = 0, cnt = 0;
    
    for (int i=0; i<closest_center.size(); i++)
        closest_center[i] = std::numeric_limits<float>::max();

    
    // Calculate inter-centroid distances
    for(curr_center=0; curr_center<centroids.size(); curr_center++){

        center_dist_mat[curr_center][curr_center] = 0;
        center_center_dist[curr_center][curr_center].first = 0;
        center_center_dist[curr_center][curr_center].second = curr_center;
        
        for (ot_center=curr_center+1; ot_center<centroids.size(); 
        ot_center++){    
                
                dist = alg_utils.calc_euclidean(centroids[curr_center], centroids[ot_center], dist_counter);

                // Enter the center-center distance into the structure for sorting
                center_center_dist[curr_center][ot_center].first = dist;
                center_center_dist[curr_center][ot_center].second = ot_center;

                center_center_dist[ot_center][curr_center].first = dist;
                center_center_dist[ot_center][curr_center].second = curr_center;

                // calculate half the distance between centers and update the center dist matrix
                dist /= 2;
                center_dist_mat[curr_center][ot_center] = dist;
                center_dist_mat[ot_center][curr_center] = center_dist_mat[curr_center][ot_center];

                if (dist < closest_center[curr_center])
                    closest_center[curr_center] = dist;
                
                if (dist < closest_center[ot_center])
                    closest_center[ot_center] = dist;
        }

        std::sort(center_center_dist[curr_center].begin(), center_center_dist[curr_center].end());
    }

    
}


template <typename Tfloat, typename Tint>
inline void determine_data_expression(vector<vector<Tfloat> > &dataset, 
vector<vector <Tfloat> > &centroids, vector<vector<Tfloat> > &cluster_size, 
vector<vector <Tfloat> > &center_dist_mat, vector<Tfloat> &closest_center,
vector<Tint> &assigned_clusters, 
vector<vector<Tint> > &neighbors,
vector<vector<vector <Tfloat> > > &affine_vectors, 
vector<vector<vector <Tfloat> > > &mid_points, 
unsigned long long int &dist_counter, vector<Tint> &temp){


algorithm_utils alg_utils;
Tfloat temp_dist = 0, my_dist = 0, ot_dist = 0, temp_result = 0 ;
Tint i = 0, j = 0, my_cluster = 0;


for (i = 0; i < assigned_clusters.size(); i++){

    my_dist = 0;
    ot_dist = 0;

    my_cluster = assigned_clusters[i];

    // if (cluster_size[my_cluster][2] == 0)
    //     continue;
        
    my_dist = alg_utils.calc_euclidean(dataset[i], centroids[my_cluster], dist_counter);
    temp_dist = my_dist;

    if(my_dist > cluster_size[my_cluster][1])
        cluster_size[my_cluster][1] = my_dist;

    if(my_dist <= closest_center[my_cluster]){
            continue;    
    }
    
         
    for (j=0; j<neighbors[my_cluster].size(); j++){

        if (my_dist > center_dist_mat[my_cluster][neighbors[my_cluster][j]]){

            if (find_context_direction(dataset[i], affine_vectors[my_cluster][neighbors[my_cluster][j]], 
            mid_points[my_cluster][neighbors[my_cluster][j]], temp_result)){

                ot_dist = alg_utils.calc_euclidean(dataset[i], centroids[neighbors[my_cluster][j]], dist_counter);

                if(ot_dist > cluster_size[neighbors[my_cluster][j]][1])
                    cluster_size[neighbors[my_cluster][j]][1] = ot_dist;

                if(ot_dist < temp_dist){
                    cluster_size[assigned_clusters[i]][0] -= 1;
                    assigned_clusters[i] = neighbors[my_cluster][j];
                    cluster_size[neighbors[my_cluster][j]][0] += 1;
                    temp_dist = ot_dist;
                }

            }
        
        }
                
    } 

}

}



template <typename Tfloat, typename Tint>
inline void determine_data_expression_record(vector<vector<Tfloat> > &dataset, 
vector<vector <Tfloat> > &centroids, vector<vector<Tfloat> > &cluster_size, 
vector<vector <Tfloat> > &center_dist_mat, vector<Tfloat> &closest_center,
vector<Tint> &assigned_clusters, 
vector<vector<Tint> > &neighbors,
vector<vector<vector <Tfloat> > > &affine_vectors, 
vector<vector<vector <Tfloat> > > &mid_points, 
unsigned long long int &dist_counter, 
double &le_exec, double &lhe_exec, double &he_exec,
vector<Tint> &temp){


algorithm_utils alg_utils;
Tfloat temp_dist = 0, my_dist = 0, ot_dist = 0, temp_result = 0 ;
Tint i = 0, j = 0, my_cluster = 0;


for (i = 0; i < assigned_clusters.size(); i++){

    my_dist = 0;
    ot_dist = 0;

    my_cluster = assigned_clusters[i];

    // Line 14 (check neighborhood)
    // if (cluster_size[my_cluster][2] == 0){
    //     neighbor_exec++;
    //     continue;
    // }
        
    my_dist = alg_utils.calc_euclidean(dataset[i], centroids[my_cluster], dist_counter);
    temp_dist = my_dist;

    if(my_dist > cluster_size[my_cluster][1])
        cluster_size[my_cluster][1] = my_dist;

    
    // Line 17, LE condition
    if(my_dist <= closest_center[my_cluster]){
        le_exec++;
        continue;    
    }
    
         
    for (j=0; j<neighbors[my_cluster].size(); j++){

        // Line 21, LHE condition
        if (my_dist > center_dist_mat[my_cluster][neighbors[my_cluster][j]]){

            lhe_exec++;

            if (find_context_direction(dataset[i], affine_vectors[my_cluster][neighbors[my_cluster][j]], 
            mid_points[my_cluster][neighbors[my_cluster][j]], temp_result)){

                he_exec++;

                ot_dist = alg_utils.calc_euclidean(dataset[i], centroids[neighbors[my_cluster][j]], dist_counter);

                if(ot_dist > cluster_size[neighbors[my_cluster][j]][1])
                    cluster_size[neighbors[my_cluster][j]][1] = ot_dist;

                if(ot_dist < temp_dist){
                    cluster_size[assigned_clusters[i]][0] -= 1;
                    assigned_clusters[i] = neighbors[my_cluster][j];
                    cluster_size[neighbors[my_cluster][j]][0] += 1;
                    temp_dist = ot_dist;
                }

            }
        
        }
                
    } 

}

}



template <typename Tfloat, typename Tint>
inline void update_bounds(vector<Tfloat> &upper_bound, vector<Tfloat> &lower_bound,
vector<Tfloat> &centroid_motion, vector<Tint> &assigned_clusters){

    int cen = 0;
    float secfurthest = 0, furthest = 0;

    // Furthest moving
    for (int i = 0; i < centroid_motion.size(); i++){
        
        if (centroid_motion[i] > furthest){
            secfurthest = furthest;
            furthest = centroid_motion[i];
            cen = i;
        }
        else if (centroid_motion[i] > secfurthest){
            secfurthest = centroid_motion[i];
        }
    }

    // Updating bounds
    for (int i = 0; i < assigned_clusters.size(); i++){
        upper_bound[i] += centroid_motion[assigned_clusters[i]];

        if (cen == assigned_clusters[i])
            lower_bound[i] -= secfurthest;
        else
            lower_bound[i] -= furthest;
    }

}
            
