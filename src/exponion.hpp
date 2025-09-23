/*
Exponion Kmeans code is adapted from the paper:
See: https://arxiv.org/abs/1602.02514. 
// The code is adpated to return statistics for comparison with Geo-Kmeans.
*/

#include <iostream>
#include <vector>
#include <cmath>
#include "annulus_utils.hpp"
#include <algorithm>

using namespace std;


class Exponion{
    
    template <typename Tfloat, typename Tint>
    output_data exponion(vector<vector <Tfloat> > &dataset,
    Tint num_clusters, Tfloat threshold, Tint num_iterations, 
    Tint numCols, string init_type, Tint seed);

};


// Parameters description:
// dataset: is the data (2D matrix format) where rows are data points and columns are features
// num_clusters (integer): number of clusters 
// num_iterations (integer): maximum number of iteration for which the function will execute
// threshold (float): the threshold value to check convergence.

// The execution stop, if the change in centroids is less than threshold value, or
// the number of iterations are reached, whichever happens first.

// numCols (integer): number of features (dimensionality of data)
// init_type (string): type of initialization, "random", "kmeans++" or read centroids from file
// seed (integer): seed for centroid initialization, if using randmom initialization.
// To learn about parameter usage, see example usage in src/main.cpp



template <typename Tfloat, typename Tint>
inline output_data exponion(vector<vector <Tfloat> > &dataset,
    Tint num_clusters, Tfloat threshold, Tint num_iterations, 
    Tint numCols, string init_type, Tint seed) {
    

    Tint loop_counter = 0;
    vector<vector<Tfloat> > centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<vector<Tfloat> > new_centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<Tint> assigned_clusters(dataset.size(), 0);
    
    vector<Tfloat> upper_bounds(dataset.size(), std::numeric_limits<float>::max());
    vector<Tfloat> lower_bounds(dataset.size(), 0);
    vector<Tfloat> centroid_motion(num_clusters, 0);
    vector<int> cluster_count(num_clusters, 0);  

    vector<vector <Tfloat> > center_dist_mat (num_clusters, vector<Tfloat>(num_clusters, 0));
    vector<Tfloat> closest_center_dist(num_clusters, std::numeric_limits<float>::max());
    vector<vector<pair<Tfloat, Tint>> > center_center_dist(num_clusters);

    unsigned short my_cluster = 0;
    float upper_comparison_bound = 0;
    unsigned long long int dist_counter = 0;

    output_data result;

    Tint i = 0, j = 0, k = 0, l = 0, m = 0;
    Tfloat temp_diff = 0, diff = 0, vec_sum = 0, max = 0, dist2 =0;

    algorithm_utils alg_utils;
    annulus_utils ann;

    vector<pair<Tfloat, Tint> > temp(num_clusters);

    for (i=0; i<num_clusters; i++){
        temp[i].first = 0;
        temp[i].second = 0;
    }
    for (i=0; i<num_clusters; i++){
        center_center_dist[i] = temp;
    }


    alg_utils.init_centroids(centroids, dataset, num_clusters, seed, init_type);

    int tt = 1;
    int my_point = 424;

    alg_utils.calculate_distances_test(dataset, centroids,
    num_clusters, assigned_clusters, dist_counter);


    while (loop_counter < num_iterations) {

        loop_counter++;
        
        // Debugging prints
        // cout << "Counter: " << loop_counter << endl;

        find_exp_neighbors(centroids, center_dist_mat, closest_center_dist, center_center_dist, dist_counter);

        // Debugging prints
        // print_neighbors(center_center_dist, center_center_dist.size(), "MAtrix");

        // Reset cluster counts
        for(i=0; i< cluster_count.size(); i++)
            cluster_count[i] = 0;


        // Debugging prints
        // if (loop_counter == tt){
        //     cout << "Counter: " << loop_counter << endl;
        //     print_vector(closest_center_dist, closest_center_dist.size(), "Closest centers");
        //     // print_2d_vector(centroids, new_centroids.size(), "Centroids");
        //     // for(i =0; i<cOrder.size(); i++){
        //     //     cout << cOrder[i].first << ":" << cOrder[i].second << "\t";
        //     // }
        // }


        // loop over all records
        for (int i = 0; i < dataset.size(); ++i) {
            
            my_cluster = assigned_clusters[i];
            upper_comparison_bound = std::max(closest_center_dist[my_cluster], lower_bounds[i]);

            if (upper_bounds[i] <= upper_comparison_bound) {
                continue;
            }

           
            float u2 = calc_squared_dist(dataset[i], centroids[my_cluster]);
            upper_bounds[i] = sqrt(u2);
            dist_counter++;
            

            if (upper_bounds[i] <= upper_comparison_bound) {
                continue;
            }

            float R = 2*upper_bounds[i] + 2*closest_center_dist[my_cluster];
            auto end = std::lower_bound(center_center_dist[my_cluster].begin(), center_center_dist[my_cluster].end(), make_pair(R, num_clusters));

            for (auto j = center_center_dist[my_cluster].begin(); j!= end; j++) {

                if (j->second == assigned_clusters[i]) continue;

                dist2 = alg_utils.calc_euclidean(dataset[i], centroids[j->second], dist_counter);

                if (dist2 < upper_bounds[i]) {
                    lower_bounds[i] = upper_bounds[i];
                    upper_bounds[i] = dist2;
                    my_cluster = j->second;

                }
                else if (dist2 < lower_bounds[i]) {
                    lower_bounds[i] = dist2;
                }

            }

            if (assigned_clusters[i] != my_cluster) {
                assigned_clusters[i] = my_cluster;
            }
        }


        // Update cluster count
        for(int i=0; i<assigned_clusters.size(); i++){
            cluster_count[assigned_clusters[i]]++ ; 
        }
        
        // Update the centroids based on latest data assignment
        alg_utils.update_centroids_ham(dataset, new_centroids, centroids, centroid_motion,
        assigned_clusters, cluster_count, numCols, dist_counter);


        // Check Convergence
        if (alg_utils.check_convergence_ub_lb(centroid_motion, threshold, diff, temp_diff, i, j)){
            cout << "Convergence at iteration: " << loop_counter << "\n";
            break;
        }

        // Move the new centroids to older
        centroids = new_centroids;

        // Reset centroids
        alg_utils.reinit(new_centroids);  

        // Update bounds
        update_bounds(upper_bounds, lower_bounds, centroid_motion, assigned_clusters);
    }
    
    result.loop_counter = loop_counter;
    result.num_dists = dist_counter;
    // result.centroids = new_centroids;
    // result.sse = get_sse(dataset, new_centroids, cluster_size, assigned_clusters, num_clusters);
    return result;
    
}






