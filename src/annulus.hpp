/*
Annulus Kmeans code is adapted from the github page of the algorithm
See: http://cs.baylor.edu/~hamerly/software/kmeans.php. The code is adpated to
return statistics for comparison with Geo-Kmeans.

See: http://cs.baylor.edu/~hamerly/software/kmeans.php
*/

#include <iostream>
#include <vector>
#include <cmath>
#include "annulus_utils.hpp"
#include <algorithm>

using namespace std;


class AnnulusKmeans{
    
    template <typename Tfloat, typename Tint>
    output_data annulus(vector<vector <Tfloat> > &dataset,
    Tint num_clusters, Tfloat threshold, Tint num_iterations, 
    Tint numCols, string init_type, Tint seed);

    void sort_means_by_norm(vector<vector<float> > &centroids, 
    vector<std::pair<float, int> > &, int &);

};


inline void sort_means_by_norm(vector<vector<float> > &centroids, 
vector<std::pair<float, int> > &cOrder, int &num_clusters) {

    annulus_utils ann;
    
    // sort the centers by their norms
    for (int c1 = 0; c1 < num_clusters; ++c1) {
        float temp = ann.point_point_inner_product(centroids[c1]);
        cOrder[c1].first = sqrt(temp);
        cOrder[c1].second = c1;
    }
    sort(cOrder.begin(), cOrder.end());
}


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
inline output_data annulus(vector<vector <Tfloat> > &dataset,
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
    vector<vector<Tfloat> > cluster_size(num_clusters, vector<Tfloat>(3, 0)); 
    vector<vector<Tfloat> > dist_matrix(dataset.size(), vector<Tfloat>(num_clusters));

    vector<Tfloat> closest_center_dist(num_clusters, std::numeric_limits<float>::max());

    unsigned short my_cluster = 0;
    float upper_comparison_bound = 0;
    unsigned long long int dist_counter = 0;

    output_data result;

    Tint i = 0, j = 0, k = 0, l = 0, m = 0;
    Tfloat temp_diff = 0, diff = 0, vec_sum = 0, max = 0, dist2 =0;

    algorithm_utils alg_utils;
    annulus_utils ann;

    float *xNorm; 
    unsigned short *guard;

    int v = 0;

    vector<std::pair<float, int> > cOrder(num_clusters); 
    xNorm = new float[dataset.size()];
    guard = new unsigned short[dataset.size()];

    // Initialize structires used in annulus Kmeans
    for (int i = 0; i < num_clusters; ++i) {
        cOrder[i].first = 0.0;
        cOrder[i].second = i;
    }

    std::fill(guard, guard + dataset.size(), 1);
    for (int i = 0; i < dataset.size(); ++i) {
        xNorm[i] = sqrt(ann.point_point_inner_product(dataset[i]));
    }

    alg_utils.init_centroids(centroids, dataset, num_clusters, seed, init_type);


    alg_utils.calculate_distances_test(dataset, centroids,
    num_clusters, assigned_clusters, dist_counter);


    while (loop_counter < num_iterations) {

        loop_counter++;
        
        // Debugging prints
        // cout << "Counter: " << loop_counter << endl;

        // compute the inter-center distances, keeping only the closest distances
        find_closest_center(centroids, closest_center_dist, dist_counter);

        // Sort centers according to their norm
        sort_means_by_norm(centroids, cOrder, num_clusters);

        // Reset cluster counts
        for(i=0; i< cluster_count.size(); i++)
            cluster_count[i] = 0;

        // reset centroids
        alg_utils.reinit(new_centroids);  


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

            float l2 = calc_squared_dist(dataset[i], centroids[guard[i]]);
            lower_bounds[i] = sqrt(l2);
            dist_counter++;

            
            float beta = std::max(lower_bounds[i], upper_bounds[i]);
            auto start = lower_bound(cOrder.begin(), cOrder.end(), make_pair(xNorm[i] - beta, num_clusters));

            for (auto jp = start; jp != cOrder.end(); ++jp) {

                if (jp->second == assigned_clusters[i]) continue;

                dist2 = alg_utils.calc_euclidean(dataset[i], centroids[jp->second], dist_counter);

                if (dist2 < upper_bounds[i]) {
                    lower_bounds[i] = upper_bounds[i];
                    upper_bounds[i] = dist2;
                    
                    guard[i] = my_cluster;
                    my_cluster = jp->second;

                }
                else if (dist2 < lower_bounds[i]) {
                    lower_bounds[i] = dist2;
                    guard[i] = jp->second;
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
        
        // Update bounds
        update_bounds(upper_bounds, lower_bounds, centroid_motion, assigned_clusters);

    }
    
    result.loop_counter = loop_counter;
    result.num_dists = dist_counter;
    // result.centroids = new_centroids;
    // result.sse = get_sse(dataset, new_centroids, cluster_size, assigned_clusters, num_clusters);
    return result;
    
}






