/*
Hamerly Kmeans code is adapted from the github page of the algorithm
See: http://cs.baylor.edu/~hamerly/software/kmeans.php. The code is adpated to
return statistics for comparison with Geo-Kmeans.
*/

#include <iostream>
#include <vector>
#include <map>
#include "geokm_utils.hpp"
#include <chrono>
using namespace std;

class Hamerly{
    template <typename Tfloat, typename Tint>
    output_data hamerly_kmeans(vector<vector <Tfloat> > &dataset,
    Tint num_clusters, Tfloat threshold, Tint num_iterations, 
    Tint numCols, string init_type, Tint seed);
};


template <typename Tfloat, typename Tint>
output_data hamerly_kmeans(vector<vector <Tfloat> > &dataset, Tint num_clusters, 
Tfloat threshold, Tint num_iterations, Tint numCols, string init_type, Tint seed){

    Tint loop_counter = 0;
    vector<vector<Tfloat> > centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<vector<Tfloat> > new_centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<Tint> assigned_clusters(dataset.size(), 0);
    
    vector<Tfloat> upper_bound(dataset.size(), std::numeric_limits<float>::max());
    vector<Tfloat> lower_bound(dataset.size(), 0);
    vector<Tfloat> centroid_motion(num_clusters, 0);
    
    vector<int> cluster_count(num_clusters, 0);  
    vector<Tfloat> closest_center_dist(num_clusters, std::numeric_limits<float>::max());

    Tint my_cluster = 0, i = 0, j = 0, k = 0, l = 0, m = 0;
    Tfloat temp_diff = 0, diff = 0, vec_sum = 0, max123 = 0, dist =0;

    unsigned long long int dist_counter = 0;
    int tt = 12;

    output_data result;

    // Create objects
    algorithm_utils alg_utils;
    dckm_utils dc_utils;

    // Initialize centroids
    alg_utils.init_centroids(centroids, dataset, num_clusters, seed, init_type);

    alg_utils.calculate_distances_test(dataset, centroids,
    num_clusters, assigned_clusters, dist_counter);


    // Iterate until convergence
    while (loop_counter < num_iterations){

        loop_counter++;
        find_closest_center(centroids, closest_center_dist, dist_counter);

        // reset centroids
        alg_utils.reinit(new_centroids);

        for(i=0; i< cluster_count.size(); i++)
            cluster_count[i] = 0; 

        // Loop over entire data
        for(i=0; i<dataset.size(); i++){

            my_cluster = assigned_clusters[i];
            
            max123 = std::max(closest_center_dist[my_cluster], lower_bound[i]);

            if (upper_bound[i] <= max123){
                continue;
            }
            
            upper_bound[i] = alg_utils.calc_euclidean(dataset[i], centroids[my_cluster], dist_counter);

            if (upper_bound[i] <= max123){
                continue;
            }

            // Loop over all centers
            for(j=0 ; j< centroids.size(); j++){
                
                if (j == my_cluster)
                    continue;

                dist = alg_utils.calc_euclidean(dataset[i], centroids[j], dist_counter);

                if (dist < upper_bound[i]){

                    lower_bound[i] = upper_bound[i];
                    upper_bound[i] = dist;
                    assigned_clusters[i] = j;
                }
                else if (dist < lower_bound[i]){
                    lower_bound[i] = dist;
                }
 
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
        update_bounds(upper_bound, lower_bound, centroid_motion, assigned_clusters);

    }
    
    result.loop_counter = loop_counter;
    result.num_dists = dist_counter;
    // result.centroids = new_centroids;
    // result.sse = get_sse(dataset, new_centroids, cluster_size, assigned_clusters, num_clusters);
    return result;
}