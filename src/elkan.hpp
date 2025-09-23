/*
The following code is adapted from the Elkan's paper:
Url: https://cdn.aaai.org/ICML/2003/ICML03-022.pdf
*/

#include <iostream>
#include <vector>
#include "elkan_utils.hpp"
#include "geokm_utils.hpp"
#include <chrono>
using namespace std;

class Elkan{
    template <typename Tfloat, typename Tint>
    output_data elkan_kmeans(vector<vector <Tfloat> > &dataset,
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
output_data elkan_kmeans(vector<vector <Tfloat> > &dataset, Tint num_clusters, 
Tfloat threshold, Tint num_iterations, Tint numCols, 
string init_type, Tint seed){

    Tint loop_counter = 0;
    vector<vector<Tfloat> > centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<vector<Tfloat> > new_centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<Tint> assigned_clusters(dataset.size(), 0);
    
    vector<Tfloat> upper_bounds(dataset.size(), std::numeric_limits<float>::max());
    vector<vector<Tfloat> > lower_bounds(dataset.size(), vector<Tfloat>(num_clusters, 0));
    vector<Tfloat> centroid_motion(num_clusters, 0);
    
    vector<vector<Tfloat> > centroid_dists(num_clusters, vector<Tfloat>(num_clusters, 0));
    
    vector<Tfloat> closest_center_dist(num_clusters, std::numeric_limits<float>::max());
    vector<int> cluster_count(num_clusters, 0);
   

    Tint my_cluster = 0, i = 0, j = 0, k = 0, l = 0, m = 0;
    Tfloat temp_diff = 0, diff = 0, vec_sum = 0;

    unsigned long long int dist_counter = 0;

    int tt = 1;

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

        vector<bool> needCalc(dataset.size(), true);

        loop_counter++;

        // Find inter-centroid distances and store the distance to the closest centroid
        // Step-1 of Elkan's algorithm
        // For all centroids, find s = 1/2 min(d(c, c'))
        find_centroid_distances(centroids, centroid_dists, closest_center_dist, dist_counter); 

        // Reset new centroids
        alg_utils.reinit(new_centroids);

        // Reinit cluster counts
        for(i=0; i<cluster_count.size(); i++)
            cluster_count[i] = 0; 


        // Loop over entire data
        for(i=0; i<dataset.size(); i++){

            my_cluster = assigned_clusters[i];

            // Step-2: if ub(x) <= s(c(x)) then this point won't change it's membership
            if (upper_bounds[i] <= closest_center_dist[my_cluster]){
                continue;
            }

            
            else{ 

            // Loop over all centers
            for(j=0 ; j< centroids.size(); j++){
                
                // If assignment(x) == c(x) then we don't need to re-calculate the distance
                if (j == my_cluster){
                    continue;
                }
                
                if (upper_bounds[i] <= lower_bounds[i][j]){
                    continue; 
                }
                
                if (upper_bounds[i] <= centroid_dists[my_cluster][j]){
                    continue; 
                }


                // Step-3, a
                // if the distance hasn't been calculated for x before, then calculate d(x, c(x))
                // otherwise, d(x, c(x)) = ub(x) - since we update the bounds at the end of the first iteration
                float dist;
                if (needCalc[i]){
                    
                    dist = alg_utils.calc_euclidean(dataset[i], centroids[my_cluster], dist_counter);
                    
                    lower_bounds[i][my_cluster] = dist;
                    upper_bounds[i] = dist;

                    // Since the distance has been calculated, recheck if we can ignore computations
                    // with other centroids
                    if (upper_bounds[i] <= lower_bounds[i][j]){
                        continue; 
                    }

                    if (upper_bounds[i] <= centroid_dists[my_cluster][j]){
                        continue;
                    }

                    needCalc[i] = false;
                }

                else{
                    dist = upper_bounds[i];
                }

                // Step 3b: if d(x, c(x)) > l(x, c) or d(x, c(x)) > 0.5 d(c(x), c)...
                if ((dist > lower_bounds[i][j]) || (dist > centroid_dists[my_cluster][j]))
                {
                
                    // Compute the distance of the point with the other centroids and update
                    // membership accordingly
                    const float dist2 = alg_utils.calc_euclidean(dataset[i], centroids[j], dist_counter);
                    
                    // Update lower bounds with the distance
                    lower_bounds[i][j] = dist2;

                    // Update membership
                    if (dist2 < dist)
                        {
                            upper_bounds[i] = dist2;
                            assigned_clusters[i] = j;
                        }
                }

            }

        }
        // else
        
        
        }

        // Update cluster count
        for(int i=0; i<assigned_clusters.size(); i++){
            cluster_count[assigned_clusters[i]]++; 
        }

        // if (loop_counter == tt){
        //     for(i =0; i< cluster_count.size(); i++)
        //         cout <<  cluster_count[i] << "\t" ;
        //     // print_2d_vector(cluster_count, cluster_count.size(), "Size of clusters");
        // }

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
        update_elk_bounds(upper_bounds, lower_bounds, centroid_motion, assigned_clusters);

    }
    
    result.loop_counter = loop_counter;
    result.num_dists = dist_counter;
    // result.centroids = new_centroids;
    // result.sse = get_sse(dataset, new_centroids, cluster_size, assigned_clusters, num_clusters);
    return result;
}