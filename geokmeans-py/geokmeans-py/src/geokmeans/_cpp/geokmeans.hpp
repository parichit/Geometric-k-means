/*
Auhtor: Parichit Sharma
Deptartment of Computer Science, Indiana University, Bloomington.
github: https://github.com/parichit/Geometric-k-means
*/

#include <iostream>
#include <vector>
#include "geokm_utils.hpp"
#include <chrono>
using namespace std;

class GeoKmeans{
    template <typename Tfloat, typename Tint>
    output_data geokmeans(vector<vector <Tfloat> > &dataset,
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
inline output_data geokmeans(vector<vector <Tfloat> > &dataset, Tint num_clusters, 
Tfloat threshold, Tint num_iterations, Tint numCols, 
string init_type, Tint seed){

    Tint loop_counter = 0;
    vector<vector<Tfloat> > centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<vector<Tfloat> > new_centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<Tint> assigned_clusters(dataset.size());
    
    vector<vector<Tfloat> > cluster_size(num_clusters, vector<Tfloat>(3));  
    
    vector<vector<Tfloat> > center_dist_mat (num_clusters, vector<Tfloat>(num_clusters, 0));
    vector<vector<Tint> > neighbors(num_clusters);
    
    vector<vector<vector <Tfloat> > > mid_points(num_clusters, vector<vector<Tfloat> >(num_clusters, vector<Tfloat>(numCols, 0)));
    vector<vector<vector <Tfloat> > > affine_vectors(num_clusters, vector<vector<Tfloat> >(num_clusters, vector<Tfloat>(numCols, 0)));
    vector<Tfloat> closest_center(num_clusters, std::numeric_limits<float>::max());
    vector<Tint> temp1;

    vector<Tfloat> center_mean(num_clusters, 0);

    Tint my_cluster = 0, i = 0, j = 0, k = 0, l = 0, m = 0;
    Tfloat temp_diff = 0, diff = 0, vec_sum = 0;

    unsigned long long int dist_counter = 0;

    output_data result;

    // Create objects
    algorithm_utils alg_utils;


    // Initialize centroids
    alg_utils.init_centroids(centroids, dataset, num_clusters, seed, init_type);


    // Assign data to nearest center
    alg_utils.calculate_distances(dataset, centroids,
    num_clusters, assigned_clusters, cluster_size, dist_counter);


    while (loop_counter < num_iterations){

        loop_counter++;

        // cout << "Counter " << loop_counter << endl;
        alg_utils.update_centroids(dataset, new_centroids, assigned_clusters, cluster_size, numCols);

        // Check Convergence
        if (alg_utils.check_convergence(new_centroids, centroids, threshold, diff, temp_diff, i, j)){
                cout << "Convergence at iteration: " << loop_counter << "\n";
                break;
        }

        
        find_neighbors(new_centroids, center_dist_mat, closest_center, cluster_size, neighbors, 
        mid_points, affine_vectors, dist_counter);

        
        // Reset the radius
        for(int i=0; i<cluster_size.size(); i++)
            cluster_size[i][1] = 0;

        determine_data_expression(dataset, new_centroids, cluster_size, center_dist_mat, closest_center, 
        assigned_clusters, neighbors, affine_vectors, mid_points, 
        dist_counter, temp1);

        // Move the new centroids to older
        centroids = new_centroids;
        
        // reset centroids
        alg_utils.reinit(new_centroids);

    }

    result.loop_counter = loop_counter;
    result.num_dists = dist_counter;
    result.centroids = new_centroids;
    // result.sse = get_sse(dataset, new_centroids, cluster_size, assigned_clusters, num_clusters);

    return result;
}