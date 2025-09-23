#include <iostream>
#include <vector>
#include "ball_km_utils.hpp"
#pragma once

using namespace std;


// class BallKMeans{
//     template <typename Tfloat, typename Tint>
//     inline output_data ball_km(vector<vector <Tfloat> > &dataset,
//     Tint num_clusters, Tfloat threshold, Tint num_iterations, 
//     Tint numCols, Tint time_limit, string init_type, Tint seed);
// };



template <typename Tfloat, typename Tint>
inline output_data ball_km(vector<vector <Tfloat> > &dataset,
    Tint num_clusters, Tfloat threshold, Tint num_iterations, 
    Tint numCols, Tint time_limit, string init_type, Tint seed) {
    

    Tint loop_counter = 0;
    vector<vector<Tfloat> > centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<vector<Tfloat> > new_centroids(num_clusters, vector<Tfloat>(numCols, 0));
    vector<Tint> assigned_clusters(dataset.size(), 0);
    
    vector<Tfloat> centroid_motion(num_clusters, 0);  
    vector<vector<Tfloat> > cluster_size(num_clusters, vector<Tfloat>(3, 0)); 
    
    vector<vector<Tfloat> > dist_matrix(dataset.size(), vector<Tfloat>(num_clusters));
    vector<vector <Tfloat> > center_dist_mat (num_clusters, vector<Tfloat>(num_clusters, 0));
    vector<Tfloat> closest_center_dist(num_clusters, std::numeric_limits<float>::max());

    vector<vector<pair<Tfloat, Tint>> > neighbors(num_clusters);

    vector<bool> stable_cluster(num_clusters, false);


    unsigned long long int he_counter = 0;

    
    Tint i = 0, j = 0, my_cluster = 0;
    Tfloat temp_diff = 0, diff = 0, dist2 =0;

    int tt = 1;

    
    algorithm_utils alg_utils;
    output_data result;

    alg_utils.init_centroids(centroids, dataset, num_clusters, seed, init_type);

    alg_utils.calculate_distances(dataset, centroids, dist_matrix, num_clusters, 
    assigned_clusters, cluster_size, he_counter);


    while (loop_counter < num_iterations) {

        loop_counter++;
        // cout << "Counter: " << loop_counter << endl;

        // Reset cluster sizes
        for(i=0; i<cluster_size.size(); i++)
            cluster_size[i][0] = 0;

        // reset new centroids
        alg_utils.reinit(new_centroids);

        // print_2d_vector(centroids, centroids.size(), "Centroids");
        // print_2d_vector(new_centroids, new_centroids.size(), "Centroids");

        // Update the centroids based on latest data assignment
        update_centroids_ballkm(dataset, new_centroids, centroids, centroid_motion, dist_matrix,
        stable_cluster, assigned_clusters, cluster_size, numCols, he_counter);


        // Check Convergence
        if (alg_utils.check_convergence_ub_lb(centroid_motion, threshold, diff, temp_diff, i, j)){
            cout << "Convergence at iteration: " << loop_counter << "\n";
            break;
        }

        print_2d_vector(centroids, centroids.size(), "Centroids");
        print_2d_vector(new_centroids, new_centroids.size(), "New Centroids");

        centroids = new_centroids;

        // Establish neighborhood
        find_ball_neighbors(centroids, center_dist_mat, cluster_size, neighbors, he_counter); 

        // print_neighbors(neighbors, neighbors.size(), "neighbors");


        // loop over all records
        for (int i = 0; i < dataset.size(); ++i) {

            my_cluster = assigned_clusters[i];

            // dist2 = dist_matrix[i][my_cluster];
            dist2 = alg_utils.calc_euclidean(dataset[i], centroids[my_cluster], he_counter);
            

            if (cluster_size[my_cluster][2] !=0 && dist2 < neighbors[my_cluster][0].first) {
                continue;
            }


            // Get the list of neighbors based on active areas partitioning
            // auto end = std::lower_bound(neighbors[my_cluster].begin(), neighbors[my_cluster].end(), 
            // make_pair(dist_matrix[i][my_cluster], num_clusters));
            // int t = end - neighbors[my_cluster].begin();
            // cout << "Loop Counter" << loop_counter << "\t" << my_cluster << "\t" << t << "\t" << dist2 << endl;
            
        if (cluster_size[my_cluster][2] !=0){

            for(j = 0; j<neighbors[my_cluster].size(); j++){

                // dist2 = alg_utils.calc_euclidean(dataset[i], new_centroids[j->second], he_counter);
                // dist_matrix[i][j->second] = dist2;

                dist2 = alg_utils.calc_euclidean(dataset[i], centroids[neighbors[my_cluster][j].second], he_counter);
                dist_matrix[i][neighbors[my_cluster][j].second] = dist2;

                // if the distance is less than the currently assigned cluster
                // Update the new cluster and distance
                if (dist2 < dist_matrix[i][assigned_clusters[i]]){
                    // assigned_clusters[i] = j->second;
                    assigned_clusters[i] = neighbors[my_cluster][j].second;
                }
            }

        }  
    }


        // Update cluster sizes
        for(int i=0; i<assigned_clusters.size(); i++){
            cluster_size[assigned_clusters[i]][0]++;

            // Update the cluster radius
            if (dist_matrix[i][assigned_clusters[i]] > cluster_size[assigned_clusters[i]][1])
                cluster_size[assigned_clusters[i]][1] = dist_matrix[i][assigned_clusters[i]];
        }

        // if (tt == loop_counter){
        //     print_2d_vector(cluster_size, cluster_size.size(), "Cluster Size");
        // }


    }
    
    result.loop_counter = loop_counter;
    result.num_he = he_counter;
    result.centroids = new_centroids;
    // result.sse = get_sse(dataset, new_centroids, cluster_size, assigned_clusters, num_clusters);
    return result;
    
}






