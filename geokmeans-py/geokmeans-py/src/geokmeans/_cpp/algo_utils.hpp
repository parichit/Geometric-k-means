#include <iostream>
#include <vector>
#include "math.h"
#include <map>
#include <algorithm>
#include <cmath>
#include <random>
#pragma once

using namespace std;

class algorithm_utils{
    public:

    template <typename T1, typename T2>
    void calculate_distances(const vector<vector<T1> > &dataset, 
    vector<vector<T1> > &centroids,
    T2 num_clusters, vector<T2> &assigned_clusters, 
    vector<vector<T1> > &cluster_size, unsigned long long int &he_counter);

    template <typename T1, typename T2>
    void update_centroids(vector<vector <T1> > &dataset, 
    vector<vector<T1> > &new_centroids, vector<T2> &assigned_clusters, 
    vector<vector<T1> > &cluster_size, T2 numCols);

    template <typename T1>
    bool check_convergence(vector<vector <T1> > &new_centroids, 
    vector<vector <T1> > &centroids, T1 threshold, float &diff, float &temp_diff, int &i, int &j);

    template <typename T1, typename T2>
    void extract_data(vector<vector <T1> > &dataset, vector<vector <T1> > &extracted_data, 
    T2 num_points, T2 seed);
    
    template <typename T1>
    void reinit(vector<vector<T1> > &);

    template <typename T1, typename T2>
    void init_centroids(vector<vector <T1> > &centroids, 
    vector<vector <T1> > &dataset, T2 num_cluster, T2 seed, string init_type);

    template <typename T1, typename T2>
    string get_centroids_index(vector<vector <T1> > &dataset, 
    T2 num_cluster, T2 seed, string init_type);
    
    void get_ranodm_indices(int *arr, int size, int seed);

    template <typename T1, typename T2>
    float get_sse(vector<vector <T1> > &dataset, vector<vector <T1> > &centroids, vector<vector<T1> > &cluster_size, 
    vector<T2> assigned_clusters, T2 num_cluster);

    template <typename T1>
    float calc_euclidean(const vector<T1> &, const vector<T1> &, unsigned long long int &he_counter);

    template <typename T1>
    float calc_squared_dist(const vector<T1> &point, const vector<T1> &center); 

    template <typename T1, typename T2>
    inline void calculate_distances_ub(const vector<vector<T1> > &dataset, 
    vector<vector<T1> > &centroids, vector<vector<T1> > &dist_mat, vector<T1> &upper_bound,
    T2 num_clusters, vector<T2> &assigned_clusters, vector<vector<T1> > &cluster_size, 
    unsigned long long int &he_counter);

    template <typename T1, typename T2>
    inline void calculate_distances_ub_lb(const vector<vector<T1> > &dataset, 
    vector<vector<T1> > &centroids, vector<vector<T1> > &dist_mat, vector<T1> &upper_bound,
    vector<T1> &lower_bound,
    T2 num_clusters, vector<T2> &assigned_clusters, 
    vector<vector<T1> > &cluster_size, 
    unsigned long long int &he_counter);

    template <typename T1, typename T2>
    inline void update_centroids_ub(vector<vector <T1> > &dataset, 
    vector<vector<T1> > &new_centroids, vector<vector<T1> > &old_centroids, 
    vector<T1> &centroid_motion,
    vector<T2> &assigned_clusters, vector<vector<T1> > &cluster_size, 
    T2 numCols,
    unsigned long long int &he_counter);

    template <typename T1>
    inline bool check_convergence_ub_lb(vector <T1> &centroids_movement, 
    T1 threshold, float &diff, float &temp_diff, int &i, int &j);

    template <typename T1, typename T2>
    inline void update_centroids_ham(vector<vector <T1> > &dataset, 
    vector<vector<T1> > &new_centroids, vector<vector<T1> > &old_centroids, 
    vector<T1> &centroid_motion,
    vector<T2> &assigned_clusters, vector<int> &cluster_count, T2 numCols,
    unsigned long long int &he_counter);

    template <typename T1, typename T2>
    inline void calculate_distances_test(const vector<vector<T1> > &dataset, 
    vector<vector<T1> > &centroids,
    T2 num_clusters, vector<T2> &assigned_clusters,  
    unsigned long long int &he_counter);

    template <typename T1, typename T2>
    inline void calculate_distances_elk(const vector<vector<T1> > &dataset, 
    vector<vector<T1> > &centroids, T2 num_clusters, vector<T2> &assigned_clusters, 
    vector<vector<T1> > &cluster_size, 
    unsigned long long int &he_counter);

    template <typename Tfloat>
    void calc_center_mean(vector<vector<Tfloat> > &centroids, vector<Tfloat> &center_mean);

    template <typename T1>
    T1 calc_pearson(const vector<T1> &, const vector<T1> &, vector<T1> &, unsigned long long int &,
    int &, int &);
    
};


template <typename T1>
void algorithm_utils::reinit(vector<vector<T1> > &container){

    for(int i=0;i<container.size(); i++){
        container[i].assign(container[i].size(), 0);
    }
}

// The following code fragmemt is taken from the following question thread.
// https://stackoverflow.com/questions/20734774/random-array-generation-with-no-duplicates
void shuffle(int *arr, size_t size, int seed)
{
    if (size > 1) 
    {
        size_t i;
        srand(seed);
        for (i = 0; i < size - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (size - i) + 1);
          int t = arr[j];
          arr[j] = arr[i];
          arr[i] = t;
        }
    }
}


// The following function is adapted from the following question thread.
// https://stackoverflow.com/questions/20734774/random-array-generation-with-no-duplicates
void algorithm_utils::get_ranodm_indices(int *arr, int size, int seed)
{
    if (size > 1) 
    {
        int i = 0, j = 0, t = 0;
        srand(seed);
        
        for (i = 0; i < size - 1; i++) 
        {
          j = i + rand() / (RAND_MAX / (size - i) + 1);
          t = arr[j];
          arr[j] = arr[i];
          arr[i] = t;
        }
    }
}


template <typename T1, typename T2>
void algorithm_utils::init_centroids(vector<vector <T1> > &centroids, 
vector<vector <T1> > &dataset, T2 num_cluster, T2 seed, string init_type){

    int i = 0, j = 0, size = 0;
    
    if (dataset.size() >= 1000000){
        size = 10000;
    }else{
        size = dataset.size();
    }
    
    if (init_type == "random"){

        int test_array[size];

        for (i = 0; i<size ; i++){
            test_array[i] = i;
        }

        shuffle(test_array, size, seed);

        for(i=0; i<num_cluster; i++){  
            for(j=0; j <dataset[i].size(); j++){
                centroids[i][j] = dataset[test_array[i]][j];
            }   
        }
    }
    else if (init_type == "sequential"){
        for(i=0; i<num_cluster; i++){  
            for(j=0; j<dataset[0].size(); j++){
                centroids[i][j] = dataset[i][j];
            }   
        }
    }
    else{
        read_kplus_plus_centroids(init_type, centroids, num_cluster);
    }
}




template <typename T1, typename T2>
string algorithm_utils::get_centroids_index(vector<vector <T1> > &dataset, 
T2 num_cluster, T2 seed, string init_type){

    int i = 0, j = 0, size = dataset.size();
    string centroid_str = "";
    
    if (init_type == "random"){

        int test_array[size];

        for (i = 0; i<size ; i++){
            test_array[i] = i;
        }

        shuffle(test_array, size, seed);

        for(i=0; i<num_cluster; i++){

            if (i < num_cluster-1){
                centroid_str = centroid_str + std::to_string(test_array[i]) + "+";
            }
            else if (i == num_cluster-1){
                centroid_str = centroid_str + std::to_string(test_array[i]);
            }
        }
    }
    return centroid_str;
}




template <typename T1, typename T2>
void algorithm_utils::extract_data(vector<vector <T1> > &dataset, vector<vector <T1> > &extracted_data, 
T2 num_points, T2 seed){

    int i = 0, j = 0, size = dataset.size();
    int test_array[size];

    for (i = 0; i<size ; i++){
        test_array[i] = i;
    }

    get_ranodm_indices(test_array, size, seed);
    
    for(i=0; i < num_points; i++){ 
        for(j=0; j<dataset[0].size(); j++){
            extracted_data[i][j] = dataset[test_array[i]][j];
        }   
    }
}


template <typename T1>
inline float algorithm_utils::calc_euclidean(const vector<T1> &point, 
const vector<T1> &center, unsigned long long int &he_counter){
    
    T1 dist = 0.0;
    T1 temp = 0.0;
    
    for (int i=0; i < point.size(); i++){
        temp = point[i] - center[i];
        dist = dist + (temp*temp);
    }
    
    dist = sqrt(dist);
    he_counter++;
    
    return dist;
}


template <typename T1>
float calc_squared_dist(const vector<T1> &point, 
const vector<T1> &center){
    
    T1 dist = 0.0;
    T1 temp = 0.0;
    
    // cout << point.size() << "\n";
    for (int i=0; i < point.size(); i++){
        temp = point[i] - center[i];
        dist = dist + (temp*temp);
    }
    return dist;
}


template <typename T1, typename T2>
float get_sse(vector<vector <T1> > &dataset, vector<vector <T1> > &centroids, vector<vector<T1> > &cluster_size, 
vector<T2> assigned_clusters, T2 num_cluster){

float total_sse = 0;
vector<float> sse_vec(num_cluster, 0);
int i = 0;

for (i = 0; i<dataset.size(); i++){
    sse_vec[assigned_clusters[i]] += calc_squared_dist(dataset[i], centroids[assigned_clusters[i]]);
}

for(i = 0; i< num_cluster;i++){
    sse_vec[i] /= cluster_size[i][0];
    total_sse += sse_vec[i];
}

return total_sse;

}


template <typename T1, typename T2>
inline void algorithm_utils::calculate_distances(const vector<vector<T1> > &dataset, 
vector<vector<T1> > &centroids,
T2 num_clusters, vector<T2> &assigned_clusters, vector<vector<T1> > &cluster_size, 
unsigned long long int &he_counter){

    T2 current_center = 0;
    vector<T1> temp_dist (num_clusters);
    float temp = 0.0;
    float shortestDist2 = 0.0;
    int i =0, j =0;

    assigned_clusters.assign(assigned_clusters.size(), 0);
    algorithm_utils::reinit(cluster_size);

    // Calculate the distance of points to nearest center
    for (i=0; i < dataset.size(); i++){

        shortestDist2 = std::numeric_limits<float>::max();
        
        for (j=0; j < centroids.size(); j++){ 
            
            temp = calc_euclidean(dataset[i], centroids[j], he_counter);
            temp_dist[j] = temp;
            
            if (temp < shortestDist2){
                shortestDist2 = temp;
                current_center = j;
            }
        }
        
        assigned_clusters[i] = current_center;

        // Increase the size of the cluster
        cluster_size[current_center][0] = cluster_size[current_center][0] + 1;
        
        // Store the max so far
        if (shortestDist2 > cluster_size[current_center][1])
            cluster_size[current_center][1] = shortestDist2;
    }
}


template <typename T1, typename T2>
inline void algorithm_utils::calculate_distances_ub(const vector<vector<T1> > &dataset, 
vector<vector<T1> > &centroids, vector<vector<T1> > &dist_mat, vector<T1> &upper_bound,
T2 num_clusters, vector<T2> &assigned_clusters, vector<vector<T1> > &cluster_size, 
unsigned long long int &he_counter){

    T2 current_center = 0;
    vector<T1> temp_dist (num_clusters);
    float temp = 0.0;
    float shortestDist = 0.0;
    int i = 0, j = 0;

    assigned_clusters.assign(assigned_clusters.size(), 0);
    algorithm_utils::reinit(cluster_size);

    // Calculate the distance of points to nearest center
    for (i=0; i < dataset.size(); i++){

        shortestDist = std::numeric_limits<float>::max();

        for (j=0; j < centroids.size(); j++){ 

            temp = calc_euclidean(dataset[i], centroids[j], he_counter);
            temp_dist[j] = temp;
            
            if (temp < shortestDist){
                shortestDist = temp;
                current_center = j;
            }
        }   

        dist_mat[i] = temp_dist;
        assigned_clusters[i] = current_center;
        upper_bound[i] = shortestDist;

        // Increase the size of the cluster
        cluster_size[current_center][0] = cluster_size[current_center][0] + 1;
        
        // Store the max so far
        if (shortestDist > cluster_size[current_center][1])
            cluster_size[current_center][1] = shortestDist;
    }
}


template <typename T1, typename T2>
inline void algorithm_utils::calculate_distances_ub_lb(const vector<vector<T1> > &dataset, 
vector<vector<T1> > &centroids, vector<vector<T1> > &dist_mat, vector<T1> &upper_bound,
vector<T1> &lower_bound,
T2 num_clusters, vector<T2> &assigned_clusters, 
vector<vector<T1> > &cluster_size, 
unsigned long long int &he_counter){

    T2 current_center = 0;
    vector<T1> temp_dist (num_clusters, 0);
    float temp = 0.0;
    float shortestDist = 0.0;
    int i = 0, j = 0;

    assigned_clusters.assign(assigned_clusters.size(), 0);
    algorithm_utils::reinit(cluster_size);

    // Calculate the distance of points to nearest center
    for (i=0; i < dataset.size(); i++){

        shortestDist = std::numeric_limits<float>::max();
        current_center = 0;

        for (j=0; j < centroids.size(); j++){ 

            temp = calc_euclidean(dataset[i], centroids[j], he_counter);
            temp_dist[j] = temp;
            
            if (temp < shortestDist){
                shortestDist = temp;
                current_center = j;
            }
        } 

        dist_mat[i] = temp_dist;
        assigned_clusters[i] = current_center;
        
        // Store the max so far
        if (shortestDist > cluster_size[current_center][1])
            cluster_size[current_center][1] = shortestDist;
    }
}



template <typename T1, typename T2>
inline void algorithm_utils::calculate_distances_test(const vector<vector<T1> > &dataset, 
vector<vector<T1> > &centroids,
T2 num_clusters, vector<T2> &assigned_clusters,  
unsigned long long int &he_counter){

    T2 current_center = 0;
    vector<T1> temp_dist (num_clusters, 0);
    float temp = 0.0;
    float shortestDist = 0.0;
    int i = 0, j = 0;

    assigned_clusters.assign(assigned_clusters.size(), 0);

    // Calculate the distance of points to nearest center
    for (i=0; i < dataset.size(); i++){

        shortestDist = std::numeric_limits<float>::max();
        current_center = 0;

        for (j=0; j < centroids.size(); j++){ 

            temp = calc_euclidean(dataset[i], centroids[j], he_counter);
            temp_dist[j] = temp;
            
            if (temp < shortestDist){
                shortestDist = temp;
                current_center = j;
            }
        } 
        assigned_clusters[i] = current_center;
    }
}



template <typename T1, typename T2>
inline void algorithm_utils::calculate_distances_elk(const vector<vector<T1> > &dataset, 
vector<vector<T1> > &centroids, T2 num_clusters, vector<T2> &assigned_clusters, 
vector<vector<T1> > &cluster_size, 
unsigned long long int &he_counter){

    T2 current_center = 0;
    vector<T1> temp_dist (num_clusters, 0);
    float temp = 0.0;
    float shortestDist = 0.0;
    int i = 0, j = 0;

    assigned_clusters.assign(assigned_clusters.size(), 0);
    algorithm_utils::reinit(cluster_size);

    // Calculate the distance of points to nearest center
    for (i=0; i < dataset.size(); i++){

        shortestDist = std::numeric_limits<float>::max();
        current_center = 0;

        for (j=0; j < centroids.size(); j++){ 

            temp = calc_euclidean(dataset[i], centroids[j], he_counter);
            temp_dist[j] = temp;
            
            if (temp < shortestDist){
                shortestDist = temp;
                current_center = j;
            }
        } 

        assigned_clusters[i] = current_center;

        // Increase the size of the cluster
        cluster_size[current_center][0] = cluster_size[current_center][0] + 1;
        
        // Store the max so far
        if (shortestDist > cluster_size[current_center][1])
            cluster_size[current_center][1] = shortestDist;
    }
}





template <typename T1, typename T2>
inline void algorithm_utils::update_centroids(vector<vector <T1> > &dataset, 
vector<vector<T1> > &new_centroids, vector<T2> &assigned_clusters, 
vector<vector<T1> > &cluster_size, T2 numCols){

    int curr_center = 0, index = 0, k = 0, j =0;

    for (index=0; index<dataset.size(); index++){
        curr_center = assigned_clusters[index];
        
        for (j = 0; j<numCols; j++){
            new_centroids[curr_center][j] = new_centroids[curr_center][j] + dataset[index][j];
        }
    }

    for(int i=0; i<new_centroids.size();i++){
        k = cluster_size[i][0];

        for (j = 0; j < new_centroids[i].size(); j++){
            if (k > 0)
                new_centroids[i][j] = new_centroids[i][j]/k;
        }
    } 

}


template <typename T1, typename T2>
inline void algorithm_utils::update_centroids_ub(vector<vector <T1> > &dataset, 
vector<vector<T1> > &new_centroids, vector<vector<T1> > &old_centroids, 
vector<T1> &centroid_motion,
vector<T2> &assigned_clusters, vector<vector<T1> > &cluster_size,
T2 numCols,
unsigned long long int &he_counter){

    int curr_center = 0, index = 0, k = 0, j =0;

    for (index=0; index<dataset.size(); index++){
        curr_center = assigned_clusters[index];
        
        for (j = 0; j<numCols; j++){
            new_centroids[curr_center][j] += dataset[index][j];
        }
    }

    for(int i=0; i<new_centroids.size();i++){
        k = cluster_size[i][0];

        for (j = 0; j < new_centroids[i].size(); j++){
            if (k > 0)
                new_centroids[i][j] = new_centroids[i][j]/k;
        }
        centroid_motion[i] = calc_euclidean(new_centroids[i], old_centroids[i], he_counter);
    } 

}


template <typename T1, typename T2>
inline void algorithm_utils::update_centroids_ham(vector<vector <T1> > &dataset, 
vector<vector<T1> > &new_centroids, vector<vector<T1> > &old_centroids, 
vector<T1> &centroid_motion,
vector<T2> &assigned_clusters, vector<int> &cluster_count, T2 numCols,
unsigned long long int &he_counter){

    int curr_center = 0, index = 0, k = 0, j =0;

    for (index=0; index<dataset.size(); index++){
        curr_center = assigned_clusters[index];
        
        for (j = 0; j<numCols; j++){
            new_centroids[curr_center][j] += dataset[index][j];
        }
    }

    for(int i=0; i<new_centroids.size();i++){
        
        k = cluster_count[i];

        for (j = 0; j < new_centroids[i].size(); j++){
            if (k > 0)
                new_centroids[i][j] = new_centroids[i][j]/k;
        }
        centroid_motion[i] = calc_euclidean(new_centroids[i], old_centroids[i], he_counter);
    } 

}



template <typename T1>
inline bool algorithm_utils::check_convergence(vector<vector <T1> > &new_centroids, 
vector<vector <T1> > &centroids, T1 threshold, float &diff, float &temp_diff, int &i, int &j){

    temp_diff = 0;
    diff = 0;
        
        for (i=0; i<new_centroids.size(); i++){
            for (j=0; j< new_centroids[i].size(); j++)
                temp_diff = new_centroids[i][j] - centroids[i][j];
                diff += (temp_diff * temp_diff);
        }
        diff = diff/new_centroids.size();
        diff = sqrt(diff);

    if (diff <= threshold)
        return true;
    return false;
}



template <typename T1>
inline bool algorithm_utils::check_convergence_ub_lb(vector <T1> &centroids_movement, 
T1 threshold, float &diff, float &temp_diff, int &i, int &j){

    temp_diff = 0;
    diff = 0;

    for (i=0; i<centroids_movement.size(); i++){  
        diff += (centroids_movement[i] * centroids_movement[i]);
    }
    
    diff = diff/centroids_movement.size();
    diff = sqrt(diff);

    // cout << "Thres: " << diff << endl;

    if (diff <= threshold)
        return true;
    return false;
}


template <typename Tfloat>
void calc_center_mean(vector<vector<Tfloat> > &centroids, vector<Tfloat> &center_mean){

    float temp = 0;
    for(int i = 0; i<centroids.size(); i++){
        temp = 0;
        for(int j = 0 ; j<centroids[0].size(); j++){
            temp = temp + centroids[i][j];
        }
        center_mean[i] = temp/centroids[0].size();
    }
}


template <typename T1>
T1 algorithm_utils::calc_pearson(const vector<T1> &point1, const vector<T1> &point2, vector<T1> &center_mean, 
unsigned long long int &he_counter, int &curr_center, int &ot_center){

    float temp = 0, temp1 = 0, temp2 = 0, acc1 = 0, acc2 = 0;

    // Find mean
    for(int i = 0; i<point1.size(); i++){
        
        temp1 = (point1[i] - center_mean[curr_center]);
        // temp1 *= temp1;
        
        temp2 = (point2[i] - center_mean[ot_center]);
        // temp2 *= temp2;
        
        temp += (temp1 * temp2);
        
        // acc1 += temp1;
        // acc2 += temp2;
    }

    temp = temp/point1.size();
    // temp = (temp / (sqrt(acc1) * sqrt(acc2)) );
    // temp = sqrt(1 -  (temp / (sqrt(acc1) * sqrt(acc2)) ));
    
    return temp;
}