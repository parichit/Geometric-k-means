#include <vector>
#include <Eigen/Dense>
using namespace std;

#pragma ONCE
struct output_data {
    int loop_counter = 0;
    unsigned long long int num_dists = 0;
    
    double neighbor_exec = 0;
    double le_exec = 0;
    double lhe_exec = 0;

    vector<vector<float> > centroids;
    Eigen::MatrixXf ballkm_centroids;
    bool timeout = false;
    float sse = 0;
};
