#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include "data_holder.hpp"
#include "algo_utils.hpp"
#include "IOutils.hpp"
#include "geokmeans.hpp"
#include <chrono>

using namespace std;

// Parameters are as follows
// dataset: is the data (2D matrix format)
// num_clusters, num_iterations and threshold are self-explanatory
// note that if you want the algorithm to stop only when
// centroids are same then pass 0 as the value for threshold.

// time_limit helps to stop the program execution if it's taking
// too long to execute. Pass a value in ms, for example, to allow the
// program to run only for 1 minute, pass 60000 as the value of time_limit.
// if the program timeouts, then it will return a large value for SSE so it can be
// ignored from the calculations

// Added num_restarts and bint as mentioned as the last two parameters.



best_results Geokmeans_rr(vector<vector <float> > &dataset, int num_clusters, 
float threshold, int num_iterations, int numCols, 
int time_limit, int num_restarts, int bint){

    float best_score = std::numeric_limits<float>::max();

    // We will update the following structure
    // everytime we get a better clustering i.e., SSE etc.
    output_data results;
    best_results best;

    // vector<vector<float> > centroids(num_clusters, vector<float>(dataset[0].size(), 0));
    unsigned long long int best_calcs = 0;

    // Run the program iteratively
    for (int i=0; i < num_restarts; i++){

        // centroids = init_centroids(dataset, num_clusters, bint+i);
        results = geokmeans(dataset, num_clusters, threshold, num_iterations, numCols, time_limit, "random", bint+i);

        // cout << "run :" << i << " SSE: " << results.sse << endl;
        if (results.sse < best_score){
            best_score = results.sse;
            best.best_calcs = results.num_dists;
            best.centroids = results.centroids;
            best.sse = best_score;
        }
    }

    // Now this structure called best contains 4 things and it can be
    // accessed by keys
    // 1. To access best score (SSE), best.best_score
    // 3. best calculations -> best.best_calcs
    // 4. best assignments -> best.centroids
    // You can see the structure on line 15-21 in this file.
    return best;
}


// Ensure that file path is always the first parameter

int main(int argc, char* argv[]){

// Read file path
string filePath = argv[1];
string index = filePath.substr(filePath.find_last_of('_') + 1);

// clusters (convert string to integer)
int num_clusters = atoi(argv[2]);

// Hard coding based on RR file that was shared
float threshold = 0;

int num_iterations = atoi(argv[3]);
int num_restarts = atoi(argv[4]);
int bint = atoi(argv[5]);

// A path where you want to write the files
string outPath = argv[6];

best_results best;
int time_limit = 7200000;

std::vector<vector <float> > dataset;
vector<int> labels;

// Load the data into dataset
std::pair<int, int> p = readSimulatedData(filePath, dataset, labels, false, false);
int numCols = p.second;

// Start timer
auto t1 = std::chrono::high_resolution_clock::now();

best = Geokmeans_rr(dataset, num_clusters, threshold, num_iterations, numCols,
time_limit, num_restarts, bint);

// Calculate elapsed time in minutes
auto t2 = std::chrono::high_resolution_clock::now();
auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
float runtime = elapsed_time.count()/1000.0;

// Save to files (phewwwww - this is not that smooth to do in CPP)
// Write the best assignments (each assignment is seperated by a comma)


ofstream outfile;
outfile.open(outPath+"file2_" + index, ios::trunc);

for (int i = 0; i<best.centroids.size(); i++){
    for(int j = 0; j< numCols; j++){
        if (j != numCols-1){
            outfile << best.centroids[i][j] << ",";
        }
        else{
            outfile << best.centroids[i][j];
        }
    }
    outfile << endl;
}
outfile.close();

// Write the best calculations
// Since calculation is a single number so writing it
// in the text file because there is nothing to separate by commas
outfile.open(outPath+"file3_" + index, ios::trunc);
outfile << best.best_calcs << endl;
outfile.close();

// Writing the time
outfile.open(outPath+"file4_" + index, ios::trunc);
outfile << runtime << endl;
outfile.close();

// I think that's it. I think there could be minor issues
// when you run, but I will fix them in real time.

}
