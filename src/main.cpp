#include <iostream>
#include <string>
#include <tuple>
#include <chrono>
#include <filesystem>
#include "data_holder.hpp"
#include "IOutils.hpp"
#include "algo_utils.hpp"
#include "misc_utils.hpp"
#include "lloyd_kmeans.hpp"
#include "geokmeans.hpp"
#include "elkan.hpp"
#include "hamerly.hpp"
#include "annulus.hpp"
#include "exponion.hpp"
#include "ball_kmeans++_xf.hpp"

using namespace std;

// Function to extract filename from a given path
std::string extractFilename(const std::string& path) {
    return std::filesystem::path(path).filename().string();
}


int main(int argc, char* argv[]){

    // Read file path
    string filePath = argv[1];

    // clusters (convert string to integer)
    int num_clusters = atoi(argv[2]);
    float threshold = stof(argv[3]);
    int num_iterations = atoi(argv[4]);
    
    // initialization types
    // 1. "random": initialize the centroids randomly from the data 
    // 2. fileptath: absolute path to the file containing the centroids. 
    // The centroids can be extracted via any method, e.g. k++. This file 
    // should be in csv format with rows as centroids and columns as features.
    
    string init_type = "random";
    int seed = atoi(argv[5]);


    // A path where you want to write the results
    string outPath = argv[6];
    
    std::vector<vector <float> > dataset;
    vector<int> labels;

    // Read in the data
    cout << "Reading data from: " << filePath << endl;
    
    auto t0 = std::chrono::high_resolution_clock::now();
    std::pair<int, int> file_data = readSimulatedData(filePath, dataset, labels, false, false);
    auto t00 = std::chrono::high_resolution_clock::now();
    auto file_int = std::chrono::duration_cast<std::chrono::milliseconds>(t00 - t0);
    
    cout << "File reading time: " << file_int.count() << " milliseconds\n";
    cout << "Data size: " << dataset.size() << " X " << dataset[0].size() << endl;
    
    int numRows = file_data.first;
    int numCols = file_data.second;
    
    output_data res;     
    algorithm_utils alg_utils;

    ofstream resFile;
    string outFile = outPath;

    resFile.open(outFile, ios::trunc);
    resFile << "Algorithm,Data,Clusters,Iters,Runtime,Distances\n";
    resFile.close();

    // Extract filename from the given path for logging
    std::string data = extractFilename(filePath);

    // Debug - Testing
    cout << "\nAlgo: KMeans," << " Clusters: " << num_clusters << ", Threshold: " << threshold << endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    res = lloyd_kmeans(dataset, num_clusters, threshold, num_iterations, numCols, init_type, seed);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto km_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    // std::cout << "Total KMeans time: " << km_int.count() << " milliseconds\n";

    resFile.open(outFile, ios::app);
    resFile << "KMeans," << data << "," << num_clusters << "," << res.loop_counter << "," 
    << km_int.count() << "," << res.num_dists << "\n";
    resFile.close();


    cout << "\nAlgo: Annulus," << " Clusters: " << num_clusters << ", Threshold: " << threshold << endl;
    auto t3 = std::chrono::high_resolution_clock::now();
    res = annulus(dataset, num_clusters, threshold, num_iterations, numCols, init_type, seed);
    auto t4 = std::chrono::high_resolution_clock::now();
    auto ms_int1 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
    
    // std::cout << "Total Annulus time: " << ms_int1.count() << " milliseconds\n";
    // cout << "Distances: " << res.num_dists << "\t" << res.loop_counter << endl;

    // Write the results to the output file
    resFile.open(outFile, ios::app);
    resFile << "Annulus," << data << "," << num_clusters << "," << res.loop_counter << "," << 
    ms_int1.count() << "," << res.num_dists << "\n";
    resFile.close();


    cout << "\nAlgo: Exponion," << " Clusters: " << num_clusters << ", Threshold: " << threshold << endl;
    auto t5 = std::chrono::high_resolution_clock::now();
    res = exponion(dataset, num_clusters, threshold, num_iterations, numCols, init_type, seed);
    auto t6 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5);
    
    // std::cout << "Total Exponion time: " << ms_int2.count() << " milliseconds\n";
    // cout << "Distances: " << res.num_dists << endl;

    // Write the results to the output file
    resFile.open(outFile, ios::app);
    resFile << "Exponion," << data << "," << num_clusters << "," << res.loop_counter << "," 
    << ms_int2.count() << "," << res.num_dists << "\n";
    resFile.close();


    cout << "\nAlgo: Geo-Kmeans," << " Clusters: " << num_clusters << ", Threshold: " << threshold << endl;
    auto t7 = std::chrono::high_resolution_clock::now();
    res = geokmeans(dataset, num_clusters, threshold, num_iterations, numCols, init_type, seed);
    auto t8 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t8 - t7);
    // std::cout << "Total Geo-Kmeans time: " << ms_int3.count() << " milliseconds\n";
    // cout << "Distances: " << res.num_dists << endl;

    // Write the results to the output file
    resFile.open(outFile, ios::app);
    resFile << "GeoKmeans," << data << "," << num_clusters << "," << res.loop_counter << "," << 
    ms_int3.count() << "," << res.num_dists << "\n";
    resFile.close();

    
    // Load data in Eigen format for Ball KMeans
    MatrixOur BallK_dataset = load_data(filePath);

    cout << "\nAlgo: Ball KM," << " Clusters: " << num_clusters << ", Threshold: " << threshold << endl;
    auto t9 = std::chrono::high_resolution_clock::now();
    res = ball_k_means_Ring(BallK_dataset, true, num_clusters, threshold, num_iterations, init_type, seed);
    auto t10 = std::chrono::high_resolution_clock::now();
    auto ms_int4 = std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9);
    // std::cout << "Total Ball KM time: " << ms_int4.count() << " milliseconds\n";
    // cout << "Distances: " << res.num_dists << endl;

    // Write the results to the output file
    resFile.open(outFile, ios::app);
    resFile << "Ball-KM," << data << "," << num_clusters << "," << res.loop_counter << "," << 
    ms_int4.count() << "," << res.num_dists << "\n";
    resFile.close();

    
    cout << "\nAlgo: Elkan," << " Clusters: " << num_clusters << ", Threshold: " << threshold << endl;
    auto t11 = std::chrono::high_resolution_clock::now();
    res = elkan_kmeans(dataset, num_clusters, threshold, num_iterations, numCols, init_type, seed);

    auto t12 = std::chrono::high_resolution_clock::now();
    auto ms_int5 = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11);
    // std::cout << "Total Elkan time: " << ms_int5.count() << " milliseconds\n";
    
    resFile.open(outFile, ios::app);
    resFile << "Elkan," << data << "," << num_clusters << "," << res.loop_counter << "," 
    << ms_int5.count() << "," << res.num_dists << "\n";
    resFile.close();
    
    return 0;

}

