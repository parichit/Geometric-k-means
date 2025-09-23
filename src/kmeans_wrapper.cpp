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

/*
If you are using this code, please also cite the paper and give a star to our github repo. 
It helps us to continue our open-source development for the community.

@article{sharma2025geometrickmeansboundfreeapproach,
      title={Geometric-k-means: A Bound Free Approach to Fast and Eco-Friendly k-means}, 
      author={Parichit Sharma and Marcin Stanislaw and Hasan Kurban and Oguzhan Kulekci and Mehmet Dalkilic},
      year={2025},
      eprint={2508.06353},
      archivePrefix={arXiv},
      primaryClass={cs.LG},
      url={https://arxiv.org/abs/2508.06353}, 
}

GitHub: https://github.com/parichit/Geometric-k-means

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

This is the wrapper script to run all the algorithms supported by this library. 
The code below provides an example of how to call each of the algorithms.

Use this library to access the C++ implementations of various K-means algorithms.

The Python wrappers are available as pip pakage at: https://pypi.org/project/DataCentricKMeans/

Currently supported algorithms include:
Elkan, Hamelry, Annulus, Exponion, Geo-Kmeans, Ball K-means++ and Lloyd's K-means.
Note that the Ball K-means++ implementation requires the header files from the Eigen library. We provide the library in the github repository. 

For any questions, please contact: parishar@iu.edu

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

Parameters:
1. filePath: absolute path to the input data file in csv format
2. num_clusters: number of clusters (k)
3. threshold: convergence threshold (e.g., 0.001)
4. num_iterations: maximum number of iterations (e.g., 100)
5. init_type: initialization type for centroids ("random" or absolute path to the file containing the centroids)
6. seed: random seed for centroid initialization (e.g., 42)
7. outPath: absolute path to the output file where results will be written 

Example command to run the code:
./kmeans_wrapper /path/to/data.csv 10 0.001 100 random 12 /path/to/output.csv

or, with custom centroids:
./kmeans_wrapper /path/to/data.csv 10 0.001 100 /path/to/centroids.csv 12 /path/to/output.csv

How to compile the code (assuming all the files are in the src directory):
g++ --std=c++17 -I"/path/to/eigen" kmeans_wrapper.cpp -o kmeans_wrapper

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
*/



int main(int argc, char* argv[]){

    // Read file path
    string filePath = argv[1];

    // clusters (convert string to integer)
    int num_clusters = atoi(argv[2]);
    float threshold = stof(argv[3]);
    int num_iterations = atoi(argv[4]);
    
    // Initialization types
    // 1. "random": initialize the centroids randomly from the data 
    // 2. fileptath: absolute path to the file containing the centroids. 
    // The centroids can be extracted via any method, e.g. k++. This file 
    // should be in csv format with rows as centroids and columns as features.
    
    string init_type = argv[5];
    int seed = atoi(argv[6]);


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

    // Calling the Lloyd's k-means
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


    // Calling the Annulus k-means
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


    // Calling the Exponion k-means
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


    // Calling the Geo-k-means
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

    
    // Calling the Ball-k-means
    // To run Ball-k-means: Make sure to include the Eigen library header files while compiling
    
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

    
    // Calling the Elkan's k-means
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

