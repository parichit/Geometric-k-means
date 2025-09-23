#include <iostream>
#include <string>
#include <iomanip>
#include <chrono>
#include <vector>
#include "IOutils.hpp"
#include "algo_utils.hpp"
#include "misc_utils.hpp"
#include "data_holder.hpp"
// #include "mlpack/methods/kmeans/kmeans.hpp"
#include "lloyd_kmeans.hpp"
#include "geokmeans.hpp"
#include "hamerly.hpp"
#include "annulus.hpp"
#include "exponion.hpp"
#include "ball_kmeans++_xf.hpp"
// #include "doubling_proportion.hpp"
// #include "doubling_clusters.hpp"
// #include "benchmark_clus.hpp"
// #include "benchmark_dims.hpp"
// #include "benchmark_scal.hpp"
// #include "benchmark_real_data.hpp"
// #include "benchmark_real_kplus.hpp"
// #include "benchmark_highdim_data.hpp"
// #include "benchmark_scalar_projection.hpp"
// #include "match_accuracy.hpp"
// #include "ablation_study.hpp"

// #include "geokmeans_record_algo_savings.hpp"
// #include "benchmark_algo_savings.hpp"
// #include "run_onlygk_kpp.hpp"
#include "test_random_uniform.hpp"

using namespace std;
// using namespace mlpack;

int main(int argc, char* argv[]){

    // string basePath = "/u/parishar/scratch/DATASETS/real_data/";
    // string benchmark_type = "benchmark_real_data";

    // string basePath = "/u/parishar/scratch/DATASETS/real_data/";
    // string benchmark_type = "benchmark_real_data_kplus";
    
    // string basePath = "/u/parishar/scratch/DATASETS/clustering_data/";
    // string benchmark_type = "benchmark_clus";

    // string basePath = "/u/parishar/scratch/DATASETS/dims_data/";
    // string benchmark_type = "benchmark_dims";

    // string basePath = "/u/parishar/scratch/DATASETS/scal_data/";
    // string benchmark_type = "benchmark_scal";

    // string basePath = "/u/parishar/scratch/DATASETS/real_data/";
    // string benchmark_type = "doubling_clusters";

    // string basePath = "/u/parishar/scratch/DATASETS/real_data/";
    // string benchmark_type = "doubling_proportion";
    // string vecFlag = "0";

    string basePath = argv[1];
    string benchmark_type = argv[2];
    string vecFlag = argv[3];

    // if (benchmark_type == "doubling_clusters"){
    //     double_clusters(basePath);
    // }
    // else if (benchmark_type == "doubling_proportion"){
    //     double_prop(basePath);
    // }
    // if (benchmark_type == "benchmark_real_data"){
    //     benchmark_on_real_data(basePath);
    // }
    // else if (benchmark_type == "benchmark_real_data_kplus"){
    //     benchmark_on_real_kplus(basePath);
    // }
    // else if (benchmark_type == "benchmark_highdim"){
    //     benchmark_on_highdim(basePath);
    // }
    // else if (benchmark_type == "benchmark_highdim_kplus"){
    //     benchmark_on_highdim_kplus(basePath);
    // }
    // else if (benchmark_type == "benchmark_scal"){
    //     benchmark_scal(basePath);
    // }
    // else if (benchmark_type == "benchmark_clus"){
    //     benchmark_clus(basePath);
    // }
    // else if (benchmark_type == "benchmark_dims"){
    //     benchmark_dims(basePath);
    // }
    // else if (benchmark_type == "ablation"){
    //     ablation(basePath, vecFlag);
    // }
    // else if (benchmark_type == "scal_proj"){
    //     benchmark_scal_proj(basePath);
    // }
    // else if (benchmark_type == "match_acc"){
    //     match_accuracy(basePath);
    // }
    
    
    // if (benchmark_type == "record_savings"){
    //     record_savings(basePath);
    // }

    if (benchmark_type == "random_uni"){
        benchmark_on_test(basePath);
    }

    // if (benchmark_type == "only_gk"){
    //     benchmark_on_real_kplus_only_gk(basePath);
    // }

    return 0;
}