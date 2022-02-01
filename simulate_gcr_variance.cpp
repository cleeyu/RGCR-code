#include "RGCR.h"


int main(int argc, char **argv) {
  reset_random_seed();
  std::cout << "Program start..." << std::endl;
  std::cout << rand() << ' ' << rand() << ' ' << rand() << std::endl << std::endl;

  std::string usage_msg = "Usage: ./" + std::string(argv[0]);
  usage_msg += " -g path_graph_name[SW32]";
  usage_msg += " -c clustering_method[r_net-3]";
  usage_msg += " -w clustering_weight[uniform]";
  usage_msg += " -h estimator_type[Hajek]";
  usage_msg += " -n n_samples_per_clustering[1000]";
  usage_msg += " -k n_clusterings[100]";
  usage_msg += " -a base_response[1]";
  usage_msg += " -b drift_coef[0.5]";
  usage_msg += " -e noise_coef[0.1]";
  usage_msg += " -t GATE[1]";
  usage_msg += " -m multiplicative_ATE[false]";
  usage_msg += " -o output_file_suffix[""]";

  // Default run parameters.
  std::string path_graph_name = "SW32";
  std::string clustering_method = "r_net-3";
  std::string clustering_node_w_opt = "uniform";
  std::string est_type_str = "Hajek";
  int n_samples_per_clustering = 1000;
  int n_clusterings = 100;
  double a = 1;
  double b = 0.5;
  double e = 0.1;
  double GATE = 1.0;
  bool additive_ATE = true;
  std::string output_file_suffix = "";

  extern char *optarg;
  int opt;
  while ((opt = getopt(argc, argv, "g:c:w:h:n:k:a:b:e:t:mo:")) != -1) {
    switch (opt) {
      case 'g':
        path_graph_name = optarg;
        break;
      case 'c':
        clustering_method = optarg;
        break;
      case 'w':
        clustering_node_w_opt = optarg;
        break;
      case 'h':
        est_type_str = optarg;
        break;
      case 'n':
        n_samples_per_clustering = atoi(optarg);
        break;
      case 'k':
        n_clusterings = atoi(optarg);
        break;
      case 'a':
        a = atof(optarg);
        break;
      case 'b':
        b = atof(optarg);
        break;
      case 'e':
        e = atof(optarg);
        break;
      case 't':
        GATE = atof(optarg);
        break;
      case 'm':
        additive_ATE = false;
        break;
      case 'o':
        output_file_suffix = optarg;
        break;
      default:
        std::cout << usage_msg << std::endl;
        return -1;
    }
  }

  std::string run_name = path_graph_name + ',' + clustering_method + '-' + clustering_node_w_opt + ',';
  run_name += est_type_str + ',';
  run_name += std::to_string(n_samples_per_clustering) + ',';
  run_name += std::to_string(n_clusterings) + ',';  
  run_name += std::to_string(a) + ',' + std::to_string(b) + ',' + std::to_string(e);
  if (additive_ATE) {
    run_name += ",additive_TE,";    
  } else {
    run_name += ",multipli_TE,";
  }
  run_name += std::to_string(GATE);
  std::cout << run_name << std::endl;

  std::cout << get_time_str() << ": Experiment starts..."<< std::endl;
  PUNGraph g = read_undirected_graph(DATA_PATH + path_graph_name);
  RGCR rgcr(g, path_graph_name, false);
  rgcr.load_node_response(a, b, e, GATE, additive_ATE);

  std::string output_file_name = "results/simulation-GCR-" + est_type_str;
  if (output_file_suffix != "") {
    output_file_name += "-" + output_file_suffix;
  }
  output_file_name += ".txt";
  std::ofstream file_out(output_file_name, std::ofstream::app);

  file_out << run_name << std::endl;
  RandomClustering random_clustering(g, path_graph_name, clustering_method, clustering_node_w_opt);
  rgcr.simulate_GCR_variance_distribution(random_clustering, n_clusterings, n_samples_per_clustering, est_type_str, file_out);
  file_out.close();

  std::cout << get_time_str() << ": Experiment finishes..."<< std::endl;
  return 0;
}
