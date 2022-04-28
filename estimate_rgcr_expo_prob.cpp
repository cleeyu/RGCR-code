#include "RGCR.h"


int main(int argc, char **argv) {
  std::cout << "Program start..." << std::endl;
  reset_random_seed();
  std::cout << rand() << ' ' << rand() << ' ' << rand() << std::endl << std::endl;

  std::string usage_msg = "Usage: ./" + std::string(argv[0]);
  usage_msg += " -g path_graph_name[SW32]";
  usage_msg += " -c clustering_method[r_net-3]";
  usage_msg += " -w clustering_node_w_opt[uniform]";
  usage_msg += " -n n_samples[100]";
  usage_msg += " -i run_id[0]";
  usage_msg += " -t(if using complete randomization [false])";
  usage_msg += " -s(if using stratified clustering sampling [false])";
  usage_msg += " -z treatment_prob";

  // Default run parameters.
  std::string path_graph_name = "SW32";
  std::string clustering_method = "r_net-3";
  std::string clustering_node_w_opt = "uniform";
  int n_samples = 100;
  int run_id = 0;
  bool use_complete_randomization = false;
  bool use_stratified_sampling = false;
  int treatment_prob = 0.5;

  extern char *optarg;
  int opt;
  while ((opt = getopt(argc, argv, "g:c:w:n:i:tsp:z:")) != -1) {
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
      case 'n':
        n_samples = atoi(optarg);
        break;
      case 'i':
        run_id = atoi(optarg);
        break;
      case 't':
        use_complete_randomization = true;
        break;
      case 's':
        use_stratified_sampling = true;
        break;
      case 'p':
        N_THREADS = atoi(optarg);
        break;
      case 'z':
        treatment_prob = atoi(optarg);
        break;
      default:
        std::cout << usage_msg << std::endl;
        return -1;
    }
  }

  std::cout << path_graph_name << ',' << clustering_method << ',' << clustering_node_w_opt << ',' << n_samples << ',' << run_id;
  if (use_complete_randomization) {
    std::cout << ",complete_randomization";
  } else {
    std::cout << ",independent_randomization";    
  }
  if (use_stratified_sampling) {
    std::cout << ",stratified_sampling";
  } else {
    std::cout << ",independent_sampling";    
  }
  std::cout << std::endl;

  std::cout << get_time_str() << ": Experiment starts..."<< std::endl;

  PUNGraph g = read_undirected_graph(DATA_PATH + path_graph_name);
  RGCR rgcr(g, path_graph_name);
  RandomClustering random_clustering(g, path_graph_name, clustering_method, clustering_node_w_opt);

  if (use_stratified_sampling) {
    rgcr.stratified_mixing_analysis(random_clustering, n_samples, use_complete_randomization, run_id, treatment_prob);
  } else {
    rgcr.independent_mixing_analysis(random_clustering, n_samples, use_complete_randomization, run_id, treatment_prob);
  }

  std::cout << get_time_str() << ": Experiment finishes..."<< std::endl;
  return 0;
}
