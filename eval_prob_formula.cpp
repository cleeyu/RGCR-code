#include "random_clustering/random_clustering.h"
#include "RGCR.h"



int main(int argc, char **argv) {
  reset_random_seed();
  std::cout << "Program start..." << std::endl;
  std::cout << rand() << ' ' << rand() << ' ' << rand() << std::endl << std::endl;

  std::string usage_msg = "Usage: ./" + std::string(argv[0]);
  usage_msg += " -g path_graph_name[synthetic/small_world/SW-2D_32_1_PL-id]";
  usage_msg += " -c clustering_method[nbr_max-1-uniform]";
  usage_msg += " -s file_suffix[3-0.txt]";
  usage_msg += " -r use_total_rand[false]";
  usage_msg += " -h use_hajek[false]";
  usage_msg += " -a base_response[1]";
  usage_msg += " -b drift_std[0.5]";
  usage_msg += " -e noise_std[0.1]";
  usage_msg += " -t GATE[1]";
  usage_msg += " -m multiplicative_ATE[false]";

  // Default run parameters.
  std::string path_graph_name = "synthetic/small_world/SW-2D_32_1_PL-id";
  std::string clustering_method = "nbr_max-1-uniform";
  std::string file_suffix = "3-0.txt";
  bool use_total_rand = false;
  bool use_hajek = false;
  double a = 1;
  double b = 0.5;
  double e = 0.1;
  double GATE = 1.0;
  bool additive_ATE = true;

  extern char *optarg;
  int opt;
  while ((opt = getopt(argc, argv, "g:c:s:rha:b:e:t:m")) != -1) {
    switch (opt) {
      case 'g':
        path_graph_name = optarg;
        break;
      case 'c':
        clustering_method = optarg;
        break;
      case 's':
        file_suffix = optarg;
        break;
      case 'r':
        use_total_rand = true;
        break;
      case 'h':
        use_hajek = true;
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
      default:
        std::cout << usage_msg << std::endl;
        return -1;
    }
  }

  std::string run_name = path_graph_name + ',' + clustering_method + ',' + file_suffix + ',';
  if (use_total_rand) {
    run_name += "tot_randomization,";
  } else {
    run_name += "ind_randomization,";
  }
  if (use_hajek) {
    run_name += "Hajek,";    
  } else {
    run_name += "HT,";
  }
  run_name += std::to_string(a) + ',' + std::to_string(b) + ',' + std::to_string(e);
  if (additive_ATE) {
    run_name += ",additive_TE,";    
  } else {
    run_name += ",multipli_TE,";
  }
  run_name += std::to_string(GATE);
  std::cout << run_name << std::endl;

  std::cout << get_time_str() << ": Experiment starts..."<< std::endl;
  PUNGraph g = read_undirected_graph(GRAPH_CORE_PATH + path_graph_name);
  RGCR rgcr(g, path_graph_name, false);
  rgcr.load_node_response(a, b, e, GATE, additive_ATE);

  std::string output_file_name = "variances-HT.txt";
  if (use_hajek) {
    output_file_name = "variances-Hajek.txt";
  }
  std::ofstream file_out(output_file_name, std::ofstream::app);
  file_out << run_name << std::endl;
  rgcr.eval_expo_prob_formula(clustering_method, use_total_rand, file_suffix, use_hajek, file_out);
  file_out.close();

  std::cout << get_time_str() << ": Experiment finishes..."<< std::endl;
  return 0;
}
