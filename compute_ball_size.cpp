#include "ball_size.h"

std::string ERR_USAGE_MSG = "Usage: ./compute_ball_size <path/graph>";





int main(int argc, char **argv) {
  std::string usage_msg = "Usage: ./" + std::string(argv[0]) + " [path_graph_name]";

  std::string PGName = "SW32";
  if (argc >= 2) {
    PGName = argv[1];
  }

  NetworkStatistics ns(PGName);
  ns.get_all_nodes_ball_size();

  return 0;
}
