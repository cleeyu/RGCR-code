/* 
 * This file provides the function to compute network statistics.
 */


#ifndef ball_size
#define ball_size

#include <utils/utils_snap.h>

// static const int MAX_SAMPLE_SIZE = 100000;


class BallSize {
public:
	BallSize(const std::string&  path_graph_name)
		: _path_graph_name(path_graph_name) {
		_g = read_undirected_graph(DATA_PATH + _path_graph_name);
		parse_path_file_name(_path_graph_name, _path, _graph_name);
		_max_nid = _g->GetMxNId();
	}

	void get_all_nodes_ball_size() const {
		// @pre: the graph _g is 
		//    (1) connected;
		//    (2) node_id's are consecutive: [0, 1, ..., max_nid-1];
		//    (3) contains no self-loops.
		check_graph_validity(_g);

		// Compute the size of each ball centered at every node.
		MatInt ball_sizes(_max_nid);
		size_t diameter = 0;
// 		double sample_prob = 1;
// 		if (_max_nid > MAX_SAMPLE_SIZE) {
// 			sample_prob = 1.0 * MAX_SAMPLE_SIZE / _max_nid;
// 		}
		#pragma omp parallel num_threads(N_THREADS) 
		{
		#pragma omp for reduction (max:diameter)
		for (int node_id = 0; node_id < _max_nid; node_id++) {
// 			std::bernoulli_distribution bern_rv(sample_prob);
// 			std::default_random_engine random_eng(rand());
// 			if (bern_rv(random_eng)) {
				get_node_ball_size(node_id, ball_sizes[node_id]);
				diameter = std::max(diameter, ball_sizes[node_id].size());
// 			}
		}
		}

		// Print the size of each ball centered at every node.
		const std::string OUTPUT_PATH = "results/ball_size/";
		int status = system(("mkdir -p " + OUTPUT_PATH + _path).c_str());
		std::string output_file = OUTPUT_PATH + _path_graph_name + "-ball_size.txt";
		std::ofstream fout(output_file);

		fout << "node_id";
		for (int r = 1; r <= diameter; r++) {
			fout << "\t|B_" << r << "|";
		}
		fout << std::endl;

		for (int node_id = 0; node_id < _max_nid; node_id++) {
			if (ball_sizes[node_id].empty()) continue;
			fout << node_id;
			int r = 0;
			for (; r < ball_sizes[node_id].size(); r++) {
				fout << '\t' << ball_sizes[node_id][r];
			}
			for (; r < diameter; r++) {
				fout << '\t' << ball_sizes[node_id].back();
			}
			fout << std::endl;
		}
	}

	void get_node_ball_size(int node_id, VecInt& ball_sizes) const {
		// Compute the size of each ball centering around node_id.
		// @post: ball_sizes = [|B_1(node_id)|, |B_2(node_id)|, ...].
		ball_sizes.clear();
		std::queue<int> boundary;
		std::unordered_set<int> visited;
		boundary.push(node_id);
		visited.insert(node_id);
		while (!boundary.empty()) {
			int boundary_size = boundary.size();
			for (int i = 0; i < boundary_size; i++) {
				int nid = boundary.front();
				boundary.pop();
				TUNGraph::TNodeI NI = _g->GetNI(nid);
				for (int j = 0; j < NI.GetOutDeg(); j++) {
					int nbr_id = NI.GetOutNId(j);
					if (!visited.count(nbr_id)) {
						boundary.push(nbr_id);
						visited.insert(nbr_id);
					}
				}
			}
			ball_sizes.push_back(visited.size());
		}
		ball_sizes.pop_back();
	}

private:
	PUNGraph _g;
	std::string _path_graph_name, _path, _graph_name;
	int _max_nid;
};



#endif  // ball_size
