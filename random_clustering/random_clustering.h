#ifndef random_clustering_h
#define random_clustering_h

#include "../utils/utils.h"
#include "../utils/utils_snap.h"



class RandomClustering {
 public:
	RandomClustering(const PUNGraph g, std::string path_graph_name, const std::string& method, const std::string& w_opt) 
	  : _g(g), _mx_nid(g->GetMxNId()), _method(method + '-' + w_opt) {
		check_graph_validity(g);
		parse_alg_param(method);
		load_node_weight(DATA_PATH + path_graph_name + "-node_meta.txt", w_opt);
	}

	const std::string& method() const {
		return _method;
	}

	double node_weight(int node_id) const {
		assert(node_id >= 0 && node_id < _mx_nid);
		if (_node_weight.empty()) {
			return 1.0;
		} else {
			return _node_weight[node_id];
		}
	}

	void gen_partition(VecFlt& partition, int priority_node = -1) const {
		if (_alg == "r_hop_max") {
			gen_nbr_max(partition, priority_node);
		} else if (_alg == "r_net") {
			gen_r_net(partition, priority_node);
		} else if (_alg == "reLDG") {
			gen_reLDG(partition, priority_node);
		} else {
			throw(std::invalid_argument("clustering method not supported: " + _alg));
		}
	}

 protected:
	const PUNGraph _g;
	const int _mx_nid;

	std::string _method, _alg;
	int _param;

	VecFlt _node_weight;

	enum TieBreakingRule {
		UNSPECIFIED,
		FAVOR_FIRST_SELECTED,
		FAVOR_SMALL_DEG,
		FAVOR_LARGE_DEG,
		RANDOM_ORDERING
	} _tie_breaking_rule;

	VecInt _node_deg;

	void parse_alg_param(const std::string& method) {
		VecStr string_pair = split_string(method, '-');
		assert(string_pair.size() == 2);
		_alg = string_pair[0];
		_param = std::stoi(string_pair[1]);

		if (_alg.substr(0, 5) == "r_net") {
			if (_alg.size() == 5) {
				_tie_breaking_rule = RANDOM_ORDERING;
				return;
			}
			if (_alg.size() != 7) {
				throw(std::invalid_argument("clustering method not supported: " + _alg));
			}
			switch (_alg[6]) {
				case 's':
					_tie_breaking_rule = FAVOR_SMALL_DEG;
					get_graph_out_deg_seq(_g, _node_deg);
					break;
				case 'l':
					_tie_breaking_rule = FAVOR_LARGE_DEG;
					get_graph_out_deg_seq(_g, _node_deg);
					break;
				case 'f':
					_tie_breaking_rule = FAVOR_FIRST_SELECTED;
					break;
				default:
					throw(std::invalid_argument("clustering method not supported: " + _alg));
			}
			_alg = "r_net";
		}
	}

	void obtain_tie_breaking(VecInt& seeds) const {
		if (_tie_breaking_rule == FAVOR_FIRST_SELECTED) return;
		if (_tie_breaking_rule == RANDOM_ORDERING) {
			std::shuffle(seeds.begin(), seeds.end(), std::default_random_engine(rand()));
			return;
		}
		if (_tie_breaking_rule == FAVOR_SMALL_DEG) {
			auto deg_smaller = [this] (int i, int j) {
				if (_node_deg[i] < _node_deg[j]) {
					return true;
				} else if (_node_deg[i] == _node_deg[j] && i < j) {
					return true;
				} else {
					return false;
				}
			};
			std::sort(seeds.begin(), seeds.end(), deg_smaller);
			return;
		}
		if (_tie_breaking_rule == FAVOR_LARGE_DEG) {
			auto deg_smaller = [this] (int i, int j) {
				if (_node_deg[i] > _node_deg[j]) {
					return true;
				} else if (_node_deg[i] == _node_deg[j] && i < j) {
					return true;
				} else {
					return false;
				}
			};
			std::sort(seeds.begin(), seeds.end(), deg_smaller);
			return;
		}
	}

	void load_node_weight(const std::string& input_file, const std::string& w_opt) {
		int col;
		if (w_opt == "uniform") col=0;
		else if (w_opt == "d1") col=1;
		else if (w_opt == "d2") col=2;
		else if (w_opt == "eig1") col=6;
		else if (w_opt == "eig2") col=7;
		else throw(std::invalid_argument("w_opt invalid!"));

		if (col == 0) {
			return;
		}
		load_vec_from_file(_node_weight, input_file, col, 1, _mx_nid);
		if (IS_DEBUG) {
			std::cout << "node weights:" << std::endl;
			print_vector(_node_weight);
		}
	}

	void gen_nbr_max(VecFlt& partition, int priority_node = -1) const {
		partition = VecFlt(_mx_nid);
		VecFlt node_num;
		gen_beta_rnd_vector(node_num, _node_weight, _mx_nid);
		if (priority_node != -1) {
			assert(_g->IsNode(priority_node));
			node_num[priority_node] = DBL_MAX;
		}
		if (IS_DEBUG) {
			std::cout << "generated node numbers:" << std::endl;
			print_vector(node_num);
		}

		assign_nbr_max_1_hop(partition, node_num);
		for (int hop_i = 1; hop_i < _param; hop_i++) {
			for (int i = 0; i < _mx_nid; i ++) {
				node_num[i] = partition[i];
			}
			assign_nbr_max_1_hop(partition, node_num);
		}
	}

	void assign_nbr_max_1_hop(VecFlt& partition, const VecFlt& node_num) const {
		#pragma omp parallel num_threads(N_THREADS)
		{
		#pragma omp for
		for (int i = 0; i < _mx_nid; i ++) {
			partition[i] = node_num[i];
			TUNGraph::TNodeI NI = _g->GetNI(i);
			for (int d = 0; d < NI.GetOutDeg(); d++) {
				partition[i] = std::max(partition[i], node_num[NI.GetOutNId(d)]);
			}
		}
		}
	}

	void gen_r_net(VecFlt& partition, int priority_node = -1) const {
		if (priority_node != -1) {
			assert(_g->IsNode(priority_node));
		}

		// Step 1: generate a seed set from random maximal distance-r independent set.
		VecInt seeds;
		get_max_r_independent_set(seeds, _param, priority_node);
		if (IS_DEBUG) {
			std::cout << "Seeds: ";
			print_vector(seeds);
		}
		obtain_tie_breaking(seeds);
		if (IS_DEBUG) {
			std::cout << "Seeds: ";
			print_vector(seeds);
		}

		// Step 2: broadcast around the seeds.
		partition = VecFlt(_mx_nid, _mx_nid);
		for (int k = 0; k < seeds.size(); k++) {
			broadcast_bfs(partition, seeds[k], (_param-1)/2, k);
		}

		// Step 3: find unassigned nodes.
		VecFlt node_label(_mx_nid, _mx_nid);
		std::unordered_set<int> unfinished;
		for (int i = 0; i < _mx_nid; i ++) {
			if (partition[i] < _mx_nid) {
				node_label[i] = partition[i];
			} else {
				unfinished.insert(i);
			}
		}

		// Step 4: iteratively assign the unassigned nodes.
		for (int k = (_param-1)/2 + 1; k < _param; k++) {
			std::unordered_set<int> justFinished;
			for (int i : unfinished) {
				TUNGraph::TNodeI NI = _g->GetNI(i);
				for (int d = 0; d < NI.GetOutDeg(); d++) {
					partition[i] = std::min(partition[i], node_label[NI.GetOutNId(d)]);
				}
				if (partition[i] < _mx_nid) {
					justFinished.insert(i);
				}
			}
			for (int i : justFinished) {
				node_label[i] = partition[i];
				unfinished.erase(i);
			}
		}

		if (!unfinished.empty()) {
			std::cout << "Seeds: ";
			print_vector(seeds);
			std::cout << "hop = " << _param << std::endl;
			print_vector(partition);
			throw std::logic_error("All not nodes have been labeled!");
		}
	}

	void get_max_r_independent_set(VecInt& seeds, int r, int priority_node = -1) const {
		// Select a random maximal (@a r)-independent set: the seeds selected
		//   are at least @a r-hop away from each other, and such node set is maximal
		//   (when @a r = 2, this is the regular maximal independent set).
		// @post: _seeds_ is a vector consisting of node_id of seeds.

		// Step 1. Generate a random ordering of all nodes.
		VecInt node_order;
		gen_rnd_ordering(node_order, _mx_nid, _node_weight, priority_node);
		
		// Step 2. Traverse all nodes in the generated ordering: if a node hasn't
		//   been selected as a seed or marked as being a (@a r-1)-hop neighbor 
		//   of a seed, we select it as a seed and mark all its (@a r-1)-hop 
		//   neighbors.
		seeds.clear();
		VecFlt node_label(_mx_nid, -1);
		for (int k = 0; k < _mx_nid; k ++) {
			if (node_label[node_order[k]] == -1) {
				seeds.push_back(node_order[k]);
				broadcast_bfs(node_label, node_order[k], r-1, k);
			}
		}
	}

	void broadcast_bfs(VecFlt& labels, const int node_id, const int hop, const double label) const {
		// This function marks all nodes within (@a hop)-neighbors of @a node_id
		//   with label @a label. 

		if (hop < 0) {
			throw std::invalid_argument("@a hop must be non-negative!");
		}
		if (!_g->IsNode(node_id)) {
			throw std::invalid_argument("@a node_id must be a valid node of the graph!");
		}

		std::unordered_set<int> visited = {node_id};
		std::queue<int> to_visit;
		to_visit.push(node_id);
		int to_visit_size = to_visit.size();
		for (int k = 0; k < hop; k++) {
			for (int k2 = 0; k2 < to_visit_size; k2++) {
				int i = to_visit.front();
				to_visit.pop();
				labels[i] = label;
				TUNGraph::TNodeI NI = _g->GetNI(i);
				for (int d = 0; d < NI.GetOutDeg(); d++) {
					int nbr_id = NI.GetOutNId(d);
					if (!visited.count(nbr_id)) {
						to_visit.push(nbr_id);
						visited.insert(nbr_id);
					}
				}
			}
			to_visit_size = to_visit.size();
		}
		while (!to_visit.empty()) {
			labels[to_visit.front()] = label;
			to_visit.pop();
		}
	}

	void gen_reLDG(VecFlt& partition, int priority_node = -1) const {
		partition = VecFlt(_mx_nid);

		// Step 1. Generate a random ordering of all nodes.
		VecInt node_order;
		gen_rnd_ordering(node_order, _mx_nid, _node_weight, priority_node);

		double cluster_target_sz = 1.0 * _mx_nid / _param;
		VecInt partition_int(_mx_nid, -1);
		// Step 2: iterate over the node stream 10 passes.
		for (int stream_i = 0; stream_i < 10; stream_i ++) {
			VecFlt cluster_size_budget(_param, cluster_target_sz);
			for (int k = 0; k < node_order.size(); k++) {
				int node_id = node_order[k];
				// Get the intersection size of node_id's 1-hop neighbors with each cluster.
				VecInt intersec_size(_param);
				TUNGraph::TNodeI NI = _g->GetNI(node_id);
				for (int d = 0; d < NI.GetOutDeg(); d ++) {
					int nbr_id = NI.GetOutNId(d);
					if (partition_int[nbr_id] >= 0) {
						intersec_size[partition_int[nbr_id]] ++;
					}
				}

				// Find the optimal cluster based on LDG.
				int opt_cluster_idx = find_opt_reLDG(intersec_size, cluster_size_budget);
				partition_int[node_id] = opt_cluster_idx;
				cluster_size_budget[opt_cluster_idx] --;
			}
		}
		for (int i = 0; i < _mx_nid; i++) {
			partition[i] = partition_int[i];
		}
	}

	int find_opt_reLDG(const VecInt& intersec_size, const VecFlt& cluster_size_budget) const {
		int opt_cluster_idx = 0;
		double opt_obj_val = cluster_size_budget[0] * intersec_size[0];
		for (int cluster_i = 1; cluster_i < _param; cluster_i++) {
			double obj_val = cluster_size_budget[cluster_i] * intersec_size[cluster_i];
			if (obj_val > opt_obj_val) {
				opt_cluster_idx = cluster_i;
				opt_obj_val = obj_val;
			} else if (obj_val == opt_obj_val && cluster_size_budget[cluster_i] > cluster_size_budget[opt_cluster_idx]) {
				// tie breaking favor to cluster with higher budget.
				opt_cluster_idx = cluster_i;
			}
		}
		assert(opt_cluster_idx != -1);
		return opt_cluster_idx;
	}
};


#endif  // random_clustering_h
