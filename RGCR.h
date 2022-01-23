/* 
 * This code/class is used to compute the variance of the GATE estimator, using
 * 	GCR with single clustering partition, as well as RGCR with a sample of clusterings.
 * The estimator is an Horvitz-Thompson estimator
 * We assume full neighborhood exposure model.
 */


#ifndef RGCR_h
#define RGCR_h

#include "utils/utils.h"
#include "utils/utils_snap.h"
#include "random_clustering/random_clustering.h"

static const int N_EXAMPLES_PER_PRINT = 100;

enum EstimatorType {
	HT,
	HAJEK,
	HAJEK_LINEAR,
	HAJEK_EXPONENTIAL
};

EstimatorType parse_estimator_type(const std::string& s) {
	if (s == "HT") return HT;
	else if (s == "Hajek") return HAJEK;
	else if (s == "Hajek_linear") return HAJEK_LINEAR;
	else if (s == "Hajek_exponential") return HAJEK_EXPONENTIAL;
	else {
		throw std::invalid_argument("Invalid estimator_type to parse!");		
	}
}

class RGCR {
 public:
	RGCR(PUNGraph g, std::string path_graph_name = "Sample/SampleGraph", bool load_default_response=true) : _g(g), _path_graph_name(path_graph_name) {
		check_graph_validity(_g);
		_mx_nid = _g->GetMxNId();
		if (load_default_response) {
			load_default_node_response(3, 1.0, false, true);
		}
	}

	void load_node_response(double a, double b, double e, double tau, bool is_additive) {
		load_base_response(a, b, e, true);
		load_treatment_response(tau, is_additive);
	}

	void independent_mixing_analysis_specified(RandomClustering& random_clustering, const std::vector<int>& mixing_levels, bool use_complete_rand, int run_id = 0) {
		get_output_directory(use_complete_rand, false);
		_clustering_method = random_clustering.method();

		if (mixing_levels.size() == 0) {
			throw std::invalid_argument("List of number of mixings must be nonempty!"); 
		}
		for (int k = 1; k < mixing_levels.size(); k++) {
			if (mixing_levels[k] % mixing_levels[k-1]) {
				throw std::invalid_argument("For hierarchical mixing, every next mixing must be a multiple of the previous!"); 
			}
		}
		int n_mixing_levels = mixing_levels.size();
		initialize_mixing_analysis(n_mixing_levels);

		for (int sample_i = 1; sample_i <= mixing_levels.back(); sample_i++) {
			VecFlt partition;
			random_clustering.gen_partition(partition);
			if (use_complete_rand) {
				analyze_partition_complete(partition);
			} else {
				analyze_partition_ind(partition);
			}
			int k = 0;
			while (k < n_mixing_levels && sample_i % mixing_levels[k] == 0) {
				double bias_mix = compute_mixing_bias(k);
				double var_mix = compute_mixing_var(k, mixing_levels[k]);
				_bias_mixing[k].push_back(bias_mix);
				_variance_mixing[k].push_back(var_mix);
				if (mixing_levels[k] * _bias_mixing[k].size() >= N_EXAMPLES_PER_PRINT) {
					print_result_to_file(_bias_mixing[k], _variance_mixing[k], "mix_" + std::to_string(mixing_levels[k]) + ".txt");
					_bias_mixing[k].clear();
					_variance_mixing[k].clear();
				}
				if (k < n_mixing_levels - 1) {
					add_vector(_sum_expo_prob[k+1], _sum_expo_prob[k]);
					add_matrix(_sum_co_expo_prob[k+1], _sum_co_expo_prob[k]);
					add_matrix(_sum_adv_expo_prob[k+1], _sum_adv_expo_prob[k]);
					_sum_expo_prob[k] = VecFlt(_mx_nid);
					_sum_co_expo_prob[k] = MatFlt(_mx_nid);
					_sum_adv_expo_prob[k] = MatFlt(_mx_nid);
					for (int i = 0; i < _mx_nid; i++) {
						_sum_co_expo_prob[k][i] = VecFlt(i);
						_sum_adv_expo_prob[k][i] = VecFlt(i);
					}
				} else {
					if (sample_i != mixing_levels[k]) {
						std::cout << "Weird! Max mixing level is exceeded?" << std::endl;
					}
					print_expo_prob(k, mixing_levels[k], "mix_" + std::to_string(mixing_levels[k]) + "-" + std::to_string(run_id) + ".txt");					
				}
				k ++;
			}
			if (sample_i % 100 == 0) {
				std::cout << get_time_str() << ": Partition " << sample_i << std::endl;
			}
		}

		if (_variance_single.size() != 0) {
			std::cout << "Weird! Results for single should be empty here!" << std::endl;
			print_result_to_file(_bias_single, _variance_single, "single");
		}
		for (int k = 0; k < n_mixing_levels; k++) {
			if (_variance_mixing[k].size() != 0) {
				std::cout << "Weird! Results for mixing should have all been printed!" << std::endl;
				print_result_to_file(_bias_mixing[k], _variance_mixing[k], "mix_" + std::to_string(k + 1));
			}
		}
	}

	void independent_mixing_analysis(RandomClustering& random_clustering, int n_samples, bool use_complete_rand, int run_id = 0) {
		std::vector<int> mixing_levels;
		int mix_level = 10;
		while (mix_level <= n_samples) {
			mixing_levels.push_back(mix_level);
			mix_level *= 10;
		}
		if (mixing_levels.back() != n_samples) {
			std::cout << "Input n_samples must be a power of 10, but is " << n_samples << ". Use " << mixing_levels.back() << " instead." << std::endl;
		}
		independent_mixing_analysis_specified(random_clustering, mixing_levels, use_complete_rand, run_id);
	}

	void stratified_mixing_analysis(RandomClustering& random_clustering, int n_rounds, bool use_complete_rand, int run_id = 0) {
		get_output_directory(use_complete_rand, true);
		_clustering_method = random_clustering.method();

		std::vector<int> mixing_levels;
		int mix_level = 1;
		while (mix_level <= n_rounds) {
			mixing_levels.push_back(mix_level);
			mix_level *= 2;
		}
		if (mixing_levels.back() != n_rounds) {
			std::cout << "Input n_rounds must be a power of 2, but is " << n_rounds << ". Use " << mixing_levels.back() << " instead." << std::endl;
			n_rounds = mixing_levels.back();
		}
		int n_mixing_levels = mixing_levels.size();
		initialize_mixing_analysis(n_mixing_levels);

		double sum_partition_wt = 0;
		for (int i = 0; i < _mx_nid; i++) {
			sum_partition_wt += random_clustering.node_weight(i);
		}
		for (int round_i = 1; round_i <= n_rounds; round_i++) {
			for (int i = 0; i < _mx_nid; i++) {
				if (i % 100 == 0) {
					std::cout << get_time_str() << ": round " << round_i << ", node " << i << std::endl;
				}
				VecFlt partition;
				random_clustering.gen_partition(partition, i);
				if (use_complete_rand) {
					analyze_partition_complete(partition, random_clustering.node_weight(i));
				} else {
					analyze_partition_ind(partition, random_clustering.node_weight(i));
				}
			}
			int k = 0;
			while (k < n_mixing_levels && round_i % mixing_levels[k] == 0) {
				double bias_mix = compute_mixing_bias(k);
				double var_mix = compute_mixing_var(k, mixing_levels[k] * sum_partition_wt);
				_bias_mixing[k].push_back(bias_mix);
				_variance_mixing[k].push_back(var_mix);
				print_result_to_file(_bias_mixing[k], _variance_mixing[k], "mix_" + std::to_string(mixing_levels[k]) + "n.txt");
				_bias_mixing[k].clear();
				_variance_mixing[k].clear();

				if (k < n_mixing_levels - 1) {
					add_vector(_sum_expo_prob[k+1], _sum_expo_prob[k]);
					add_matrix(_sum_co_expo_prob[k+1], _sum_co_expo_prob[k]);
					add_matrix(_sum_adv_expo_prob[k+1], _sum_adv_expo_prob[k]);
					_sum_expo_prob[k] = VecFlt(_mx_nid);
					_sum_co_expo_prob[k] = MatFlt(_mx_nid);
					_sum_adv_expo_prob[k] = MatFlt(_mx_nid);
					for (int i = 0; i < _mx_nid; i++) {
						_sum_co_expo_prob[k][i] = VecFlt(i);
						_sum_adv_expo_prob[k][i] = VecFlt(i);
					}
				} else {
					if (round_i != mixing_levels[k]) {
						std::cout << "Weird! Number of rounds is exceeded?" << std::endl;
					}
					print_expo_prob(k, mixing_levels[k] * sum_partition_wt, "mix_" + std::to_string(mixing_levels[k]) + "n-" + std::to_string(run_id) + ".txt");					
				}
				k ++;
			}
		}
		print_result_to_file(_bias_single, _variance_single, "single.txt");
	}

	void simulate_GCR_variance_distribution(RandomClustering& random_clustering, int n_clusterings, int n_samples, const std::string& est_type_str, std::ostream& out = std::cout) {
		out << "GATE = " << _mu1 - _mu0 << std::endl;
		for (int clustering_i = 0; clustering_i < n_clusterings; clustering_i++) {
			VecFlt partition;
			random_clustering.gen_partition(partition);
			simulate_GCR_variance(partition, n_samples, est_type_str, out);
		}
	}

	void simulate_GCR_variance(const VecFlt& partition, int n_samples, const std::string& est_type_str, std::ostream& out = std::cout) {
		// GCR only supports independent randomization, not complete randomization.
		initialize_mixing_analysis(1);
		analyze_partition_ind(partition);
		EstimatorType est_type = parse_estimator_type(est_type_str);
		if (est_type != HT) {
			compute_normalizer(est_type);
		}

		double tau_gt = _mu1 - _mu0;
		double sum_bias = 0;
		double sum_SE = 0;
		double sum_QE = 0;
		#pragma omp parallel num_threads(N_THREADS) 
		{
		#pragma omp for reduction (+:sum_bias,sum_SE,sum_QE)
		for (int sample_i = 0; sample_i < n_samples; sample_i++) {
			double hat_tau_1, hat_tau_2;
			if (est_type == HT) {
				simulate_HT(partition, false, hat_tau_1, hat_tau_2);	// use_complete_rand=false
			} else {
				simulate_Hajek(partition, false, hat_tau_1, hat_tau_2);	// use_complete_rand=false
			}
			sum_bias += hat_tau_1 + hat_tau_2 - 2 * tau_gt;
			double SE_1 = (hat_tau_1 - tau_gt) * (hat_tau_1 - tau_gt);
			double SE_2 = (hat_tau_2 - tau_gt) * (hat_tau_2 - tau_gt);
			sum_SE += SE_1 + SE_2;
			sum_QE += SE_1*SE_1 + SE_2*SE_2;
		}
		}
		sum_bias /= n_samples * 2;
		sum_SE /= n_samples * 2;
		sum_QE /= n_samples * 2;

		out << sum_bias << '\t' << sqrt((sum_SE - sum_bias*sum_bias) / n_samples) << '\t';
		out << sum_SE << '\t' << sqrt((sum_QE - sum_SE*sum_SE) / n_samples);
		if (est_type == HT) {
			out << '\t' << _variance_single[0];
		}
		out << std::endl;
	}

	void eval_expo_prob_formula(const std::string& clustering_method, bool use_complete_rand, const std::string& file_suffix, bool use_hajek=false, std::ostream& out = std::cout) {
		get_output_directory(use_complete_rand, true);
		std::string file_prefix = _output_file_directory + "expo_prob/" + clustering_method;
		load_expo_prob(file_prefix, file_suffix);
		double var_1 = 0, var_2 = 0;
		compute_mixing_var_12(0, 1, var_1, var_2, use_hajek);
		out << var_1 << '\t' << var_2 << '\t' << var_1 + var_2 << std::endl;
	}

	void eval_expo_prob_simulation_stratified(RandomClustering& random_clustering, int n_rounds, bool use_complete_rand, 
		const std::string& file_suffix, const std::string& est_type_str, std::ostream& out = std::cout) {
		get_output_directory(use_complete_rand, true);
		_clustering_method = random_clustering.method();
		std::string file_prefix = _output_file_directory + "expo_prob/" + _clustering_method;
		load_expo_prob(file_prefix, file_suffix);

		EstimatorType est_type = parse_estimator_type(est_type_str);
		if (est_type != HT) {
			compute_normalizer(est_type);
		}

		double tau_gt = _mu1 - _mu0;
		out << "GATE = " << _mu1 - _mu0 << std::endl;
		double sum_partition_wt = 0;
		for (int i = 0; i < _mx_nid; i++) {
			sum_partition_wt += random_clustering.node_weight(i);
		}
		double bias = 0, var_bias = 0;
		double mse = 0, var_mse = 0;
		for (int round_i = 1; round_i <= n_rounds; round_i++) {
			std::cout << get_time_str() << ": round " << round_i << std::endl;
			double bias_r = 0, var_bias_r = 0;
			double mse_r = 0, var_mse_r = 0;
			#pragma omp parallel num_threads(N_THREADS) 
			{
			#pragma omp for reduction (+:bias_r,var_bias_r,mse_r,var_mse_r)
			for (int i = 0; i < _mx_nid; i++) {
				VecFlt partition;
				random_clustering.gen_partition(partition, i);
				double hat_tau_1, hat_tau_2;
				if (est_type == HT) {
					simulate_HT(partition, use_complete_rand, hat_tau_1, hat_tau_2);
				} else {
					simulate_Hajek(partition, use_complete_rand, hat_tau_1, hat_tau_2);
				}
				double w_i = random_clustering.node_weight(i);
				double SE_1 = (hat_tau_1 - tau_gt) * (hat_tau_1 - tau_gt);
				double SE_2 = (hat_tau_2 - tau_gt) * (hat_tau_2 - tau_gt);
				bias_r += w_i * (hat_tau_1 + hat_tau_2 - 2 * tau_gt);
				var_bias_r += w_i * w_i * (SE_1 + SE_2);
				mse_r += w_i * (SE_1 + SE_2);
				var_mse_r += w_i * w_i * (SE_1*SE_1 + SE_2*SE_2);
			}
			}
			bias += bias_r;
			var_bias += var_bias_r;
			mse += mse_r;
			var_mse += var_mse_r;
		}
		bias /= sum_partition_wt * n_rounds * 2;
		var_bias /= sum_partition_wt * sum_partition_wt * n_rounds * 2 / _mx_nid;
		var_bias -= bias * bias;
		var_bias /= 2 * n_rounds * _mx_nid;
		mse /= sum_partition_wt * n_rounds * 2;
		var_mse /= sum_partition_wt * sum_partition_wt * n_rounds * 2 / _mx_nid;
		var_mse -= mse * mse;
		var_mse /= 2 * n_rounds * _mx_nid;

		out << bias << '\t' << sqrt(var_bias) << '\t';
		out << mse << '\t' << sqrt(var_mse) << std::endl;
	}

 private:
	PUNGraph _g;
	std::string _path_graph_name, _output_file_directory, _clustering_method;
	int _mx_nid;
	VecFlt _node_response_0, _node_response_1;
	double _mu0, _mu1;
	VecFlt _Q, _normalizer;
	double _Q_bar;
	EstimatorType _estimator_type;


	std::vector<VecFlt> _sum_expo_prob;	// _sum_expo_prob[k][u] = sum_{i=1} ^ mix[k] P_i(u).
	std::vector<MatFlt> _sum_co_expo_prob;	// _sum_co_expo_prob[k][u][v] = sum_{i=1} ^ mix[k] P_i(u,v). (u > v)
	std::vector<MatFlt> _sum_adv_expo_prob;	// _sum_adv_expo_prob[k][u][v] = sum_{i=1} ^ mix[k] P_i(u1,v0). (u > v)
	VecFlt _variance_single;	// List of variances of each single partition instance.
	MatFlt _variance_mixing;	// List of variances at each mixing level instance.
	VecFlt _bias_single;	// List of biases of each single partition instance.
	MatFlt _bias_mixing;	// List of biases at each mixing level instance.

	void compute_Q() {
		_Q = VecFlt(_mx_nid, 1.0);
		#pragma omp parallel num_threads(N_THREADS)
		{
		#pragma omp for
		for (int i = 0; i < _mx_nid; i++) {
			for (int j = 0; j < i; j++) {
				_Q[i] += _sum_co_expo_prob[0][i][j] / _sum_expo_prob[0][j];
			}
			for (int j = i+1; j < _mx_nid; j++) {
				_Q[i] += _sum_co_expo_prob[0][j][i] / _sum_expo_prob[0][j];
			}
			_Q[i] /= _mx_nid * _sum_expo_prob[0][i];
		}
		}
		_Q_bar = std::accumulate(_Q.begin(), _Q.end(), 0.0) / _mx_nid;
	}

	void compute_normalizer(EstimatorType estimator_type) {
		if (_estimator_type == estimator_type) return;

		_estimator_type = estimator_type;
		if (_estimator_type == HAJEK) {
			_normalizer = VecFlt(_mx_nid, 1.0);
			return;
		}

		if (_Q.size() != _mx_nid) {
			compute_Q();
		}
		_normalizer.clear();
		_normalizer.reserve(_mx_nid);
		if (_estimator_type == HAJEK_LINEAR) {
			for (int i = 0; i < _mx_nid; i++) {
				_normalizer.push_back(1 + _Q_bar - _Q[i]);
			}
		} else if (_estimator_type == HAJEK_EXPONENTIAL) {
			for (int i = 0; i < _mx_nid; i++) {
				_normalizer.push_back(exp(_Q_bar -_Q[i]));
			}
		} else {
			throw std::invalid_argument("Invalid estimator_type value!");
		}
	}


	void load_default_node_response(int response_opt, double tau, bool is_additive, bool multiply_deg=true) {
		std::string response_file_name = DATA_PATH + _path_graph_name + "-response.txt";
		load_vec_from_file(_node_response_0, response_file_name, response_opt, 0, _mx_nid);
		if (multiply_deg) {
			double avg_deg = 2.0 * _g->GetEdges() / _g->GetNodes();
			for (TUNGraph::TNodeI NI = _g->BegNI(); NI != _g->EndNI(); NI++) {
				_node_response_0[NI.GetId()] *= NI.GetOutDeg() / avg_deg;
			}
		}
		_mu0 = std::accumulate(_node_response_0.begin(), _node_response_0.end(), 0.0) / _mx_nid;

		load_treatment_response(tau, is_additive);

		if (IS_DEBUG) {
			std::cout << "response_0:" << std::endl;
			print_vector(_node_response_0);
			std::cout << "mu(0) = " << _mu0 << std::endl;
			std::cout << "response_1:" << std::endl;
			print_vector(_node_response_1);
			std::cout << "mu(1) = " << _mu1 << std::endl;
		}
	}

	void load_base_response(double a, double b, double e, bool multiply_deg=true) {
		std::string response_file_name = DATA_PATH + _path_graph_name + "-response.txt";
		VecFlt drift, noise;
		if (b != 0) {
			load_vec_from_file(drift, response_file_name, 1, 0, _mx_nid);
		}
		if (e != 0) {
			load_vec_from_file(noise, response_file_name, 2, 0, _mx_nid);
		}
		b *= 2;
		e *= 10;

		_node_response_0.clear();
		_node_response_0.reserve(_mx_nid);
		double min_response = a;
		for (int i = 0; i < _mx_nid; i++) {
			double response = a;
			if (b != 0) {
				response += b * (drift[i] - 1);
			}
			if (e != 0) {
				response += e * (noise[i] - 1);
			}
			_node_response_0.push_back(response);
			min_response = std::min(min_response, response);
		}
		if (min_response < 0) {
			std::cout << "Negative node response encountered: ";
			std::cout << ", a = " << a;
			std::cout << ", b = " << b;
			std::cout << ", e = " << e << std::endl;
		}
		if (multiply_deg) {
			double avg_deg = 2.0 * _g->GetEdges() / _g->GetNodes();
			for (TUNGraph::TNodeI NI = _g->BegNI(); NI != _g->EndNI(); NI++) {
				_node_response_0[NI.GetId()] *= NI.GetOutDeg() / avg_deg;
			}
		}
		_mu0 = std::accumulate(_node_response_0.begin(), _node_response_0.end(), 0.0) / _mx_nid;
	}

	void load_treatment_response(double tau, bool is_additive) {
		_node_response_1.clear();
		_node_response_1.reserve(_mx_nid);
		if (is_additive) {
			for (int i = 0; i < _mx_nid; i++) {
				_node_response_1.push_back(_node_response_0[i] + tau);
			}
		} else {
			double multiplier = 1.0 + tau / _mu0;
			for (int i = 0; i < _mx_nid; i++) {
				_node_response_1.push_back(_node_response_0[i] * multiplier);
			}
		}
		_mu1 = std::accumulate(_node_response_1.begin(), _node_response_1.end(), 0.0) / _mx_nid;
	}

	void get_output_directory(bool use_complete_rand, bool use_stratified_sampling) {
		_output_file_directory = "results/";
		if (use_complete_rand) {
			_output_file_directory += "com";
		} else {
			_output_file_directory += "ind";
		}
		if (use_stratified_sampling) {
			_output_file_directory += "_stratified";
		}
		_output_file_directory += "/" + _path_graph_name;
		// if (specified) {
		// 	_output_file_directory += "-specified";
		// }
		_output_file_directory += '/';

		int status;
		status = system(("mkdir -p " + _output_file_directory + "bias").c_str());
		status = system(("mkdir -p " + _output_file_directory + "variance").c_str());
		status = system(("mkdir -p " + _output_file_directory + "mse").c_str());
		status = system(("mkdir -p " + _output_file_directory + "expo_prob").c_str());
		status = system(("mkdir -p " + _output_file_directory + "partition_info").c_str());
	}

	VecFlt eval_partition_meta(const VecFlt& partition) const {
		std::unordered_map<double, int> cluster_size_map;
		long num_edge_cut = 0;
		for (TUNGraph::TNodeI NI = _g->BegNI(); NI != _g->EndNI(); NI++) {
			int node_id = NI.GetId();
			cluster_size_map[partition[node_id]] ++;
			for (int d = 0; d < NI.GetOutDeg(); d++) {
				if (partition[node_id] != partition[NI.GetOutNId(d)]) {
					num_edge_cut ++;
				}
			}
		}

		std::vector<int> cluster_size_vec;
		cluster_size_vec.reserve(cluster_size_map.size());
		for (auto const& p : cluster_size_map) {
			cluster_size_vec.push_back(p.second);
		}
		std::sort(cluster_size_vec.begin(), cluster_size_vec.end(), std::greater<int>());
		while (cluster_size_vec.size() < 5) {
			cluster_size_vec.push_back(0);
		}

		VecFlt summary;	// [n_clusters, ratio_edge_cut, ratio_size_1, ratio_size_2, ..., ratio_size_5]
		summary.push_back(cluster_size_vec.size());
		summary.push_back(0.5 * num_edge_cut / _g->GetEdges());
		for (int k = 0; k < 5; k ++) {
			summary.push_back(1.0 * cluster_size_vec[k] / _mx_nid);
		}
		return summary;
	}

	void initialize_mixing_analysis(int n_mixing_levels) {
		_sum_expo_prob = std::vector<VecFlt>(n_mixing_levels);
		_sum_co_expo_prob = std::vector<MatFlt>(n_mixing_levels);
		_sum_adv_expo_prob = std::vector<MatFlt>(n_mixing_levels);
		for (int k = 0; k < n_mixing_levels; k ++) {
			_sum_expo_prob[k] = VecFlt(_mx_nid);
			_sum_co_expo_prob[k] = MatFlt(_mx_nid);
			_sum_adv_expo_prob[k] = MatFlt(_mx_nid);
			for (int i = 0; i < _mx_nid; i ++) {
				_sum_co_expo_prob[k][i] = VecFlt(i);
				_sum_adv_expo_prob[k][i] = VecFlt(i);
			}
		}
		_variance_single.clear();
		_variance_single.reserve(N_EXAMPLES_PER_PRINT);
		_variance_mixing = std::vector<VecFlt>(n_mixing_levels);
		_bias_single.clear();
		_bias_single.reserve(N_EXAMPLES_PER_PRINT);
		_bias_mixing = std::vector<VecFlt>(n_mixing_levels);
	}

	void analyze_partition_ind(const VecFlt& partition, double partition_wt = 1) {
		// Step 1: populate the following data structures:
		// cluster_to_adjacent_nodes[f] is the set of nodes whose B_1(i) intersects with cluster f.
		std::unordered_map<double, std::unordered_set<int>> cluster_to_adjacent_nodes;
		// node_to_adjencent_clusters[u] is the set of clusters intersecting with B_1(i) of node i.
		std::vector<std::unordered_set<double>> node_to_adjencent_clusters(_mx_nid);
		for (int i = 0; i < _mx_nid; i++) {
			cluster_to_adjacent_nodes[partition[i]].insert(i);
			node_to_adjencent_clusters[i].insert(partition[i]);
			TUNGraph::TNodeI NI = _g->GetNI(i);
			for (int d = 0; d < NI.GetOutDeg(); d++) {
				int j = NI.GetOutNId(d);
				cluster_to_adjacent_nodes[partition[j]].insert(i);
				node_to_adjencent_clusters[i].insert(partition[j]);
			}
		}

		// Step 2: compute (1) the exposure probability of each node; 
		// 	(2) compute the number of common clusters for each pair of nodes if > 0, and 
		// 	store this information in the data structure below.
		// n_common_clusters[u][v] is the number of clusters intersects with both B_1(u) and B_1(v) (u > v)
		std::vector<std::vector<int>> n_common_clusters(_mx_nid);
		VecFlt expo_prob(_mx_nid);
		#pragma omp parallel num_threads(N_THREADS)
		{
		#pragma omp for
		for (int i = 0; i < _mx_nid; i++) {
			expo_prob[i] = std::pow(0.5, node_to_adjencent_clusters[i].size());
			_sum_expo_prob[0][i] += expo_prob[i] * partition_wt;

			n_common_clusters[i] = std::vector<int>(i);
			for (double cluster_id : node_to_adjencent_clusters[i]) {
				for (int j : cluster_to_adjacent_nodes.at(cluster_id)) {
				 	if (j < i) {
						n_common_clusters[i][j] ++;
				 	}
				}
			}
		}
		}
		double bias_p = 0;

		// Step 3: compute the dependency score of each pair of nodes, and update the variance.
		double variance_p = 0;
		#pragma omp parallel num_threads(N_THREADS)
		{
		#pragma omp for reduction (+:variance_p)
		for (int i = 0; i < _mx_nid; i++) {
			double variance_i = 0;
			variance_i += (1 / expo_prob[i] - 1) * _node_response_0[i] * _node_response_0[i];
			variance_i += (1 / expo_prob[i] - 1) * _node_response_1[i] * _node_response_1[i];
			variance_i += 2 * _node_response_0[i] * _node_response_1[i];

			for (int j = 0; j < i; j++) {
				double prod_prob = expo_prob[i] * expo_prob[j];
				if (n_common_clusters[i][j] > 0) {
					double score = std::pow(2, n_common_clusters[i][j]);
					variance_i += 2 * (score - 1) * _node_response_0[i] * _node_response_0[j];
					variance_i += 2 * (score - 1) * _node_response_1[i] * _node_response_1[j];
					variance_i += 2 * _node_response_0[i] * _node_response_1[j];
					variance_i += 2 * _node_response_1[i] * _node_response_0[j];
					_sum_co_expo_prob[0][i][j] += prod_prob * score * partition_wt;
				} else {
					_sum_co_expo_prob[0][i][j] += prod_prob * partition_wt;
					_sum_adv_expo_prob[0][i][j] += prod_prob * partition_wt;
				}
			}
			variance_p += variance_i;
		}
		}
		variance_p /= _mx_nid * _mx_nid;

		// Step 4: push the analysis result.
		summarize_partition_analysis(partition, bias_p, variance_p);
	}

	void summarize_partition_analysis(const VecFlt& partition, double bias, double variance) {
		_bias_single.push_back(bias);
		_variance_single.push_back(variance);
		if (_variance_single.size() >= N_EXAMPLES_PER_PRINT) {
			print_result_to_file(_bias_single, _variance_single, "single.txt");
			_bias_single.clear();
			_variance_single.clear();
		}

		std::bernoulli_distribution bern_rv(0.01);
		std::default_random_engine random_eng(rand());
		if (bern_rv(random_eng)) {
		// if (true) {
			VecFlt partition_meta = eval_partition_meta(partition);

			std::string output_file_name = _output_file_directory + "partition_info/" + _clustering_method + "-" + "-partition_info.txt";
			std::ofstream file_output;
			file_output.open(output_file_name, std::ofstream::app);
			file_output << bias << '\t' << variance;
			for (double v : partition_meta) {
				file_output << '\t' << v;
			}
			file_output << std::endl;
			file_output.close();
		}
	}

	void pair_clusters(const std::unordered_map<double, int>& cluster_sz, 
		std::unordered_map<double, double>& pairs) const {
		pairs.clear();
		pairs.reserve(cluster_sz.size());
		std::vector<double> cluster_id;
		cluster_id.reserve(cluster_sz.size());
		for (auto const& p : cluster_sz) {
			cluster_id.push_back(p.first);
		}
		auto size_greater = [&cluster_sz](double c1, double c2) {return cluster_sz.at(c1) > cluster_sz.at(c2);};
		std::sort(cluster_id.begin(), cluster_id.end(), size_greater);

		for (int k = 0; k < cluster_id.size(); k += 2) {
			if (k + 1 < cluster_id.size()) {
				pairs[cluster_id[k]] = cluster_id[k+1];
				pairs[cluster_id[k+1]] = cluster_id[k];
			} else {
				pairs[cluster_id[k]] = cluster_id[k];
			}
		}
	}

	double compute_expo_prob_complete_rand(const std::unordered_set<double>& adjacent_clusters, 
		const std::unordered_map<double, double>& pairs) const {
		std::unordered_map<double, int> assignments;
		return compute_expo_prob_complete_rand(adjacent_clusters, pairs, assignments);
	}

	double compute_expo_prob_complete_rand(const std::unordered_set<double>& adjacent_clusters, 
		const std::unordered_map<double, double>& pairs, 
		std::unordered_map<double, int>& assignments) const {
		assignments.clear();
		for (const double c : adjacent_clusters) {
			if (assignments.count(c) && assignments.at(c) == -1) {
				return 0.0;
			} else {
				assignments[c] = 1;
				if (pairs.at(c) != c) {
					assignments[pairs.at(c)] = -1;
				}
			}
		}
		return std::pow(0.5, (assignments.size()+1)/2 );
	}

	void analyze_partition_complete(const VecFlt& partition, double partition_wt = 1) {
		// Step 0: populate the following data structures
		// cluster_sz[f] is the size of cluster f (number of nodes).
		std::unordered_map<double, int> cluster_sz;
		// node_to_adjencent_clusters[u] is the set of clusters intersecting with B_1(u).
		std::vector<std::unordered_set<double>> node_to_adjencent_clusters(_mx_nid);
		for (int i = 0; i < _mx_nid; i++) {
			cluster_sz[partition[i]]++;
			node_to_adjencent_clusters[i].insert(partition[i]);
			TUNGraph::TNodeI NI = _g->GetNI(i);
			for (int d = 0; d < NI.GetOutDeg(); d++) {
				node_to_adjencent_clusters[i].insert(partition[NI.GetOutNId(d)]);
			}
		}

		// Step 1: pair the clusters for complete randomization.
		std::unordered_map<double, double> pairs;
		pair_clusters(cluster_sz, pairs);

		// Step 2: compute the exposure probability of each node, and the bias.
		double bias_p = 0;
		VecFlt expo_prob(_mx_nid);
		#pragma omp parallel num_threads(N_THREADS) 
		{
		#pragma omp for reduction (+:bias_p)
		for (int i = 0; i < _mx_nid; i++) {
			expo_prob[i] = compute_expo_prob_complete_rand(node_to_adjencent_clusters[i], pairs);
			_sum_expo_prob[0][i] += expo_prob[i] * partition_wt;
			if (expo_prob[i] == 0) {
				bias_p += _node_response_1[i] - _node_response_0[i];
			}
		}
		}
		bias_p /= _mx_nid;

		// Step 3: compute the co-exposure and adv-exposure probability of each node pair, and the variance.
		double variance_p = 0;
		#pragma omp parallel num_threads(N_THREADS)
		{
		#pragma omp for reduction (+:variance_p)
		for (int i = 0; i < _mx_nid; i++) {
			if (expo_prob[i] == 0) continue;

			double variance_i = 0;
			variance_i += (1 / expo_prob[i] - 1) * _node_response_0[i] * _node_response_0[i];
			variance_i += (1 / expo_prob[i] - 1) * _node_response_1[i] * _node_response_1[i];
			variance_i += 2 * _node_response_0[i] * _node_response_1[i];

			std::unordered_map<double, int> assignments;
			compute_expo_prob_complete_rand(node_to_adjencent_clusters[i], pairs, assignments);
			for (int j = 0; j < i; j++) {
				if (expo_prob[j] == 0) continue;

				int n_common = 0, n_conflict = 0;
				for (const double c : node_to_adjencent_clusters[j]) {
					if (assignments.count(c)) {
						if (assignments.at(c) == 1) {
							n_common++;
						} else if (assignments.at(c) == -1) {
							n_conflict++;
						} else {
							throw std::logic_error("Cluster assignment should be either +1/-1!");
						}
					}
				}
				double prod_prob = expo_prob[i] * expo_prob[j];
				if (n_conflict == 0) {
					// no conflicts for co-exposure
					double dep_coef = std::pow(2, n_common);
					variance_i += 2 * (dep_coef - 1) * _node_response_0[i] * _node_response_0[j];
					variance_i += 2 * (dep_coef - 1) * _node_response_1[i] * _node_response_1[j];
					_sum_co_expo_prob[0][i][j] += prod_prob * dep_coef * partition_wt;
				} else {
					// Not possible for both nodes exposed to treatment at the same time.
					variance_i -= 2 * _node_response_0[i] * _node_response_0[j];
					variance_i -= 2 * _node_response_1[i] * _node_response_1[j];
				}
				if (n_common == 0) {
					// no conflicts for adv-exposure
					double dep_coef = std::pow(2, n_conflict);
					variance_i -= 2 * (dep_coef - 1) * _node_response_0[i] * _node_response_1[j];
					variance_i -= 2 * (dep_coef - 1) * _node_response_0[j] * _node_response_1[i];
					_sum_adv_expo_prob[0][i][j] += prod_prob * dep_coef * partition_wt;
				} else {
					// Not possible for one exposed to treatment and the other to control.
					variance_i += 2 * _node_response_0[i] * _node_response_1[j];
					variance_i += 2 * _node_response_0[j] * _node_response_1[i];
				}
			}
			variance_p += variance_i;
		}
		}
		variance_p /= _mx_nid * _mx_nid;

		// Step 4: push the analysis result.
		summarize_partition_analysis(partition, bias_p, variance_p);
	}

	double compute_mixing_bias(int vec_id) const {
		double bias = 0;
		for (int i = 0; i < _mx_nid; i++) {
			if (_sum_expo_prob[vec_id][i] == 0) {
				bias += _node_response_1[i] - _node_response_0[i];
			}
		}
		return bias / _mx_nid;
	}

	double compute_mixing_var(int vec_id, double n_mix) const {
		double var_1 = 0, var_2 = 0;
		compute_mixing_var_12(vec_id, n_mix, var_1, var_2);
		return var_1 + var_2; 
	}

	void compute_mixing_var_12(int vec_id, double n_mix, double& var_1, double& var_2, bool use_hajek=false) const {
		if (use_hajek) {
			// Evaluate the Taylor linear approximation of the Hajek estimator.
			VecFlt node_response_0_hajek, node_response_1_hajek;
			node_response_0_hajek.reserve(n_mix);
			node_response_1_hajek.reserve(n_mix);
			for (int i = 0; i < _mx_nid; i++) {
				node_response_0_hajek.push_back(_node_response_0[i] - _mu0);
				node_response_1_hajek.push_back(_node_response_1[i] - _mu1);
			}
			compute_mixing_var_12_helper(vec_id, n_mix, var_1, var_2, node_response_0_hajek, node_response_1_hajek);
		} else {
			compute_mixing_var_12_helper(vec_id, n_mix, var_1, var_2, _node_response_0, _node_response_1);
		}
	}

	void compute_mixing_var_12_helper(int vec_id, double n_mix, double& var_1, double& var_2,
		const VecFlt& node_response_0, const VecFlt& node_response_1) const {
		var_1 = 0;
		var_2 = 0;
		#pragma omp parallel num_threads(N_THREADS) 
		{
		#pragma omp for reduction (+:var_1,var_2)
		for (int i = 0; i < _mx_nid; i++) {
			if (_sum_expo_prob[vec_id][i] == 0) continue;
			double variance_1i = 0;
			double expo_prob_i = _sum_expo_prob[vec_id][i] / n_mix;
			variance_1i += (1 / expo_prob_i - 1) * node_response_0[i] * node_response_0[i];
			variance_1i += (1 / expo_prob_i - 1) * node_response_1[i] * node_response_1[i];
			variance_1i += 2 * node_response_1[i] * node_response_0[i];
			var_1 += variance_1i;
			double variance_2i = 0;
			for (int j = 0; j < i; j++) {
				if (_sum_expo_prob[vec_id][j] == 0) continue;
				double prod_prob = expo_prob_i * _sum_expo_prob[vec_id][j];
				double score = _sum_co_expo_prob[vec_id][i][j] / prod_prob;
				variance_2i += 2 * (score - 1) * node_response_0[i] * node_response_0[j];
				variance_2i += 2 * (score - 1) * node_response_1[i] * node_response_1[j];
				score = _sum_adv_expo_prob[vec_id][i][j] / prod_prob;
				variance_2i -= 2 * (score - 1) * node_response_0[i] * node_response_1[j];
				variance_2i -= 2 * (score - 1) * node_response_1[i] * node_response_0[j];
			}
			var_2 += variance_2i;
		}
		}
		var_1 /= _mx_nid * _mx_nid;
		var_2 /= _mx_nid * _mx_nid;
	}

	void print_result_to_file(const VecFlt& biases, const VecFlt& variances, const std::string& file_suffix) const {
		if (biases.size() != variances.size()) {
			throw std::invalid_argument("Bias and variance should have the same length");
		}
		if (variances.size() == 0) return;

		std::ofstream file_output;
		file_output.open(_output_file_directory + "bias/" + _clustering_method + "-bias_" + file_suffix, std::ofstream::app);
		for (double v : biases) {
			file_output << v << std::endl;
		}
		file_output.close();

		file_output.open(_output_file_directory + "variance/" + _clustering_method + "-var_" + file_suffix, std::ofstream::app);
		for (double v : variances) {
			file_output << v << std::endl;
		}
		file_output.close();

		file_output.open(_output_file_directory + "mse/" + _clustering_method + "-mse_" + file_suffix, std::ofstream::app);
		for (int k = 0; k < variances.size(); k++) {
			file_output << variances[k] + biases[k] * biases[k] << std::endl;
		}
		file_output.close();
	}

	void print_expo_prob(int vec_id, double n_mixing, const std::string& file_suffix) const {
		#pragma omp parallel sections
		{
		{
		std::string output_file_name = _output_file_directory + "expo_prob/" + _clustering_method + "-expo_prob_" + file_suffix;
		std::ofstream file_output(output_file_name);
		for (double v : _sum_expo_prob[vec_id]) {
			file_output << v / n_mixing << std::endl;
		}
		file_output.close();
		}
		#pragma omp section
		{
		std::string output_file_name = _output_file_directory + "expo_prob/" + _clustering_method + "-co_expo_prob_" + file_suffix;
		std::ofstream file_output(output_file_name);
		for (const VecFlt& row : _sum_co_expo_prob[vec_id]) {
			for (double v : row) {
				file_output << v / n_mixing << '\t';
			}
			file_output << std::endl;
		}
		file_output.close();
		}
		#pragma omp section
		{
		std::string output_file_name = _output_file_directory + "expo_prob/" + _clustering_method + "-adv_expo_prob_" + file_suffix;
		std::ofstream file_output(output_file_name);
		for (const VecFlt& row : _sum_adv_expo_prob[vec_id]) {
			for (double v : row) {
				file_output << v / n_mixing << '\t';
			}
			file_output << std::endl;
		}
		file_output.close();
		}
		}
	}

	void load_expo_prob(const std::string& file_prefix, const std::string& file_suffix) {
		#pragma omp parallel sections
		{
		{
		_sum_expo_prob = std::vector<VecFlt>(1);
		load_vec_from_file(_sum_expo_prob[0], file_prefix+"-expo_prob_"+file_suffix, 0, 0, _mx_nid);
		}
		#pragma omp section
		{
		_sum_co_expo_prob = std::vector<MatFlt>(1);
		load_lower_triad_mat(_sum_co_expo_prob[0], file_prefix+"-co_expo_prob_"+file_suffix, _mx_nid);
		}
		#pragma omp section
		{
		_sum_adv_expo_prob = std::vector<MatFlt>(1);
		load_lower_triad_mat(_sum_adv_expo_prob[0], file_prefix+"-adv_expo_prob_"+file_suffix, _mx_nid);
		}
		}
	}

	void assign_treatment_control(const VecFlt& partition, bool use_complete_rand, std::unordered_map<double, int>& assignment, int priority_node=-1) const {
		// cluster_sz[f] is the size of cluster f (number of nodes).
		std::unordered_map<double, int> cluster_sz;
		for (int i = 0; i < _mx_nid; i++) {
			cluster_sz[partition[i]]++;
		}
		std::bernoulli_distribution bern_rv(0.5);
		std::default_random_engine random_eng(rand());
		assignment.clear();
		if (use_complete_rand) {
			if (priority_node != -1) {
				throw std::invalid_argument("Priority Node unsupported for complete randomization!"); 
			}
			std::vector<double> cluster_id;
			cluster_id.reserve(cluster_sz.size());
			for (auto const& p : cluster_sz) {
				cluster_id.push_back(p.first);
			}
			auto size_greater = [&cluster_sz](double c1, double c2) {return cluster_sz.at(c1) > cluster_sz.at(c2);};
			std::sort(cluster_id.begin(), cluster_id.end(), size_greater);

			for (int k = 0; k < cluster_id.size(); k += 2) {
				int a = bern_rv(random_eng);
				assignment[cluster_id[k]] = a;
				if (k + 1 < cluster_id.size()) {
					assignment[cluster_id[k+1]] = 1 - a;
				}
			}
		} else {
			for (auto const& c : cluster_sz) {
				assignment[c.first] = bern_rv(random_eng);
			}
			if (priority_node == -1) return;
			TUNGraph::TNodeI NI = _g->GetNI(priority_node);
			for (int k = 0; k < NI.GetOutDeg(); k++) {
				assignment[partition[NI.GetOutNId(k)]] = 1;
			}
		}
	}

	void simulate_Hajek(const VecFlt& partition, bool use_complete_rand, double& hat_tau_1, double& hat_tau_2) const {
		std::unordered_map<double, int> assignment;
		assign_treatment_control(partition, use_complete_rand, assignment);

		double sum_inv_prob_0 = 0, sum_inv_prob_1 = 0;
		double sum_inv_prob_response_0_1 = 0, sum_inv_prob_response_1_1 = 0;
		double sum_inv_prob_response_0_2 = 0, sum_inv_prob_response_1_2 = 0;
		for (int i = 0; i < _mx_nid; i++) {
			int a = assignment.at(partition[i]);
			TUNGraph::TNodeI NI = _g->GetNI(i);
			bool exposed = true;
			for (int k = 0; k < NI.GetOutDeg(); k++) {
				if (assignment.at(partition[NI.GetOutNId(k)]) != a) {
					exposed = false;
					break;
				}
			}
			if (exposed) {
				if (a) {
					sum_inv_prob_1 += 1 / _sum_expo_prob[0][i];
					sum_inv_prob_response_1_1 += _node_response_1[i] / _sum_expo_prob[0][i] / _normalizer[i];
					sum_inv_prob_response_0_2 += _node_response_0[i] / _sum_expo_prob[0][i] / _normalizer[i];
				} else {
					sum_inv_prob_0 += 1 / _sum_expo_prob[0][i];
					sum_inv_prob_response_0_1 += _node_response_0[i] / _sum_expo_prob[0][i] / _normalizer[i];
					sum_inv_prob_response_1_2 += _node_response_1[i] / _sum_expo_prob[0][i] / _normalizer[i];
				}
			}
		}
		double hat_mu_0_1 = 0, hat_mu_1_1 = 0;
		double hat_mu_0_2 = 0, hat_mu_1_2 = 0;
		if (sum_inv_prob_1 > 0) {
			hat_mu_1_1 = sum_inv_prob_response_1_1 / sum_inv_prob_1;
			hat_mu_0_2 = sum_inv_prob_response_0_2 / sum_inv_prob_1;
		}
		if (sum_inv_prob_0 > 0) {
			hat_mu_0_1 = sum_inv_prob_response_0_1 / sum_inv_prob_0;
			hat_mu_1_2 = sum_inv_prob_response_1_2 / sum_inv_prob_0;
		}
		hat_tau_1 = hat_mu_1_1 - hat_mu_0_1;
		hat_tau_2 = hat_mu_1_2 - hat_mu_0_2;
	}

	void simulate_HT(const VecFlt& partition, bool use_complete_rand, double& hat_tau_1, double& hat_tau_2) const {
		std::unordered_map<double, int> assignment;
		assign_treatment_control(partition, use_complete_rand, assignment);

		double sum_inv_prob_response_0_1 = 0, sum_inv_prob_response_1_1 = 0;
		double sum_inv_prob_response_0_2 = 0, sum_inv_prob_response_1_2 = 0;
		for (int i = 0; i < _mx_nid; i++) {
			int a = assignment.at(partition[i]);
			TUNGraph::TNodeI NI = _g->GetNI(i);
			bool exposed = true;
			for (int k = 0; k < NI.GetOutDeg(); k++) {
				if (assignment.at(partition[NI.GetOutNId(k)]) != a) {
					exposed = false;
					break;
				}
			}
			if (exposed) {
				if (a) {
					sum_inv_prob_response_1_1 += _node_response_1[i] / _sum_expo_prob[0][i];
					sum_inv_prob_response_0_2 += _node_response_0[i] / _sum_expo_prob[0][i];
				} else {
					sum_inv_prob_response_0_1 += _node_response_0[i] / _sum_expo_prob[0][i];
					sum_inv_prob_response_1_2 += _node_response_1[i] / _sum_expo_prob[0][i];
				}
			}
		}
		double hat_mu_0_1 = sum_inv_prob_response_0_1 / _mx_nid;
		double hat_mu_1_1 = sum_inv_prob_response_1_1 / _mx_nid;
		double hat_mu_0_2 = sum_inv_prob_response_0_2 / _mx_nid;
		double hat_mu_1_2 = sum_inv_prob_response_1_2 / _mx_nid;
		hat_tau_1 = hat_mu_1_1 - hat_mu_0_1;
		hat_tau_2 = hat_mu_1_2 - hat_mu_0_2;
	}
};



#endif	// RGCR_h
