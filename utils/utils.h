#ifndef utils_h
#define utils_h


#include <iostream>
#include <fstream>
#include <iomanip>

#include <vector>
#include <unordered_map>
#include <map>
#include <unordered_set>
// #include <set>
#include <queue>
#include <string>
#include <sstream>

#include <random>
#include <stdlib.h>		 /* srand, rand */

#include <time.h>			 /* time */
#include <chrono>

#include <cmath>
#include <algorithm>

#include <cassert>

#include "meta.h"




	


// Type definition shorthands
	typedef std::vector<std::string> VecStr;
	typedef std::vector<long> VecInt;
	typedef std::vector<double> VecFlt;
	typedef std::vector<VecInt> MatInt;
	typedef std::vector<VecFlt> MatFlt;

	typedef std::unordered_map<long, long> MapInt;
	typedef std::unordered_map<long, double> MapFlt;
	typedef std::pair<long, long> PairInt;


// Parse string
	VecStr split_string(const std::string& str, const char& delim) {
		VecStr output;
		std::stringstream ss;
		ss.str(str);
		std::string item;
		while (std::getline(ss, item, delim)) {
			output.push_back(item);
		}

		return output;
	}

	VecStr split_string(const std::string& str) {
		VecStr output;

		std::stringstream ss;
		ss.str(str);
		std::string item;
		while (ss >> item) {
			output.push_back(item);
		}
		return output;
	}

	void parse_path_file_name(const std::string& path_file_name, std::string& path, std::string& filename) {
		size_t posSlash = path_file_name.rfind('/');
		if (posSlash == std::string::npos) {
			filename = path_file_name;
			path = "";
		} else {
			path = path_file_name.substr(0, posSlash + 1);
			filename = path_file_name.substr(posSlash + 1);
		}
	}


// Print ADT
	template <typename T>
	void print_vector(const std::vector<T>& vec, std::ostream& out = std::cout, char delim = ' ') {
		for (int i = 0; i < vec.size(); i ++) {
			out << vec[i] << delim;
		}
		out << std::endl;
	}	

	template <typename T1, typename T2>
	void print_map(const std::unordered_map<T1, T2>& mapping, std::ostream& out = std::cout, std::string delim = ", ") {
		out << "A hashmap of size " << mapping.size() << ": [ ";
		for (auto const& p : mapping) {
			out << p.first << ":" << p.second << delim;
		}
		out << "]" << std::endl;
	}	

	template <typename T>
	void print_set(const std::unordered_set<T>& set, std::ostream& out = std::cout, char delim = ',') {
		out << "A hashset of size " << set.size() << ": {";
		for (auto const& n : set) {
			out << n << delim;
		}
		out << '}' << std::endl;
	}

	template <typename T>
	void print_matrix(const std::vector<std::vector<T> >& mat, std::ostream& out = std::cout) {
		for (size_t i = 0; i < mat.size(); i ++) {
			print_vector(mat[i], out, '\t');
		}
	}

	void load_vec_from_file(VecFlt& vec, const std::string& file_name, const int col, const int skip_row = 0, const int expected_size = 0) {
		vec.clear();
		if (expected_size > 0) {
			vec.reserve(expected_size);
		}

		std::ifstream file_input(file_name);
		std::string line;
		for (int i = 0; i < skip_row; i++) {
			std::getline(file_input, line);
		}

		while (std::getline(file_input, line)) {
			VecStr elements = split_string(line, '\t');
			assert(elements.size() > col);
			vec.push_back(std::stod(elements[col]));
		}
		file_input.close();

		if (expected_size > 0 && vec.size() != expected_size) {
			std::cout << "Unexpected size of loaded vector!!" << std::endl;
			std::cout << vec.size() << ' ' << expected_size << ' ' << file_name << std::endl;
		}
	}

	void load_lower_triad_mat(MatFlt& mat, const std::string& file_name, int expected_size = 0) {
		mat.clear();
		if (expected_size > 0) {
			mat.reserve(expected_size);
		}

		std::ifstream file_input(file_name);
		std::string line;
		while (std::getline(file_input, line)) {
			mat.push_back(VecFlt(0));
			VecFlt& row = mat.back();
			row.reserve(mat.size() - 1);
			VecStr line_vec = split_string(line, '\t');
			if (line_vec.size() > 0 && line_vec.back() == "") {
				line_vec.pop_back();
			}
			if (line_vec.size() != mat.size() - 1) {
				std::cout << "Unexpected size of loaded matrix row!!" << std::endl;
				std::cout << mat.size() << ' ' << line_vec.size() << ' ' << line << std::endl;
			}
			for (const std::string& s : line_vec) {
				row.push_back(std::stod(s));
			}
		}
		if (expected_size > 0 && mat.size() != expected_size) {
			std::cout << "Unexpected size of loaded matrix!!" << std::endl;
			std::cout << mat.size() << ' ' << expected_size << ' ' << file_name << std::endl;
		}
		file_input.close();
	}


// Random number generation
	// To reset the seed for rng
	void reset_random_seed() {
		std::srand(std::time(NULL));
		int num_rounds = rand()%5;
		for (int i = 0; i < num_rounds; i ++) {
			rand();
		}
	}

	// To generate a random floating-point number in [LB, UB)
	double rand_double(double LB = 0, double UB = 1) {
		return LB + rand() / ((double)RAND_MAX+1) * (UB - LB);
	}

	void gen_beta_rnd_vector(VecFlt& output, const VecFlt& weight, const int vec_size) {
		output.clear();
		output.reserve(vec_size);

		if (weight.empty()) {
			for (int i = 0; i < vec_size; i ++) {
				output.push_back( rand_double() );
			}			
		} else {
			assert(weight.size() == vec_size);
			for (int i = 0; i < vec_size; i ++) {
				output.push_back( std::log(rand_double()) / weight[i]);
			}
		}
	}

	/**
	 * @brief	A function to generate a random ordering of numbers {0, 1, ..., n-1}.
	 *
	 * @param[out] ordering	A placeholder to store the random ordering output.
	 * @param[in]	n				 The number of consecutive numbers in the random ordering output.
	 * @param[in]	weights	 If nonempty, then Pr[i comes before j] = _weights_[i] / (_weights_[i] + _weights_[j]),
	 *		 i.e., this variable adjusts which number is more likely to be in the front.
	 * @param[in]	priority_i	If specified, then the number _priority_i_ is guaranteed to be the first. 
	 *
	 * @pre	_n_ >= 0.
	 * @pre	If nonempty, then Size of _weights_ == n.
	 * @pre	If specified, then 0 <= _priority_i_ < n.
	 *
	 * @post	ordering[k] is the number at the k-th order.
	 *
	 * Complexity: if _weight_ is empty, then O(n); otherwise, O(n(log n)).
	 */
	void gen_rnd_ordering(VecInt& ordering, const int n, const VecFlt& weights, const int priority_i = -1) {
		ordering.clear();
		ordering.reserve(n);
		for (int k = 0; k < n; k ++) {
			ordering.push_back(k);
		}
		if (weights.empty()) {
			if (priority_i != -1) {
				assert(priority_i >= 0 && priority_i < n);
				ordering[0] = priority_i;
				ordering[priority_i] = 0;
				std::shuffle(ordering.begin()+1, ordering.end(), std::default_random_engine(rand()));
			} else {
				std::shuffle(ordering.begin(), ordering.end(), std::default_random_engine(rand()));
			}
		} else {
			VecFlt rnd_beta_num;
			gen_beta_rnd_vector(rnd_beta_num, weights, n);
			if (priority_i != -1) {
				rnd_beta_num[priority_i] = 1;
			}
			auto num_greater = [&rnd_beta_num](int i, int j) {return rnd_beta_num[i] > rnd_beta_num[j];};
			std::sort(ordering.begin(), ordering.end(), num_greater);
		}
	}

	void get_rnd_N01(VecFlt& output, const int vec_size) {
		std::default_random_engine generator;
		std::normal_distribution<double> distribution(0, 1);
		output.clear();
		output.reserve(vec_size);
		for (int i = 0; i < vec_size; i ++) {
			output.push_back(distribution(generator));
		}
	}


// Process ADT
	template <typename T1, typename T2>
	void add_vector(std::vector<T1>& addTo, const std::vector<T2>& toAdd) {
		if (addTo.size() != toAdd.size()) {
			std::cout << "Vector size does not match!!" << std::endl;
			assert(false);
		}
		for (int i = 0; i < addTo.size(); i ++) {
			if (!std::isnan(toAdd[i])) {
				addTo[i] += toAdd[i];
			}
		}
	}

	template <typename T1, typename T2>
	void add_matrix(std::vector< std::vector<T1> >& addTo, const std::vector< std::vector<T2> >& toAdd) {
		if (addTo.size() != toAdd.size()) {
			std::cout << "Number of rows does not match!!" << std::endl;
			assert(false);
		}
		for (int i = 0; i < addTo.size(); i ++) {
			add_vector(addTo[i], toAdd[i]);
		}
	}

	template <typename T>
	void get_vector_prefix_sum(const std::vector<T>& vec, std::vector<T>& prefix_sum) {
		prefix_sum.clear();
		prefix_sum.reserve(vec.size());
		if (vec.empty()) return;
		prefix_sum.push_back(vec[0]);
		for (int i = 1; i < vec.size(); i++) {
			prefix_sum.push_back(prefix_sum.back() + vec[i]);
		}
	}

	template <typename T>
	bool is_vector_sorted(const std::vector<T>& vec) {
		if (vec.empty()) return true;
		for (int i = 1; i < vec.size(); i++) {
			if (vec[i-1] > vec[i]) return false;
		}
		return true;
	}

	template <typename T>
	size_t arg_upper_bound(const std::vector<T>& vec, int val) {
		assert(is_vector_sorted(vec));
		if (vec.empty()) return 0;
		if (vec[0] > val) return 0;
		if (vec.back() <= val) return vec.size();
		size_t lb = 0, ub = vec.size() - 1;
		while (ub - lb > 1) {
			size_t m = lb + (ub - lb) / 2;
			if (vec[m] > val) {
				ub = m;
			}
			else {
				lb = m;
			}
		}
		return ub;
	}

	template <typename T>
	T get_max_val(const std::vector<T>& vec) {
		if (vec.size() == 0) {
			throw std::invalid_argument("get_max_val: vector must be nonempty!");
		}
		T max_val = vec[0];
		for (int i = 0; i < vec.size(); i++) {
			max_val = std::max(max_val, vec[i]);
		}
		return max_val;
	}

	template <typename T>
	T get_min_val(const std::vector<T>& vec) {
		if (vec.size() == 0) {
			throw std::invalid_argument("get_min_val: vector must be nonempty!");
		}
		T min_val = vec[0];
		for (int i = 0; i < vec.size(); i++) {
			min_val = std::min(min_val, vec[i]);
		}
		return min_val;
	}

	template <typename T>
	T get_max_abs_val(const std::vector<T>& vec) {
		T max_val = get_max_val(vec);
		T min_val = get_min_val(vec);
		return std::max(max_val, -min_val);
	}


// Other MISC
	std::string get_time_str() {
		std::time_t t1 = std::time(nullptr);
		std::string time_str = std::asctime(std::localtime(&t1));
		time_str.pop_back();
		return time_str;
	}



#endif	// utils_h
