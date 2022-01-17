#ifndef utils_snap_h
#define utils_snap_h


#include "Snap.h"
#include "utils.h"


// Graph IO
	// Read an undirected graph from a binary file and remove self-loops
	// If such a binary file does not exist, will read from an edgelist.
	PUNGraph read_undirected_graph(const std::string& PGName) {
		PUNGraph G;
		std::ifstream f((PGName + ".ungraph").c_str());
	  if (f.good()) {
	  	TFIn FIn((PGName + ".ungraph").c_str());
	  	G = TUNGraph::Load(FIn);
	  	TSnap::DelSelfEdges(G);
	  } else {
	  	G = TSnap::LoadEdgeList<PUNGraph>((PGName + ".txt").c_str(), 0, 1);
	  	TSnap::DelSelfEdges(G);
		  TFOut out((PGName + ".ungraph").c_str());
		  G->Save(out);
		  out.Flush();
	  }
	  f.close();
	  return G;
	}

// Check
	// Sanity check on graph to ensure node id's are consecutive.
	bool is_node_id_consecutive(PUNGraph G) {
	  for (int i = 0; i < G->GetMxNId(); i ++) {
	    if (!G->IsNode(i)) return false;
	  }
	  return true;
	}

	// Sanity check on graph to ensure (1) node id's are consecutive, (2) graph is connected (3) loopless.
	void check_graph_validity(PUNGraph G) {
	  if (!is_node_id_consecutive(G)) {
	    throw std::invalid_argument("Graph must have consecutive nodeIds [0, 1, ..., n-1]!");
	  }
	  if (!TSnap::IsConnected(G)) {
	    throw std::invalid_argument("Graph must be connected!");
	  }
	  if (TSnap::CntSelfEdges(G)) {
	    throw std::invalid_argument("Graph must not contain self edges!");
	  }
	}

// Degree related
	template <typename GraphType>
	void get_graph_out_deg_seq(GraphType g, VecInt& deg_seq) {
		deg_seq.clear();
		int mx_nid = g->GetMxNId();
		deg_seq.reserve(mx_nid);
		for (int i = 0; i < mx_nid; i++) {
			if (g->IsNode(i)) {
				deg_seq.push_back(g->GetNI(i).GetOutDeg());
			} else {
				deg_seq.push_back(0);
			}
		}
	}

// Counting
	template<typename NIType>
	bool higher_out_degree(const NIType& NI1, const NIType& NI2) {
		int deg_1 = NI1.GetOutDeg(), deg_2 = NI2.GetOutDeg();
		if (deg_1 > deg_2) {
			return true;
		} else if (deg_1 == deg_2 && NI1.GetId() > NI2.GetId()) {
			return true;
		}
		return false;
	}

	void count_triangle(PUNGraph g, VecInt& count) {
	  if (TSnap::CntSelfEdges(g)) {
	    throw std::invalid_argument("Graph must not contain self edges!");
	  }
		count = VecInt(g->GetMxNId());
		for (TUNGraph::TNodeI NI = g->BegNI(); NI != g->EndNI(); NI++) {
			int node_id = NI.GetId();
			std::unordered_set<int> nbr_ids;
			for (int d = 0; d < NI.GetOutDeg(); d++) {
				int nbr_id = NI.GetOutNId(d);
				if (higher_out_degree(NI, g->GetNI(nbr_id))) {
					nbr_ids.insert(nbr_id);
				}
			}
			for (int nbr_id : nbr_ids) {
				TUNGraph::TNodeI nbr_NI = g->GetNI(nbr_id);
				for (int d2 = 0; d2 < nbr_NI.GetOutDeg(); d2++) {
					int nbr_nbr_id = nbr_NI.GetOutNId(d2);
					if (nbr_ids.count(nbr_nbr_id) && nbr_nbr_id > nbr_id) {
						count[node_id]++;
						count[nbr_id]++;
						count[nbr_nbr_id]++;
					}
				}
			}
		}
	}

// Clustering Coefficients
	VecFlt get_clustering_coef(PUNGraph g) {
		VecInt count;
		count_triangle(g, count);
		int n_triangles = 0, n_wedges = 0;
		double avg_clustering = 0;
		double avg_closure = 0;
		for (TUNGraph::TNodeI NI = g->BegNI(); NI != g->EndNI(); NI++) {
			int node_id = NI.GetId();
			int out_deg = NI.GetOutDeg();
			int n_wedges_c = out_deg * (out_deg - 1);
			int n_wedges_h = -out_deg;
			for (int d = 0; d < out_deg; d++) {
				int nbr_id = NI.GetOutNId(d);
				n_wedges_h += g->GetNI(nbr_id).GetOutDeg();
			}
			n_wedges += n_wedges_c;
			n_triangles += count[node_id];
			if (n_wedges_c > 0) {
				avg_clustering += 2.0 * count[node_id] / n_wedges_c;
			}
			if (n_wedges_h > 0) {
				avg_closure += 2.0 * count[node_id] / n_wedges_h;
			}
		}
		avg_clustering /= g->GetNodes();
		avg_closure /= g->GetNodes();
		double gcc = 2.0 * n_triangles / n_wedges;
		return {gcc, avg_clustering, avg_closure};
	}

#endif  // utils_snap_h
