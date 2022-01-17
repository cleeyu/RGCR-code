include Makefile.inc

default: all

all: simulate_rgcr eval_prob_formula eval_prob_simulation simulate_GCR_variance

simulate_rgcr: simulate_rgcr.cpp RGCR.h random_clustering/random_clustering.h utils/utils.h utils/utils_snap.h
	$(CXX) $(CFLAGS) $< $(LDFLAGS) $(SNAP_OBJ) -o simulate_rgcr

eval_prob_formula: eval_prob_formula.cpp RGCR.h random_clustering/random_clustering.h utils/utils.h utils/utils_snap.h
	$(CXX) $(CFLAGS) $< $(LDFLAGS) $(SNAP_OBJ) -o eval_prob_formula

eval_prob_simulation: eval_prob_simulation.cpp RGCR.h random_clustering/random_clustering.h utils/utils.h utils/utils_snap.h
	$(CXX) $(CFLAGS) $< $(LDFLAGS) $(SNAP_OBJ) -o eval_prob_simulation

simulate_GCR_variance: simulate_GCR_variance.cpp RGCR.h random_clustering/random_clustering.h utils/utils.h utils/utils_snap.h
	$(CXX) $(CFLAGS) $< $(LDFLAGS) $(SNAP_OBJ) -o simulate_GCR_variance

.PHONY: clean
clean:
	rm -rf *.o *~ simulate_rgcr eval_prob_formula eval_prob_simulation simulate_GCR_variance
