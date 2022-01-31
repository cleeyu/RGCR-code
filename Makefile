include Makefile.inc

default: all

all: estimate_rgcr_expo_prob compute_rgcr_variance simulate_rgcr_variance simulate_gcr_variance compute_ball_size

estimate_rgcr_expo_prob: estimate_rgcr_expo_prob.cpp RGCR.h random_clustering/random_clustering.h utils/utils.h utils/utils_snap.h
	$(CXX) $(CFLAGS) $< $(LDFLAGS) $(SNAP_OBJ) -o estimate_rgcr_expo_prob

compute_rgcr_variance: compute_rgcr_variance.cpp RGCR.h random_clustering/random_clustering.h utils/utils.h utils/utils_snap.h
	$(CXX) $(CFLAGS) $< $(LDFLAGS) $(SNAP_OBJ) -o compute_rgcr_variance

simulate_rgcr_variance: simulate_rgcr_variance.cpp RGCR.h random_clustering/random_clustering.h utils/utils.h utils/utils_snap.h
	$(CXX) $(CFLAGS) $< $(LDFLAGS) $(SNAP_OBJ) -o simulate_rgcr_variance

simulate_gcr_variance: simulate_gcr_variance.cpp RGCR.h random_clustering/random_clustering.h utils/utils.h utils/utils_snap.h
	$(CXX) $(CFLAGS) $< $(LDFLAGS) $(SNAP_OBJ) -o simulate_gcr_variance

compute_ball_size: compute_ball_size.cpp ball_size.h utils/utils.h utils/utils_snap.h
	$(CXX) $(CFLAGS) $< $(LDFLAGS) $(SNAP_OBJ) -o compute_ball_size


.PHONY: clean
clean:
	rm -rf *.o *~ estimate_rgcr_expo_prob compute_rgcr_variance simulate_rgcr_variance simulate_gcr_variance compute_ball_size
