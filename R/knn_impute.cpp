// knn_impute.cpp — Rcpp KNN imputer matching MATLAB's knnimpute behaviour
//
// Input:  X  = n x p numeric matrix (rows = samples, cols = features)
//         k  = number of nearest neighbours
//
// Algorithm (matching MATLAB's knnimpute):
//   1. Find globally-complete features (no NA in any sample).
//   2. Compute full n x n Euclidean distance matrix on those features.
//   3. Per sample: sort all other samples by distance (precomputed once).
//   4. For each NA at X(i,j):
//      a. Take k nearest neighbours (include ties at k-th boundary).
//      b. Among those: use ones with valid value at j → 1/d weighted mean.
//      c. If none valid: 1-NN fallback — scan sorted list beyond k for
//         first sample with valid value at j → use its value.
//      d. If no sample anywhere has a valid value: feature mean.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix knn_impute_rcpp(NumericMatrix X, int k) {
    int n = X.nrow();
    int p = X.ncol();

    // --- 1. NA mask (row-major for cache efficiency) ----------------------
    std::vector<bool> na_mask(n * p);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j)
            na_mask[i * p + j] = NumericMatrix::is_na(X(i, j));

    // --- 2. Find globally-complete features -------------------------------
    std::vector<int> complete_cols;
    for (int j = 0; j < p; ++j) {
        bool complete = true;
        for (int i = 0; i < n; ++i) {
            if (na_mask[i * p + j]) { complete = false; break; }
        }
        if (complete) complete_cols.push_back(j);
    }
    if (complete_cols.empty()) {
        Rcpp::stop("All rows of the input data contain missing values. "
                    "Unable to compute distances.");
    }
    int d = (int)complete_cols.size();

    // --- 3. Full n x n Euclidean distance matrix --------------------------
    std::vector<double> dist(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int r = i + 1; r < n; ++r) {
            double sum_sq = 0.0;
            for (int ci = 0; ci < d; ++ci) {
                int col = complete_cols[ci];
                double diff = X(i, col) - X(r, col);
                sum_sq += diff * diff;
            }
            double d_val = std::sqrt(sum_sq);
            dist[i * n + r] = d_val;
            dist[r * n + i] = d_val;
        }
    }

    // --- 4. Per-sample sorted neighbour lists (precomputed once) ----------
    std::vector<std::vector<std::pair<double, int>>> neighbors(n);
    for (int i = 0; i < n; ++i) {
        neighbors[i].reserve(n - 1);
        for (int r = 0; r < n; ++r) {
            if (r == i) continue;
            neighbors[i].push_back(std::make_pair(dist[i * n + r], r));
        }
        std::sort(neighbors[i].begin(), neighbors[i].end());
    }

    // --- 5. Feature means (ultimate fallback) -----------------------------
    std::vector<double> col_mean(p, 0.0);
    for (int j = 0; j < p; ++j) {
        double sum = 0.0;
        int cnt = 0;
        for (int i = 0; i < n; ++i) {
            if (!na_mask[i * p + j]) { sum += X(i, j); ++cnt; }
        }
        col_mean[j] = (cnt > 0) ? sum / cnt : 0.0;
    }

    // --- 6. Deep-copy result (do not mutate X) ----------------------------
    NumericMatrix result(n, p);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j)
            result(i, j) = X(i, j);

    // --- 7. Impute --------------------------------------------------------
    const double MINDIST = 1e-6;

    for (int i = 0; i < n; ++i) {
        const auto& nbrs = neighbors[i];
        int nn = (int)nbrs.size(); // n - 1

        for (int j = 0; j < p; ++j) {
            if (!na_mask[i * p + j]) continue;

            // k_actual neighbours, extended to include ties at k-th boundary
            int k_actual = std::min(k, nn);
            double kth_dist = nbrs[k_actual - 1].first;
            int k_with_ties = k_actual;
            while (k_with_ties < nn &&
                   nbrs[k_with_ties].first <= kth_dist) {
                ++k_with_ties;
            }

            // Tier 1: among k (+ ties), use valid ones → 1/d weighted mean
            double w_sum = 0.0;
            double val_sum = 0.0;
            for (int c = 0; c < k_with_ties; ++c) {
                int r = nbrs[c].second;
                if (na_mask[r * p + j]) continue;
                double dv = nbrs[c].first;
                if (dv < MINDIST) dv = MINDIST;
                double w = 1.0 / dv;
                val_sum += w * X(r, j);
                w_sum += w;
            }
            if (w_sum > 0.0) {
                result(i, j) = val_sum / w_sum;
                continue;
            }

            // Tier 2: 1-NN fallback — scan beyond k for first valid neighbour
            bool found = false;
            for (int c = k_with_ties; c < nn; ++c) {
                int r = nbrs[c].second;
                if (!na_mask[r * p + j]) {
                    result(i, j) = X(r, j);
                    found = true;
                    break;
                }
            }
            if (found) continue;

            // Tier 3: no sample has a valid value → feature mean
            result(i, j) = col_mean[j];
        }
    }

    return result;
}
