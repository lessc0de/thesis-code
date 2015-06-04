#include "viterbi.hpp"
#include "matrix.hpp"
#include "seq_io.hpp"
#include "prob_spaces.hpp"
#include "timer.hpp"

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <climits>
#include <utility>
#include <string>
#include <map>
#include <iostream>
#include <cstdlib>
#include <list>

namespace zipHMM {

    double viterbi_helper(const std::vector<unsigned> &seq,
                          const Matrix &pi,
                          const Matrix &A,
                          const Matrix &B,
                          const bool compute_path,
                          std::vector<unsigned> &viterbi_path) {

        size_t no_states = A.get_height();
        size_t length = seq.size();

        // Convert matrices to log space.
        Matrix logpi;
        Matrix::copy(pi, logpi);
        Matrix logA;
        Matrix::copy(A, logA);
        Matrix logB;
        Matrix::copy(B, logB);

        for (size_t i = 0; i < no_states; ++i) {
            logpi(i, 0) = std::log(logpi(i, 0));
        }
        for (size_t i = 0; i < no_states; ++i) {
            for (size_t j = 0; j < no_states; ++j) {
                logA(i, j) = std::log(logA(i, j));
            }
        }
        for (size_t i = 0; i < no_states; ++i) {
            for (size_t c = 0; c < B.get_width(); ++c) {
                logB(i, c) = std::log(logB(i, c));
            }
        }

        // Set up matrices for computation.
        Matrix viterbi_table;
        Matrix tmp(no_states, 1);
        Matrix res(no_states, 1);
        if (compute_path) {
            viterbi_table.reset(no_states, length);
        }

        // Initial state.
        for(size_t r = 0; r < no_states; ++r) {
            viterbi_table(r, 0) = logpi(r, 0) + logB(r, seq[0]);
            res(r, 0) = logpi(r, 0) + logB(r, seq[0]);
        }

        // Recursion.
        for(size_t c = 1; c < length; ++c) {
            for(size_t r = 0; r < no_states; ++r) {
                double max_value = -std::numeric_limits<double>::max();
                size_t max_state = 0;
                for(size_t prev = 0; prev < no_states; ++prev) {
                    double value = res(prev, 0) + logA(prev, r) + logB(r, seq[c]);
                    if (value > max_value) {
                        max_value = value;
                        max_state = prev;
                    }
                }
                tmp(r,0) = max_value;
                if (compute_path) {
                    viterbi_table(r, c) = max_state;
                }
            }
            Matrix::copy(tmp, res);
        }

        // Log-likelihood and end state.
        double path_ll = res(0, 0);
        size_t end_state = 0;
        for(size_t r = 1; r < no_states; ++r) {
            if(res(r, 0) > path_ll) {
                path_ll = res(r, 0);
                end_state = r;
            }
        }

        // Backtracking.
        if (compute_path) {
            viterbi_path.resize(length);
            viterbi_path[length - 1] = end_state;

            for (size_t c = length - 1; c > 0; --c) {
                viterbi_path[c - 1] = viterbi_table(viterbi_path[c], c);
            }
        }
        return path_ll;
    }

    double viterbi(const std::vector<unsigned> &seq,
                   const Matrix &pi,
                   const Matrix &A,
                   const Matrix &B,
                   std::vector<unsigned> &viterbi_path) {
        return viterbi_helper(seq, pi, A, B, true, viterbi_path);
    }

    double viterbi(const std::vector<unsigned> &seq,
                   const Matrix &pi,
                   const Matrix &A,
                   const Matrix &B) {
        std::vector<unsigned> path;
        return viterbi_helper(seq, pi, A, B, false, path);
    }

    double viterbi(const std::string &seq_filename,
                   const Matrix &pi,
                   const Matrix &A,
                   const Matrix &B,
                   std::vector<unsigned> &viterbi_path) {
        std::vector<unsigned> seq;
        readSeq(seq, seq_filename);
        return viterbi(seq, pi, A, B, viterbi_path);
    }

    double viterbi(const std::string &seq_filename,
                   const Matrix &pi,
                   const Matrix &A,
                   const Matrix &B) {
        std::vector<unsigned> seq;
        readSeq(seq, seq_filename);
        return viterbi(seq, pi, A, B);
    }
}
