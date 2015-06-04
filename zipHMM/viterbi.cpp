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
            for (size_t c = 0; c < 4; ++c) {
                logB(i, c) = std::log(logB(i, c));
            }
        }
        Matrix viterbi_table;
        Matrix tmp;
        if (compute_path) {
            viterbi_table.reset(no_states, length);
        } else {
            viterbi_table.reset(no_states, 1);
            tmp.reset(no_states, 1);
        }

        //init
        for(size_t r = 0; r < no_states; ++r) {
            viterbi_table(r, 0) = logpi(r, 0) + logB(r, seq[0]);
        }

        // recursion
        for(size_t c = 1; c < length; ++c) {
            for(size_t r = 0; r < no_states; ++r) {
                double max_value = -std::numeric_limits<double>::max();
                if (compute_path) {
                    for(size_t prev_state = 0; prev_state < no_states; ++prev_state) {
                        max_value = std::max(max_value, viterbi_table(prev_state, c - 1) + logA(prev_state, r) + logB(r, seq[c]));
                    }
                    viterbi_table(r,c) = max_value;
                } else {
                    for(size_t prev_state = 0; prev_state < no_states; ++prev_state) {
                        max_value = std::max(max_value, viterbi_table(prev_state, 0) + logA(prev_state, r) + logB(r, seq[c]));
                    }
                    tmp(r, 0) = max_value;
                }
            }

            if (!compute_path) {
                Matrix::copy(tmp, viterbi_table);
            }
        }

        double path_ll;
        size_t end_point;
        if (compute_path) {
            path_ll = viterbi_table(0, length - 1);
            end_point = 0;
            for(size_t r = 1; r < no_states; ++r) {
                if(viterbi_table(r, length - 1) > path_ll) {
                    path_ll = viterbi_table(r, length - 1);
                    end_point = r;
                }
            }
        } else {
            path_ll = viterbi_table(0, 0);
            end_point = 0;
            for(size_t r = 1; r < no_states; ++r) {
                if(viterbi_table(r, 0) > path_ll) {
                    path_ll = viterbi_table(r, 0);
                    end_point = r;
                }
            }
        }

        // backtrack
        if (compute_path) {
            viterbi_path.resize(length);
            viterbi_path[length - 1] = unsigned(end_point);

            size_t current_state = end_point;
            for(size_t c = length - 1; c > 0; --c) {
                double max_value = viterbi_table(0, c - 1) + logA(0, current_state) + logB(current_state, seq[c]);
                size_t max_state = 0;
                for(size_t prev_state = 1; prev_state < no_states; ++prev_state) {
                    double prev_state_ll = viterbi_table(prev_state, c - 1) + logA(prev_state, current_state) + logB(current_state, seq[c]);
                    if(prev_state_ll > max_value) {
                        max_value = prev_state_ll;
                        max_state = prev_state;
                    }
                }
                viterbi_path[c - 1] = unsigned(max_state);
                current_state = max_state;
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
