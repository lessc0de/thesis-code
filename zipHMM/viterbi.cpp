#include "viterbi.hpp"
#include "matrix.hpp"
#include "seq_io.hpp"
#include "prob_spaces.hpp"

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

namespace zipHMM {

    double Viterbi::viterbi_seq(const Matrix &pi, const Matrix &A, const Matrix &B,
                                const std::vector<unsigned> &sequence,
                                const Matrix *symbol2matrix,
                                const Matrix *symbol2argmax_matrix,
                                const bool compute_path,
                                std::vector<unsigned> &viterbi_path) const {
        size_t no_states = A.get_height();
        Matrix res(no_states, 1);

        Matrix viterbi_table(sequence.size(), no_states);

        init_apply_em_log_prob(res, pi, B, sequence[0]);

        // // Copy to Viterbi table.
        // for (size_t i = 0; i < no_states; ++i) {
        //     viterbi_table(i, 0) = res(i, 0);
        // }

        Matrix tmp;
        for(size_t c = 1; c < sequence.size(); ++c) {
            Matrix::maxMult<LogSpace>(symbol2matrix[sequence[c]], res, tmp);

            // For each state i, from which state j did we come from?
            Matrix where;
            where.reset(symbol2matrix[sequence[c]].get_height(), res.get_width());
            for(size_t l = 0; l < symbol2matrix[sequence[c]].get_height(); ++l) {
                for(size_t r = 0; r < res.get_width(); ++r) {
                    // Find the best k
                    size_t k = Matrix::argMaxMult<LogSpace>(symbol2matrix[sequence[c]], l, res, r);
                    where(l, r) = k;
                }
            }

            for (size_t i = 0; i < no_states; ++i) {
                viterbi_table(c, i) = where(i, 0);
            }

            Matrix::copy(tmp, res);
        }

        // double path_ll = res(0, 0);
        // for(size_t r = 1; r < no_states; ++r) {
        //     if(res(r, 0) > path_ll) {
        //         path_ll = res(r, 0);
        //     }
        // }
        // return path_ll;

        // backtrack
        double path_ll = res(0, 0);
        size_t end_point = 0;
        for(size_t r = 1; r < no_states; ++r) {
            if(res(r, 0) > path_ll) {
                path_ll = res(r, 0);
                end_point = r;
            }
        }

        if (compute_path) {
            viterbi_path.resize(sequence.size());
            viterbi_path[sequence.size() - 1] = unsigned(end_point);

            for (size_t c = sequence.size() - 1; c > 0; --c) {
                viterbi_path[c - 1] = viterbi_table(c, viterbi_path[c]);
            }

            // Try to recreate the original sequence and Viterbi path.
            std::vector<unsigned> orig_seq = sequence;
            std::vector<unsigned> orig_path = viterbi_path;
            while (orig_seq.size() < ds.get_orig_seq_length()) {
                std::pair<std::vector<unsigned>, std::vector<unsigned> > p;
                p = deducted_path(orig_seq, orig_path, symbol2argmax_matrix);
                orig_seq = p.first;
                orig_path = p.second;
            }
            viterbi_path = orig_path;
        }

        return path_ll;
    }

    std::pair<std::vector<unsigned>, std::vector<unsigned> > Viterbi::deducted_path(const std::vector<unsigned> &sequence,
                                                                                    const std::vector<unsigned> &path,
                                                                                    const Matrix *symbol2argmax_matrix) const {
        std::vector<unsigned> orig_sequence;
        std::vector<unsigned> path_sequence;

        orig_sequence.push_back(sequence[0]);
        path_sequence.push_back(path[0]);

        for (size_t c = 1; c < sequence.size(); ++c) {
            if (sequence[c] >= ds.get_orig_alphabet_size()) {
                // This is a character added during compression.
                std::pair<size_t, size_t> p = ds.get_symbol2pair()[sequence[c]];
                orig_sequence.push_back(p.first);
                orig_sequence.push_back(p.second);

                unsigned prev_state = path[c - 1];
                unsigned next_state = path[c];
                path_sequence.push_back(symbol2argmax_matrix[sequence[c]](prev_state, next_state));
                path_sequence.push_back(path[c]);
            } else {
                orig_sequence.push_back(sequence[c]);
                path_sequence.push_back(path[c]);
            }
        }
        return std::make_pair(orig_sequence, path_sequence);
    }


    double Viterbi::viterbi_helper(const Matrix &pi, const Matrix &A, const Matrix &B, const bool compute_path, std::vector<unsigned> &viterbi_path) const {
        if(pi.get_width() != 1 || pi.get_height() != A.get_width() || A.get_height() != A.get_width() ||
           B.get_height() != A.get_width() || B.get_width() != ds.get_orig_alphabet_size()) {
            std::cerr << "Dimensions of input matrices do not match:" << std::endl;
            std::cerr << "\t" << "pi width:  " << pi.get_width()  << std::endl;
            std::cerr << "\t" << "pi height: " << pi.get_height() << std::endl;
            std::cerr << "\t" << "A width:  "  << A.get_width()   << std::endl;
            std::cerr << "\t" << "A height: "  << A.get_height()  << std::endl;
            std::cerr << "\t" << "B width:  "  << B.get_width()   << std::endl;
            std::cerr << "\t" << "B height: "  << B.get_height()  << std::endl;
            std::exit(-1);
        }

        // find alphabet and seqs for given number of states
        size_t no_states = A.get_width();
        size_t alphabet_size = 0;
        std::vector<std::vector< unsigned> > sequences;
        for(std::map<size_t, std::vector<std::vector<unsigned> > >::const_iterator it = ds.nStates2seqs.begin(); it != ds.nStates2seqs.end(); ++it) {
            if(it->first >= no_states) {
                sequences = it->second;
                alphabet_size = ds.get_nStates2alphabet_size().find(it->first)->second;
                break;
            }
        }
        if(sequences.empty()) {
            sequences = ds.get_nStates2seqs().rbegin()->second;
            alphabet_size = ds.get_nStates2alphabet_size().rbegin()->second;
        }

        Matrix *symbol2matrix = new Matrix[alphabet_size];
        Matrix *symbol2argmax_matrix = new Matrix[alphabet_size];

        compute_symbol2matrix(symbol2matrix, symbol2argmax_matrix, A, B, alphabet_size);

        double ll = 0.0;
        for(std::vector<std::vector<unsigned> >::const_iterator it = sequences.begin(); it != sequences.end(); ++it) {
            const std::vector<unsigned> &sequence = (*it);
            ll += viterbi_seq(pi, A, B, sequence, symbol2matrix, symbol2argmax_matrix, compute_path, viterbi_path);
        }

        delete[] symbol2matrix;
        delete[] symbol2argmax_matrix;

        return ll;
    }

    double Viterbi::viterbi(const Matrix &pi, const Matrix &A, const Matrix &B) const {
        std::vector<unsigned> v;
        return Viterbi::viterbi_helper(pi, A, B, false, v);
    }

    double Viterbi::viterbi(const Matrix &pi, const Matrix &A, const Matrix &B,
                            std::vector<unsigned> &viterbi_path) const {
        return Viterbi::viterbi_helper(pi, A, B, true, viterbi_path);
    }

    void Viterbi::compute_symbol2matrix(Matrix *symbol2matrix,
                                        Matrix *symbol2argmax_matrix,
                                        const Matrix &A,
                                        const Matrix &B,
                                        const size_t alphabet_size) const{
        // compute C matrices for each symbol in the original alphabet
        make_em_trans_log_probs_array(symbol2matrix, A, B);

        // compute C matrices for each symbol in the extended alphabet
        for(size_t i = ds.get_orig_alphabet_size(); i < alphabet_size; ++i) {
            const s_pair symbol_pair = ds.get_symbol2pair().find(unsigned(i))->second;
            const unsigned left_symbol = symbol_pair.second; // the multiplication is done in the reverse direction of the sequence
            const unsigned right_symbol  = symbol_pair.first;
            Matrix &left_matrix  = symbol2matrix[left_symbol];
            Matrix &right_matrix = symbol2matrix[right_symbol];

            Matrix::maxMult<LogSpace>(left_matrix, right_matrix, symbol2matrix[i]);

            for(size_t l = 0; l < left_matrix.get_height(); ++l) {
                for(size_t r = 0; r < right_matrix.get_width(); ++r) {
                    // Find the best k
                    symbol2argmax_matrix[i].reset(left_matrix.get_height(), right_matrix.get_width());
                    size_t k = Matrix::argMaxMult<LogSpace>(left_matrix, l, right_matrix, r);
                    symbol2argmax_matrix[i](r, l) = k;
                }
            }
        }
    }



    double viterbi(const std::vector<unsigned> &seq,
                   const Matrix &pi,
                   const Matrix &A,
                   const Matrix &B,
                   std::vector<unsigned> &viterbi_path) {

        size_t no_states = A.get_height();
        size_t length = seq.size();

        Matrix viterbi_table(no_states, length);

        //init
        for(size_t r = 0; r < no_states; ++r) {
            viterbi_table(r, 0) = std::log(pi(r, 0) * B(r, seq[0]));
        }

        // recursion
        for(size_t c = 1; c < length; ++c) {
            for(size_t r = 0; r < no_states; ++r) {
                double max_value = -std::numeric_limits<double>::max();
                for(size_t prev_state = 0; prev_state < no_states; ++prev_state) {
                    max_value = std::max(max_value, viterbi_table(prev_state, c - 1) + std::log(A(prev_state, r) * B(r, seq[c])));
                }
                viterbi_table(r,c) = max_value;
            }
        }

        double path_ll = viterbi_table(0, length - 1);
        size_t end_point = 0;
        for(size_t r = 1; r < no_states; ++r) {
            if(viterbi_table(r, length - 1) > path_ll) {
                path_ll = viterbi_table(r, length - 1);
                end_point = r;
            }
        }

        // backtrack
        viterbi_path.resize(length);
        viterbi_path[length - 1] = unsigned(end_point);

        size_t current_state = end_point;
        for(size_t c = length - 1; c > 0; --c) {
            double max_value = viterbi_table(0, c - 1) + std::log(A(0, current_state) * B(current_state, seq[c]));
            size_t max_state = 0;
            for(size_t prev_state = 1; prev_state < no_states; ++prev_state) {
                double prev_state_ll = viterbi_table(prev_state, c - 1) + std::log(A(prev_state, current_state) * B(current_state, seq[c]));
                if(prev_state_ll > max_value) {
                    max_value = prev_state_ll;
                    max_state = prev_state;
                }
            }
            viterbi_path[c - 1] = unsigned(max_state);
            current_state = max_state;
        }

        return path_ll;
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
}
