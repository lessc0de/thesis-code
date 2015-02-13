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

    double Viterbi::viterbi_seq(const Matrix &pi, const Matrix &A, const Matrix &B,
                                const std::vector<unsigned> &sequence,
                                const Matrix *symbol2matrix,
                                const Matrix *symbol2argmax_matrix,
                                const bool compute_path,
                                std::vector<unsigned> &viterbi_path) const {
        size_t no_states = A.get_height();
        Matrix res(no_states, 1);

        // Save paths in this table.
        Matrix viterbi_table(sequence.size(), no_states);

        // Init.
        init_apply_em_log_prob(res, pi, B, sequence[0]);

        // Recursion.
        Matrix tmp;
        for(size_t c = 1; c < sequence.size(); ++c) {
            // Compute log likelihood for each cell in column.
            Matrix::maxMult<LogSpace>(symbol2matrix[sequence[c]], res, tmp);

            // Compute path for each cell in column.
            if (compute_path) {
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
            }
            Matrix::copy(tmp, res);
        }

        // Find the end point with the largest log likelihood.
        double path_ll = res(0, 0);
        size_t end_point = 0;
        for(size_t r = 1; r < no_states; ++r) {
            if(res(r, 0) > path_ll) {
                path_ll = res(r, 0);
                end_point = r;
            }
        }

        // Backtrack.
        if (compute_path) {
            viterbi_path.resize(sequence.size());
            viterbi_path[sequence.size() - 1] = unsigned(end_point);

            for (size_t c = sequence.size() - 1; c > 0; --c) {
                viterbi_path[c - 1] = viterbi_table(c, viterbi_path[c]);
            }

            // Recreate the original Viterbi path.
            deduct_path(sequence, viterbi_path, symbol2argmax_matrix);
        }

        return path_ll;
    }

    void Viterbi::deduct_path(const std::vector<unsigned> &sequence,
                              std::vector<unsigned> &path,
                              const Matrix *symbol2argmax_matrix) const {
        // Convert the vectors to lists.
        std::list<unsigned> orig_seq(sequence.begin(), sequence.end());
        std::list<unsigned> orig_path(path.begin(), path.end());

        // Iterate through the list. Insert/delete symbols.
        std::map<unsigned, s_pair> symbol2pair = ds.get_symbol2pair();
        std::list<unsigned>::iterator seq_it = orig_seq.begin();
        std::list<unsigned>::iterator path_it = orig_path.begin();
        while (seq_it != orig_seq.end()) {
            // Check if character at this position is an extended character.
            if (*seq_it >= ds.get_orig_alphabet_size()) {
                // Insert the state in the list.
                --path_it;
                unsigned prev_state = *path_it;
                ++path_it;
                unsigned next_state = *path_it;
                unsigned current_state = symbol2argmax_matrix[*seq_it](prev_state, next_state);
                orig_path.insert(path_it, current_state);
                --path_it;

                // Insert the symbol in the list.
                std::pair<size_t, size_t> p = symbol2pair[*seq_it];
                orig_seq.insert(seq_it, p.first);
                orig_seq.insert(seq_it, p.second);
                orig_seq.erase(seq_it);
                --seq_it;
                --seq_it;
            } else {
                ++seq_it;
                ++path_it;
            }
        }

        // Convert the list to vectors.
        path.clear();
        path.insert(path.begin(), orig_path.begin(), orig_path.end());

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

        zipHMM::Timer stage_1_timer;
        zipHMM::Timer stage_2_timer;
        double stage_1_time = 0;
        double stage_2_time = 0;
        stage_1_timer.start();
        Matrix *symbol2matrix = new Matrix[alphabet_size];
        Matrix *symbol2argmax_matrix = new Matrix[alphabet_size];

        compute_symbol2matrix(symbol2matrix, symbol2argmax_matrix, A, B, alphabet_size);

        stage_1_timer.stop();
        stage_1_time += stage_1_timer.timeElapsed();

        stage_2_timer.start();
        double ll = 0.0;
        for(std::vector<std::vector<unsigned> >::const_iterator it = sequences.begin(); it != sequences.end(); ++it) {
            const std::vector<unsigned> &sequence = (*it);
            ll += viterbi_seq(pi, A, B, sequence, symbol2matrix, symbol2argmax_matrix, compute_path, viterbi_path);
        }
        stage_2_timer.stop();
        stage_2_time += stage_2_timer.timeElapsed();

        stage_1_timer.start();

        delete[] symbol2matrix;
        delete[] symbol2argmax_matrix;

        stage_1_timer.stop();
        stage_1_time += stage_1_timer.timeElapsed();

        std::cout << stage_1_time << " ";
        std::cout.flush();
        std::cout << stage_2_time << " ";
        std::cout.flush();

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
        std::map<unsigned, s_pair> symbol2pair = ds.get_symbol2pair();
        for(size_t i = ds.get_orig_alphabet_size(); i < alphabet_size; ++i) {
            const s_pair symbol_pair = symbol2pair.find(unsigned(i))->second;
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
