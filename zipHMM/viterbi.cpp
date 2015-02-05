#include "viterbi.hpp"
#include "matrix.hpp"
#include "seq_io.hpp"
#include "prob_spaces.hpp"

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cassert>
#include <climits>
#include <deque>

namespace zipHMM {

    double Viterbi::viterbi_seq(const Matrix &pi, const Matrix &A, const Matrix &B,
                                const std::vector<unsigned> &sequence,
                                const Matrix *symbol2matrix) const {
        size_t no_states = A.get_height();
        Matrix res(no_states, 1);

        init_apply_em_log_prob(res, pi, B, sequence[0]);

        Matrix tmp;
        for(size_t c = 1; c < sequence.size(); ++c) {
            Matrix::maxMult<LogSpace>(symbol2matrix[sequence[c]], res, tmp);
            Matrix::copy(tmp, res);
        }

        double path_ll = res(0, 0);
        for(size_t r = 1; r < no_states; ++r) {
            if(res(r, 0) > path_ll) {
                path_ll = res(r, 0);
            }
        }
        return path_ll;
    }


    double Viterbi::viterbi(const Matrix &pi, const Matrix &A, const Matrix &B) const {
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
        for(std::map<size_t, std::vector<std::vector<unsigned> > >::const_iterator it = ds.get_nStates2seqs().begin(); it != ds.get_nStates2seqs().end(); ++it) {
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

        compute_symbol2matrix(symbol2matrix, A, B, alphabet_size);

        double ll = 0.0;
        for(std::vector<std::vector<unsigned> >::const_iterator it = sequences.begin(); it != sequences.end(); ++it) {
            const std::vector<unsigned> &sequence = (*it);
            ll += viterbi_seq(pi, A, B, sequence, symbol2matrix);
        }

        delete[] symbol2matrix;

        return ll;
    }

    void Viterbi::compute_symbol2matrix(Matrix *symbol2matrix,
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
