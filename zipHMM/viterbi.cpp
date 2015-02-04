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

    void Viterbi::init_data_structure(const std::vector<std::vector<unsigned> > &sequences,
                                      const size_t alphabet_size,
                                      std::vector<size_t> &nStatesSave,
                                      const size_t min_no_eval) {
        orig_alphabet_size = alphabet_size;
        std::vector<size_t> pair_n;
        bool counted = false;
        size_t prev_max_count = UINT_MAX;
        Timer iteration_timer;
        bool nStatesSave_init_empty = nStatesSave.empty();
        size_t no_states_save = 0;
        std::map<size_t, double> nStates2t_MM;
        std::map<size_t, double> nStates2t_MV;
        double iteration_time;

        if(!nStatesSave_init_empty) {
            no_states_save = nStatesSave.back();
            for(std::vector<size_t>::iterator it = nStatesSave.begin(); it != nStatesSave.end(); ++it) {
                nStates2t_MM[*it] = Matrix::time_blas_mult((*it));
                nStates2t_MV[*it] = Matrix::time_blas_matrix_vector_mult((*it));
            }
        }

        orig_seq_length = 0;
        for(std::vector<std::vector<unsigned> >::const_iterator it = sequences.begin(); it != sequences.end(); ++it) {
            orig_seq_length += it->size();
        }

        std::vector<std::vector<unsigned> > *prev_seqs_p = new std::vector<std::vector<unsigned> >;
        std::vector<std::vector<unsigned> > *current_seqs_p = new std::vector<std::vector<unsigned> >;
        std::vector<std::vector<unsigned> > *next_seqs_p = 0;

        prev_seqs_p->assign(sequences.begin(), sequences.end()); // copy sequences - no way out of this
        current_seqs_p->assign(sequences.begin(), sequences.end()); // copy sequences - no way out of this

        // weird special cases I need to handle
        if(orig_alphabet_size == 1 || no_states_save == 1) {
            if(nStatesSave_init_empty) {
                nStates2seqs[2] = *current_seqs_p; // copy seq
                nStates2alphabet_size[2] = orig_alphabet_size;
                return;
            }
            for(std::vector<size_t>::iterator it = nStatesSave.begin(); it != nStatesSave.end(); ++it) {
                nStates2seqs[*it] = *current_seqs_p; // copy seq
                nStates2alphabet_size[*it] = orig_alphabet_size;
            }
            return;
        }

        iteration_timer.start();

        // count each pair of symbols across the original sequence
        pair_n.resize(orig_alphabet_size * orig_alphabet_size);
        for(std::vector<std::vector<unsigned> >::iterator it = current_seqs_p->begin(); it != current_seqs_p->end(); ++it) {
            std::vector<unsigned> &current_seq = (*it);
            update_pair_n(pair_n, s_pair( current_seq[1] , current_seq[2] ), orig_alphabet_size);
            counted = true;
            for(size_t i = 2; i < current_seq.size() - 1; ++i)
                add_count(pair_n, &current_seq, i, counted, orig_alphabet_size);
        }

        size_t new_alphabet_size = orig_alphabet_size;
        while(true) {
            // find pair with maximal number of non-overlapping occurrences
            s_pair max_pair;
            size_t max_count = 0;
            for(size_t i = 0; i < new_alphabet_size*new_alphabet_size; ++i) {
                if(pair_n[i] > max_count) {
                    max_count = pair_n[i];
                    unsigned left = unsigned(i / new_alphabet_size);
                    unsigned right = unsigned(i - left * new_alphabet_size);
                    max_pair = s_pair(left, right);
                }
            }

            iteration_timer.stop();
            iteration_time = iteration_timer.timeElapsed();

            // stopping criteria
            if(nStatesSave_init_empty && prev_max_count == max_count) {
                nStates2seqs[2] = *current_seqs_p; // copying sequences
                nStates2alphabet_size[2] = new_alphabet_size;
                break;
            } else if(!nStatesSave_init_empty) {
                double time_saved_per_eval = double(max_count) * double(nStates2t_MV[no_states_save]) - double(nStates2t_MM[no_states_save]);
                double time_saved_on_evals = double(min_no_eval) * time_saved_per_eval;
                double total_time_saved = time_saved_on_evals - iteration_time;

                if(total_time_saved <= 0) {
                    nStates2seqs[no_states_save] = *current_seqs_p; // copying sequences
                    nStates2alphabet_size[no_states_save] = new_alphabet_size;

                    nStatesSave.pop_back();
                    if(nStatesSave.empty())
                        break;
                    else
                        no_states_save = nStatesSave.back();
                }
            }

            iteration_timer.start();

            // save the components of max_pair in symbol2pair
            symbol2pair[unsigned(new_alphabet_size)] = max_pair;

            pair_n.clear();
            pair_n.resize( (new_alphabet_size+1) * (new_alphabet_size+1));
            next_seqs_p = new std::vector<std::vector<unsigned> >();

            for(std::vector<std::vector<unsigned> >::iterator it = current_seqs_p->begin(); it != current_seqs_p->end(); ++it) {
                std::vector<unsigned> &current_seq = (*it);
                next_seqs_p->push_back(std::vector<unsigned>());
                std::vector<unsigned> &next_seq = next_seqs_p->back();

                next_seq.push_back( current_seq[0] );
                // replace every occurrence of max_pair with a new alphabet symbol and count pairs at the same time
                size_t i = 1; // first position is special because we have not seen any pairs yet.
                if(matches(max_pair, &current_seq, i)) {
                    next_seq.push_back(unsigned(new_alphabet_size));
                    ++i;
                } else {
                    next_seq.push_back(current_seq[i]);
                }
                // do the rest of the sequence
                for(i = i+1; i < current_seq.size() - 1; i++) {
                    if(matches(max_pair, &current_seq, i)) {
                        next_seq.push_back(unsigned(new_alphabet_size));
                        ++i;
                    } else {
                        next_seq.push_back(current_seq[i]);
                    }

                    add_count(pair_n, &next_seq, next_seq.size() - 2, counted, new_alphabet_size + 1);
                }
                if(i == current_seq.size() - 1) { // push last symbol of sequence
                    next_seq.push_back(current_seq[i]);
                    add_count(pair_n, &next_seq, next_seq.size() - 2, counted, new_alphabet_size + 1);
                }
            }

            // set up for next round
            prev_max_count = max_count;
            delete prev_seqs_p;
            prev_seqs_p = current_seqs_p;
            current_seqs_p = next_seqs_p;
            new_alphabet_size++;
        }

        delete prev_seqs_p;
        delete current_seqs_p;
        // delete next_seqs_p; has already been freed!
    }

    void Viterbi::read_seq(const std::string &seq_filename,
                           const size_t alphabet_size,
                           const size_t min_no_eval) {
        std::vector<size_t> nStatesSave;
        read_seq(seq_filename, alphabet_size, nStatesSave, min_no_eval);
    }

    void Viterbi::read_seq(const std::string &seq_filename,
                           const size_t alphabet_size,
                           std::vector<size_t> &nStatesSave,
                           const size_t min_no_eval) {
        std::vector<std::vector<unsigned> > sequences;
        sequences.push_back(std::vector<unsigned>());
        readSeq(sequences[0], seq_filename);
        init_data_structure(sequences, alphabet_size, nStatesSave, min_no_eval);
    }


    void Viterbi::read_seq(const std::string &seq_filename,
                           const size_t alphabet_size,
                           const size_t no_states,
                           const size_t min_no_eval) {
        std::vector<size_t> nStatesSave;
        nStatesSave.push_back(no_states);
        read_seq(seq_filename, alphabet_size, nStatesSave, min_no_eval);
    }

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
           B.get_height() != A.get_width() || B.get_width() != orig_alphabet_size) {
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
        const std::vector<std::vector< unsigned> > *sequences = 0;
        for(std::map<size_t, std::vector<std::vector<unsigned> > >::const_iterator it = nStates2seqs.begin(); it != nStates2seqs.end(); ++it) {
            if(it->first >= no_states) {
                sequences = &(it->second);
                alphabet_size = nStates2alphabet_size.find(it->first)->second;
                break;
            }
        }
        if(sequences == 0) {
            sequences = &(nStates2seqs.rbegin()->second);
            alphabet_size = nStates2alphabet_size.rbegin()->second;
        }

        Matrix *symbol2matrix = new Matrix[alphabet_size];

        compute_symbol2matrix(symbol2matrix, A, B, alphabet_size);

        double ll = 0.0;
        for(std::vector<std::vector<unsigned> >::const_iterator it = sequences->begin(); it != sequences->end(); ++it) {
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
        for(size_t i = orig_alphabet_size; i < alphabet_size; ++i) {
            const s_pair symbol_pair = symbol2pair.find(unsigned(i))->second;
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
