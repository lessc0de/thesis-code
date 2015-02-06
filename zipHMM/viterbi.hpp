#ifndef VITERBI_HPP
#define VITERBI_HPP

#include "matrix.hpp"
#include "hmm_utils.hpp"
#include "ds.hpp"

#include <vector>
#include <map>

namespace zipHMM {

    class Viterbi {

    public:
        Viterbi() { }
        ~Viterbi() { }

        double viterbi_seq(const Matrix &pi, const Matrix &A, const Matrix &B,
                           const std::vector<unsigned> &sequence,
                           const Matrix *symbol2matrix,
                           const Matrix *symbol2matrixR,
                           const bool compute_path,
                           std::vector<unsigned> &viterbi_path) const;

        double viterbi(const Matrix &pi, const Matrix &A, const Matrix &B) const;

        double viterbi(const Matrix &pi, const Matrix &A, const Matrix &B,
                       std::vector<unsigned> &viterbi_path) const;

        // Legacy methods.
        size_t get_orig_seq_length() const {
            return ds.get_orig_seq_length();
        }

        size_t get_orig_alphabet_size() const {
            return ds.get_orig_alphabet_size();
        }

        size_t get_seq_length(size_t no_states) const {
            return ds.get_seq_length(no_states);
        };

        size_t get_alphabet_size(size_t no_states) const {
            return ds.get_alphabet_size(no_states);
        };

        s_pair get_pair(unsigned symbol) const {
            return ds.get_pair(symbol);
        };

        std::map<unsigned, s_pair> get_symbol2pair() const {
            return ds.get_symbol2pair();
        };

        std::map<size_t, size_t> get_nStates2alphabet_size() const {
            return ds.get_nStates2alphabet_size();
        };

        std::map<size_t, std::vector<std::vector<unsigned> > > get_nStates2seqs() const {
            return ds.get_nStates2seqs();
        };

        void write_to_directory(const std::string &directory) const {
            ds.write_to_directory(directory);
        };

        void init_data_structure(const std::vector<std::vector<unsigned> > &sequences,
                                 const size_t alphabet_size, std::vector<size_t> &nStatesSave,
                                 const size_t min_no_eval) {
            ds.init_data_structure(sequences, alphabet_size, nStatesSave, min_no_eval);
        };

        void read_seq(const std::string &seq_filename, const size_t alphabet_size,
                      const size_t min_no_eval = 1) {
            ds.read_seq(seq_filename, alphabet_size, min_no_eval);
        };

        void read_seq(const std::string &seq_filename, const size_t alphabet_size,
                      std::vector<size_t> &nStatesSave, const size_t min_no_eval = 1) {
            ds.read_seq(seq_filename, alphabet_size, nStatesSave, min_no_eval);
        };

        void read_seq(const std::string &seq_filename, const size_t alphabet_size,
                      const size_t no_states, const size_t min_no_eval) {
            ds.read_seq(seq_filename, alphabet_size, no_states, min_no_eval);
        }; // convenient

        void read_seq_directory(const std::string &dir_filename, const size_t alphabet_size,
                                const size_t min_no_eval = 1) {
            ds.read_seq_directory(dir_filename, alphabet_size, min_no_eval);
        };

        void read_seq_directory(const std::string &dir_filename, const size_t alphabet_size,
                                std::vector<size_t> &nStatesSave, const size_t min_no_eval = 1) {
            ds.read_seq_directory(dir_filename, alphabet_size, nStatesSave, min_no_eval);
        };

        void read_seq_directory(const std::string &dir_filename, const size_t alphabet_size,
                                const size_t no_states, const size_t min_no_eval) {
            ds.read_seq_directory(dir_filename, alphabet_size, no_states, min_no_eval);
        };

        void read_from_directory(const std::string &dirname) {
            ds.read_from_directory(dirname);
        };

        void read_from_directory(const std::string &directory, const size_t no_states) {
            ds.read_from_directory(directory, no_states);
        };

    private:
        DS ds;
        void compute_symbol2matrix(Matrix *symbol2matrix,
                                   Matrix *symbol2matrixR,
                                   const Matrix &A,
                                   const Matrix &B,
                                   size_t alphabet_size) const;

        void deduct_path(const std::vector<unsigned> &sequence,
                         std::vector<unsigned> &path,
                         const Matrix *symbol2matrixR) const;

        double viterbi_helper(const Matrix &pi, const Matrix &A, const Matrix &B,
                              const bool compute_path, std::vector<unsigned> &viterbi_path) const;

    };

    double viterbi(const std::vector<unsigned> &seq,
                   const Matrix &initProbs,
                   const Matrix &transProbs,
                   const Matrix &emProbs,
                   std::vector<unsigned> &viterbi_path);

    double viterbi(const std::string &seq_filename,
                   const Matrix &pi,
                   const Matrix &A,
                   const Matrix &B,
                   std::vector<unsigned> &viterbi_path);

}

#endif
