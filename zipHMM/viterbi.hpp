#ifndef VITERBI_HPP
#define VITERBI_HPP

#include "matrix.hpp"

#include <vector>
#include <map>

namespace zipHMM {

    class Viterbi {

    protected:
        size_t orig_seq_length;
        size_t orig_alphabet_size;
        std::map<unsigned, s_pair> symbol2pair;
        std::map<size_t, size_t> nStates2alphabet_size;
        std::map<size_t, std::vector<std::vector<unsigned> > > nStates2seqs;

    public:
        Viterbi() { }
        ~Viterbi() { }

        void read_seq(const std::string &seq_filename,
                      const size_t alphabet_size,
                      const size_t min_no_eval = 1);

        void read_seq(const std::string &seq_filename,
                      const size_t alphabet_size,
                      std::vector<size_t> &nStatesSave,
                      const size_t min_no_eval = 1);

        void read_seq(const std::string &seq_filename,
                      const size_t alphabet_size,
                      const size_t no_states,
                      const size_t min_no_eval); // convenient

        double forward(const Matrix &pi, const Matrix &A, const Matrix &B) const;


    double viterbi_orig(const std::vector<unsigned> &seq,
                        const Matrix &initProbs,
                        const Matrix &transProbs,
                        const Matrix &emProbs,
                        std::vector<unsigned> &viterbi_path);

    double viterbi_orig(const std::string &seq_filename,
                        const Matrix &pi,
                        const Matrix &A,
                        const Matrix &B,
                        std::vector<unsigned> &viterbi_path);

    double viterbi_comp(const std::vector<unsigned> &seq,
                        const Matrix &initProbs,
                        const Matrix &transProbs,
                        const Matrix &emProbs,
                        std::vector<unsigned> &viterbi_path);

    double viterbi_comp(const std::string &seq_filename,
                        const Matrix &pi,
                        const Matrix &A,
                        const Matrix &B,
                        std::vector<unsigned> &viterbi_path);

    };
}

#endif
