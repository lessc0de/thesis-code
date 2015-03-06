#ifndef HMM_SUITE_HPP
#define HMM_SUITE_HPP

#include "matrix.hpp"
#include "PThreadProcessingDevice.hpp"
#include "performance_description.hpp"
#include "hmm_utils.hpp"

#include <vector>
#include <string>

namespace zipHMM {

    class SimpleStopForwarder;

    class HMMSuite {

    protected:
        size_t orig_seq_length;
        size_t orig_alphabet_size;
        std::map<unsigned, s_pair> symbol2pair;
        std::map<size_t, size_t> nStates2alphabet_size;
        std::map<size_t, std::vector<std::vector<unsigned> > > nStates2seqs;

    public:
        HMMSuite() { }
        ~HMMSuite() { }

        double forward_seqs(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales) const;
        double forward_seq(const Matrix &pi, const Matrix &A, const Matrix &B,
                           const std::vector<unsigned> &sequence, const double *symbol2scale,
                           const Matrix *symbol2matrix, std::vector<double> &scales, bool compute_matrix, Matrix &forward_table) const;

        double forward(const Matrix &pi, const Matrix &A, const Matrix &B) const;

        double forward(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales) const;

        double forward(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales, Matrix &forward_table) const;

        double pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B,
                               const std::string &device_filename = DEFAULT_DEVICE_FILENAME) const;

        double pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B,
                               const DeviceDescriptor &device_descriptor) const;

        double mr_pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B,
                                  const DeviceDescriptor &device_descriptor) const;

        double mr_pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B,
                                  const std::string &device_filename = DEFAULT_DEVICE_FILENAME) const;

        void backward_seq(const Matrix &pi, const Matrix &A, const Matrix &B,
                          const std::vector<unsigned> &sequence, const std::vector<double> &scales,
                          const Matrix *symbol2matrix, Matrix &backward_table) const;

        void backward(const Matrix &pi, const Matrix &A, const Matrix &B, const std::vector<double> &scales, Matrix &backward_table) const;

        double viterbi_seq(const Matrix &pi, const Matrix &A, const Matrix &B,
                           const std::vector<unsigned> &sequence,
                           const Matrix *symbol2matrix,
                           const Matrix *symbol2matrixR,
                           const bool compute_path,
                           const bool memory_save,
                           std::vector<unsigned> &viterbi_path) const;

        double viterbi(const Matrix &pi, const Matrix &A, const Matrix &B) const;

        double viterbi(const Matrix &pi, const Matrix &A, const Matrix &B,
                       std::vector<unsigned> &viterbi_path) const;

        double viterbi(const Matrix &pi, const Matrix &A, const Matrix &B,
                       const bool memory_save, std::vector<unsigned> &viterbi_path) const;

        void posterior_decoding(const Matrix &pi, const Matrix &A, const Matrix &B,
                                std::vector<unsigned> &posterior_path) const;

        size_t get_orig_seq_length() const { return orig_seq_length; }
        size_t get_orig_alphabet_size() const { return orig_alphabet_size; }
        size_t get_seq_length(size_t no_states) const;
        size_t get_alphabet_size(size_t no_states) const;
        s_pair get_pair(unsigned symbol) const;

        std::map<unsigned, s_pair> get_symbol2pair() const { return symbol2pair; };
        std::map<size_t, size_t> get_nStates2alphabet_size() const { return nStates2alphabet_size; };
        std::map<size_t, std::vector<std::vector<unsigned> > > get_nStates2seqs() const { return nStates2seqs; };

        void write_to_directory(const std::string &directory) const;

        void init_data_structure(const std::vector<std::vector<unsigned> > &sequences, const size_t alphabet_size, std::vector<size_t> &nStatesSave, const size_t min_no_eval);

        void read_seq(const std::string &seq_filename, const size_t alphabet_size, const size_t min_no_eval = 1);
        void read_seq(const std::string &seq_filename, const size_t alphabet_size, std::vector<size_t> &nStatesSave, const size_t min_no_eval = 1);
        void read_seq(const std::string &seq_filename, const size_t alphabet_size, const size_t no_states, const size_t min_no_eval); // convenient

        void read_seq_directory(const std::string &dir_filename, const size_t alphabet_size, const size_t min_no_eval = 1);
        void read_seq_directory(const std::string &dir_filename, const size_t alphabet_size, std::vector<size_t> &nStatesSave, const size_t min_no_eval = 1);
        void read_seq_directory(const std::string &dir_filename, const size_t alphabet_size, const size_t no_states, const size_t min_no_eval);

        void read_from_directory(const std::string &dirname);
        void read_from_directory(const std::string &directory, const size_t no_states);

    private:
        void write_seqs(const std::string &nstates2seq_absolute_dir_name) const;
        void write_seq(const std::string &seq_filename, size_t no_states, size_t seq_no) const;
        void write_seq(std::ofstream &seq_stream, size_t no_states, size_t seq_no) const;

        void write_data_structure(const std::string &data_structure_filename) const;
        void write_data_structure(std::ofstream &data_structure_stream) const;

        void read_data_structure_from_directory(const std::string &directory);
        void read_data_structure_from_stream(std::ifstream &data_structure_stream);

        void read_all_seqs(const std::string &directory);
        void read_no_states_seqs_from_directory(const std::string &no_states_dir, size_t no_states);
        void read_seq_from_stream(std::ifstream &in, size_t no_states);

        double forward_helper(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales, bool compute_path, Matrix &forward_table) const;

        void forward_compute_symbol2scale_and_symbol2matrix(Matrix *symbol2matrix, double *symbol2scale,
                                                            const Matrix &A, const Matrix &B,
                                                            size_t alphabet_size) const;

        void pthread_compute_symbol2scale_and_symbol2matrix(double *symbol2scale,
                                                            Matrix *symbol2matrix,
                                                            const Matrix &A, const Matrix &B,
                                                            const size_t alphabet_size,
                                                            const std::vector<ProcessingDevice*> &devices) const;

        void backward_compute_symbol2scale_and_symbol2matrix(Matrix *symbol2matrix, double *symbol2scale,
                                                             const Matrix &A, const Matrix &B,
                                                             size_t alphabet_size) const;

        void viterbi_compute_symbol2matrix(Matrix *symbol2matrix,
                                           Matrix *symbol2matrixR,
                                           const Matrix &A,
                                           const Matrix &B,
                                           size_t alphabet_size) const;

        void deduct_path(const std::vector<unsigned> &sequence,
                         std::vector<unsigned> &path,
                         const Matrix *symbol2matrixR) const;

        double viterbi_helper(const Matrix &pi, const Matrix &A, const Matrix &B,
                              const bool compute_path, const bool memory_save, std::vector<unsigned> &viterbi_path) const;
    };

}

#endif
