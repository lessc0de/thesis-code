#ifndef FORWARDER_HPP
#define FORWARDER_HPP

#include "matrix.hpp"
#include "PThreadProcessingDevice.hpp"
#include "performance_description.hpp"
#include "hmm_utils.hpp"
#include "ds.hpp"

#include <vector>
#include <string>

namespace zipHMM {

    class SimpleStopForwarder;

    class Forwarder {

    public:
        Forwarder() { }
        ~Forwarder() { }

        double forward_seqs(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales) const;
        double forward_seq(const Matrix &pi, const Matrix &A, const Matrix &B,
                           const std::vector<unsigned> &sequence, const double *symbol2scale,
                           const Matrix *symbol2matrix, std::vector<double> &scales) const;

        double forward(const Matrix &pi, const Matrix &A, const Matrix &B) const;

        double forward(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales) const;

        double pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B,
                               const std::string &device_filename = DEFAULT_DEVICE_FILENAME) const;

        double pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B,
                               const DeviceDescriptor &device_descriptor) const;

        double mr_pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B,
                                  const DeviceDescriptor &device_descriptor) const;

        double mr_pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B,
                                  const std::string &device_filename = DEFAULT_DEVICE_FILENAME) const;

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

        double forward_helper(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales) const;

        void compute_symbol2scale_and_symbol2matrix(Matrix *symbol2matrix, double *symbol2scale,
                                                    const Matrix &A, const Matrix &B,
                                                    size_t alphabet_size) const;

        void pthread_compute_symbol2scale_and_symbol2matrix(double *symbol2scale,
                                                            Matrix *symbol2matrix,
                                                            const Matrix &A, const Matrix &B,
                                                            const size_t alphabet_size,
                                                            const std::vector<ProcessingDevice*> &devices) const;

    };

}

#endif
