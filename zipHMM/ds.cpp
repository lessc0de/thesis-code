#include "ds.hpp"

#include <utility>
#include <algorithm>
#include <map>
#include <vector>
#include <deque>
#include <stack>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <dirent.h>

#include "matrix.hpp"
#include "seq_io.hpp"
#include "io_utils.hpp"
#include "hmm_utils.hpp"
#include "PThreadProcessingDevice.hpp"
#include "performance_description.hpp"
#include "Stage2JobControl.hpp"
#include "debug.hpp"

namespace zipHMM {

    size_t DS::get_seq_length(size_t no_states) const {
        size_t length = 0;

        for(std::map<size_t, std::vector<std::vector<unsigned> > >::const_iterator it = nStates2seqs.begin(); it != nStates2seqs.end(); it++) {
            if(it->first >= no_states) {
                const std::vector<std::vector<unsigned> > &sequences = it->second;

                for(std::vector<std::vector<unsigned> >::const_iterator it2 = sequences.begin(); it2 != sequences.end(); ++it2)
                    length += it2->size();

                return length;
            }
        }

        const std::vector<std::vector<unsigned> > &sequences = nStates2seqs.rbegin()->second;

        for(std::vector<std::vector<unsigned> >::const_iterator it = sequences.begin(); it != sequences.end(); ++it)
            length += it->size();

        return length;
    }

    size_t DS::get_alphabet_size(size_t no_states) const {
        for(std::map<size_t, std::vector<std::vector<unsigned> > >::const_iterator it = nStates2seqs.begin(); it != nStates2seqs.end(); ++it) {
            if(it->first >= no_states)
                return nStates2alphabet_size.find(it->first)->second;
        }
        return nStates2alphabet_size.rbegin()->second;
    }

    s_pair DS::get_pair(unsigned symbol) const {
        if(symbol2pair.count(symbol) != 0)
            return symbol2pair.find(symbol)->second;
        return std::pair<unsigned, unsigned>((unsigned)-1, (unsigned)-1);
    }

    void DS::write_to_directory(const std::string &directory) const {
        std::string wd = get_working_directory();
        std::string absolute_dir_name;
        if(!directory[0] != '/')
            absolute_dir_name = wd + "/" + directory;
        else
            absolute_dir_name = directory;
        std::string data_structure_filename = absolute_dir_name + "/data_structure";
        std::string nStates2seq_absolute_dir_name = absolute_dir_name + "/nStates2seq";

        // create directory
        mk_dir(absolute_dir_name);
        mk_dir(nStates2seq_absolute_dir_name);

        // write basic data structure
        write_data_structure(data_structure_filename);

        // write nStates2seq
        write_seqs(nStates2seq_absolute_dir_name);
    }

    void DS::init_data_structure(const std::vector<std::vector<unsigned> > &sequences, const size_t alphabet_size, std::vector<size_t> &nStatesSave, const size_t min_no_eval) {
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

        // If min_no_eval is set to 0, don't do any compression at all.
        if (min_no_eval == 0) {
            nStates2seqs[2] = *current_seqs_p; // copying sequences
            nStates2alphabet_size[2] = orig_alphabet_size;
            delete prev_seqs_p;
            delete current_seqs_p;
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

    void DS::read_seq(const std::string &seq_filename, const size_t alphabet_size, std::vector<size_t> &nStatesSave, const size_t min_no_eval) {
        std::vector<std::vector<unsigned> > sequences;
        sequences.push_back(std::vector<unsigned>());
        readSeq(sequences[0], seq_filename);
        init_data_structure(sequences, alphabet_size, nStatesSave, min_no_eval);
    }

    void DS::read_seq(const std::string &seq_filename, const size_t alphabet_size, const size_t no_states, const size_t min_no_eval) {
        std::vector<size_t> nStatesSave;
        nStatesSave.push_back(no_states);
        read_seq(seq_filename, alphabet_size, nStatesSave, min_no_eval);
    }

    void DS::read_seq(const std::string &seq_filename, const size_t alphabet_size, const size_t min_no_eval) {
        std::vector<size_t> nStatesSave;
        read_seq(seq_filename, alphabet_size, nStatesSave, min_no_eval);
    }

    void DS::read_seq_directory(const std::string &dirname, const size_t alphabet_size, std::vector<size_t> &nStatesSave, const size_t min_no_eval) {
        std::vector<std::vector<unsigned> > sequences;
        DIR *dir;
        struct dirent *ent;
        std::string suffix = ".seq";

        if ((dir = opendir (dirname.c_str())) != NULL) {
            while ((ent = readdir (dir)) != NULL) {
                std::string filename = ent->d_name;
                if(filename.size() >= suffix.size() && filename.compare(filename.size() - suffix.size(), suffix.size(), suffix) == 0) {
                    std::string path = dirname + "/" + filename;
                    std::vector<unsigned> sequence;
                    zipHMM::readSeq(sequence, path);
                    sequences.push_back(sequence);
                }
            }
            closedir (dir);
        } else {
            /* could not open directory */
            perror ("could not open directory");
        }

        init_data_structure(sequences, alphabet_size, nStatesSave, min_no_eval);
    }

    void DS::read_seq_directory(const std::string &dirname, const size_t alphabet_size, const size_t no_states, const size_t min_no_eval) {
        std::vector<size_t> nStatesSave;
        nStatesSave.push_back(no_states);
        read_seq_directory(dirname, alphabet_size, nStatesSave, min_no_eval);
    }

    void DS::read_seq_directory(const std::string &dirname, const size_t alphabet_size, const size_t min_no_eval) {
        std::vector<size_t> nStatesSave;
        read_seq_directory(dirname, alphabet_size, nStatesSave, min_no_eval);
    }

    void DS::read_from_directory(const std::string &directory) {
        read_data_structure_from_directory(directory);
        read_all_seqs(directory + "/nStates2seq");
    }

    void DS::read_from_directory(const std::string &directory, const size_t no_states) {
        read_data_structure_from_directory(directory);
        std::map<size_t, size_t> new_nStates2alphabet_size;
        // remove what we read too much
        new_nStates2alphabet_size[no_states] = nStates2alphabet_size[no_states];
        nStates2alphabet_size = new_nStates2alphabet_size;

        std::stringstream ss;
        ss << directory << "/nStates2seq" << "/" << no_states;
        read_no_states_seqs_from_directory(ss.str(), no_states);
    }

    void DS::write_seqs(const std::string &nstates2seq_absolute_dir_name) const {
        for(std::map<size_t, std::vector<std::vector<unsigned> > >::const_iterator it = nStates2seqs.begin(); it != nStates2seqs.end(); ++it) {
            size_t nStates = (*it).first;
            const std::vector<std::vector<unsigned> > &sequences = (*it).second;

            std::stringstream nStatesDirname_stream;
            nStatesDirname_stream << nstates2seq_absolute_dir_name << "/" << nStates;
            std::string nStatesDirname = nStatesDirname_stream.str();
            mk_dir(nStatesDirname);

            for(size_t i = 0; i < sequences.size(); ++i) {
                std::stringstream seq_filename_stream;
                seq_filename_stream << nStatesDirname << "/" << i << ".seq";

                write_seq(seq_filename_stream.str(), nStates, i);
            }
        }
    }

    void DS::write_seq(const std::string &seq_filename, size_t no_states, size_t seq_no) const {
        std::ofstream seq_stream(seq_filename.c_str());
        if(!seq_stream) {
            std::cerr << "Unable to open \"" << seq_filename << "\"" << std::endl;
            exit(-1);
        }

        write_seq(seq_stream, no_states, seq_no);
    }

    void DS::write_seq(std::ofstream &out, size_t no_states, size_t seq_no) const {
        const std::vector<unsigned> &sequence = (nStates2seqs.find(no_states)->second)[seq_no];
        for(std::vector<unsigned>::const_iterator it = sequence.begin(); it != sequence.end(); ++it)
            out << (*it) << " ";
    }

    void DS::write_data_structure(const std::string &data_structure_filename) const {
        std::ofstream data_structure_out(data_structure_filename.c_str());
        if(!data_structure_out) {
            std::cerr << "Unable to open \"" << data_structure_filename << "\"" << std::endl;
            exit(-1);
        }

        write_data_structure(data_structure_out);
    }

    void DS::write_data_structure(std::ofstream &out) const {
        out << "orig_alphabet_size" << std::endl;
        out << orig_alphabet_size << std::endl;

        out << "orig_seq_length" << std::endl;
        out << orig_seq_length << std::endl;

        out << "nStates2alphabet_size" << std::endl;
        for(std::map<size_t, size_t>::const_iterator it = nStates2alphabet_size.begin(); it != nStates2alphabet_size.end(); ++it)
            out << (*it).first << " " << (*it).second <<std::endl;

        out << "symbol2pair" << std::endl;
        for(std::map<unsigned, s_pair>::const_iterator it = symbol2pair.begin(); it != symbol2pair.end(); ++it)
            out << (*it).first << " " << (*it).second.first << " " << (*it).second.second <<std::endl;
    }

    void DS::read_data_structure_from_directory(const std::string &directory) {
        std::string full_filename = directory + "/data_structure";

        std::ifstream data_structure_in(full_filename.c_str());
        if(!data_structure_in) {
            std::cerr << "Unable to open \"" << full_filename << "\"" << std::endl;
            exit(-1);
        }

        read_data_structure_from_stream(data_structure_in);
        data_structure_in.close();
    }

    void DS::read_data_structure_from_stream(std::ifstream &in) {
        read_token_or_exit(in, "orig_alphabet_size");
        orig_alphabet_size = read_or_exit<size_t>(in, "original alphabet size");

        read_token_or_exit(in, "orig_seq_length");
        orig_seq_length = read_or_exit<size_t>(in, "original alphabet size");

        read_token_or_exit(in, "nStates2alphabet_size");
        while(!read_token_or_tell(in, "symbol2pair")) {
            size_t no_states = read_or_exit<size_t>(in, "no_states");
            nStates2alphabet_size[no_states] = read_or_exit<size_t>(in, "alphabet_size");
        }

        for(size_t i = orig_alphabet_size; i < nStates2alphabet_size.begin()->second; ++i) {
            unsigned pair_symbol  = read_or_exit<unsigned>(in, "pair symbol");
            unsigned left_symbol  = read_or_exit<unsigned>(in, "left symbol");
            unsigned right_symbol = read_or_exit<unsigned>(in, "right symbol");
            symbol2pair[pair_symbol] = s_pair(left_symbol, right_symbol);
        }
    }

    void DS::read_all_seqs(const std::string &directory) {
        std::vector<std::string> dirnames;

        DIR *dpdf;
        struct dirent *epdf;

        dpdf = opendir(directory.c_str());
        if (dpdf != NULL) {
            while ((epdf = readdir(dpdf))) {
                if(!(std::strcmp(epdf->d_name, ".") == 0) && !(std::strcmp(epdf->d_name, "..") == 0) && !(std::strcmp(epdf->d_name, ".svn") == 0)) {
                    dirnames.push_back(epdf->d_name);
                }
            }
        }

        for(std::vector<std::string>::iterator it = dirnames.begin(); it != dirnames.end(); ++it) {
            std::string dirname = (*it);
            size_t no_states = (size_t) std::atol(dirname.substr(0, dirname.length()-4).c_str());

            std::stringstream no_states_dir_stream;
            no_states_dir_stream << directory << "/" << no_states;
            std::string no_states_dir = no_states_dir_stream.str();

            read_no_states_seqs_from_directory(no_states_dir, no_states);
        }
    }

    void DS::read_no_states_seqs_from_directory(const std::string &no_states_dir, size_t no_states) {
        std::vector<std::string> filenames;

        DIR *dpdf;
        struct dirent *epdf;

        dpdf = opendir(no_states_dir.c_str());
        if (dpdf != NULL) {
            while ((epdf = readdir(dpdf))) {
                if(!(std::strcmp(epdf->d_name, ".") == 0) && !(std::strcmp(epdf->d_name, "..") == 0) && !(std::strcmp(epdf->d_name, ".svn") == 0)) {
                    filenames.push_back(epdf->d_name);
                }
            }
        }

        for(std::vector<std::string>::const_iterator it = filenames.begin(); it != filenames.end(); ++it) {
            std::string path = no_states_dir + "/" + (*it);

            std::ifstream in(path.c_str());

            if(!in) {
                std::cerr << "Unable to open \"" << path << "\"" << std::endl;
                exit(-1);
            }
            read_seq_from_stream(in, no_states);

            in.close();
        }
    }


    void DS::read_seq_from_stream(std::ifstream &in, size_t no_states) {
        std::vector<unsigned> res;
        unsigned tmp = 0;
        while(in >> tmp)
            res.push_back(tmp);
        nStates2seqs[no_states].push_back(res);
    }
}
