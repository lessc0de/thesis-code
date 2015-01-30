#ifndef PTHREADPROCESSINGDEVICE_HPP
#define PTHREADPROCESSINGDEVICE_HPP

#include "ProcessingDevice.hpp"
#include "prob_spaces.hpp"
#include "hmm_utils.hpp"
#include "Stage2JobControl.hpp"
#include "MapReduceJobControl.hpp"

#include <pthread.h>
#include <map>
#include <stack>
#include <vector>

namespace zipHMM {

  class PThreadProcessingDevice : public ProcessingDevice {

    const Matrix *pi;
    const Matrix *A;
    const Matrix *B;
    double *symbol2scale;
    Matrix *symbol2matrix;
    const Matrix *symbol2logMatrix;
    const std::map<unsigned, s_pair> *symbol2pair;

    pthread_t thread;

    const std::vector<unsigned> *seq;
    const std::vector<std::vector<unsigned> > *seqs;

  public:
    
    unsigned id;

    PThreadProcessingDevice(unsigned iid)
      : pi(0), A(0), B(0),
	symbol2scale(0), symbol2matrix(0), symbol2logMatrix(0), symbol2pair(0), 
	thread(), 
	seq(0) {
      
      id = iid;

    }
    
    PThreadProcessingDevice(const PThreadProcessingDevice &);
    PThreadProcessingDevice &operator=(const PThreadProcessingDevice &);

    ~PThreadProcessingDevice() { clear(); }

    void likelihoodVector(Stage2JobControl &control);
    void likelihoodMatrix(Stage2JobControl &control);

    void mapReduceLoglikelihood(MapReduceJobControl &control);

    void join();

    void setParameters(const Matrix *pi, const Matrix *A, const Matrix *B, 
		       const std::map<unsigned, s_pair> *symbol2pair, double *symbol2scale, Matrix *symbol2matrix);

    void setHMM(const Matrix *pi, const Matrix *A, const Matrix *B);
    void setSymbol2scale(double *symbol2scale) { this->symbol2scale = symbol2scale; }
    void setSymbol2matrix(Matrix *symbol2matrix) { this->symbol2matrix = symbol2matrix; }
    void setSymbol2pair(const std::map<unsigned, s_pair> *symbol2pair) { this->symbol2pair = symbol2pair; }
    void setSeq(const std::vector<unsigned> *seq) { this->seq = seq; }
    void setSeqs(const std::vector<std::vector<unsigned> > *seqs) { this->seqs = seqs; };

    const Matrix &getInitProbs() const { return *pi; }
    const Matrix &getTransProbs() const { return *A; }
    const Matrix &getEmProbs() const { return *B; }

    Matrix &getMatrix(unsigned symbol) { return symbol2matrix[symbol]; }
    const Matrix &getLogMatrix(unsigned symbol) const { return symbol2logMatrix[symbol]; }
    double &getScale(unsigned symbol) { return symbol2scale[symbol]; }
    const s_pair &getPair(unsigned symbol) const { return symbol2pair->find(symbol)->second; }

    const std::vector<unsigned> &getSeq() const { return *seq; }
    const std::vector<unsigned> &getSeq(size_t seqNo) const {return seqs->at(seqNo); }

    void computeSymbol2logMatrix();

    void clearCache();
    void clearResult();
    void clear();
  };
} // end namespace

#endif
