#include "Multilign_object.h"
#include <svm.h>
//#include "../Multifind/yinghan_math.h"

//! Multifind_object Class.
/*!
    The Multifind_object class provides an entry point for the Multifind algorithm.
*/


class Multifind_object:public Multilign_object{
  
 public:
  //!Constructor:
  //!\param outputmultifind is the name of the Multifind output file to which the output is written to.
  //!\param ctfiles is a vector of strings storing the names of the ct files to which the output structures are written to.
  //!\param inputalignment is a vector of strings storing the input sequences in the alignment (with gaps).
  //!\param inputsequences is a vector of strings storing the input sequences in the alignment (without gaps).
  //!\param processors is a interger indicating the number of processors required by Multifind in smp calculations.(only applicable in smp version)
 //!\param progress is a TProgressDialog for reporting progress of the calculation to the user.  The default value of NULL means that no communication is provided.
  Multifind_object(const string &outputmultifind, const vector<string> &ctfiles,  const vector<string> &inputalignment, const vector<string> &inputsequences, const int &processors,ProgressHandler *progress=NULL);
  //! The core function doing Multilign calculation and SVM prediction.
  int Multifind_Predict();
  

 private:
  int num_processors;
  double sum_multifind(vector<double>& series);
  double average_multifind(vector<double>& series);
  double variation_multifind(vector<double>& series);
  double sum_multifind(vector<float>& series);
  double average_multifind(vector<float>& series);
  double variation_multifind(vector<float>& series);
  double sum_multifind(vector<int>& series);
  double average_multifind(vector<int>& series);
  double variation_multifind(vector<int>& series);
  double normalized_ensemble_defect(RNA* rna);
  double entropy(vector<char> column);
  double average_entropy();
  vector<double> single_z_predict(string sequence,svm_model* model_folding_average,svm_model* model_folding_std,svm_model* model_ensemble_average,svm_model* model_ensemble_std);
  int get_gap(string sequence);
  //  string compact(string& sequence);
  double common_energies();
  vector<double> predict_ncRNA_probabilities(double sci,double entropy,double single_z,double ensemble_defect_z);
};
