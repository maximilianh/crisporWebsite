#include <fstream>
#include <cstdlib>
#include <svm.h>
#include "Multilign_object.h"
#include "Multifind_object.h"


double Multifind_object::sum_multifind(vector<double>& series){
  double sum=0;
  for(int i=0;i<series.size();++i){
    sum+=series[i];
  }
  return sum;
}

double Multifind_object::average_multifind(vector<double>& series){
   return (sum_multifind(series))/(series.size());
}

double Multifind_object::variation_multifind(vector<double>& series){
  double average=average_multifind(series);
  double sum=0;
  for(int i=0;i<series.size();++i){
    sum+=(series[i]-average)*(series[i]-average);
  }
  return sqrt(sum/(series.size()-1));
}

double Multifind_object::sum_multifind(vector<float>& series){
  double sum=0;
  for(int i=0;i<series.size();++i){
    sum+=series[i];
  }
  return sum;
}

double Multifind_object::average_multifind(vector<float>& series){
  return (sum_multifind(series))/series.size();
}


double Multifind_object::variation_multifind(vector<float>& series){
  double average=average_multifind(series);
  double sum=0;
  for(int i=0;i<series.size();++i){
    sum+=(series[i]-average)*(series[i]-average);
  }
  return sqrt(sum/(series.size()-1));
}


double Multifind_object::sum_multifind(vector<int>& series){
  double sum=0;
  for(int i=0;i<series.size();++i){
    sum+=series[i];
  }
  return sum;
}

double Multifind_object::average_multifind(vector<int>& series){
  return (sum_multifind(series))/series.size();
}


double Multifind_object::variation_multifind(vector<int>& series){
  double average=average_multifind(series);
  double sum=0;
  for(int i=0;i<series.size();++i){
    sum+=(series[i]-average)*(series[i]-average);
  }
  return sqrt(sum/(series.size()-1));
}




double Multifind_object::normalized_ensemble_defect(RNA* rna){
//    cout <<rna->GetSequenceLength()<<"\n";
  double defect=0.0;
  for(int i=1;i<=rna->GetSequenceLength();++i){
    if(rna->GetPair(i)==0){
      for(int j=1;j<=rna->GetSequenceLength();++j){
	if(i<j)defect+=rna->GetPairProbability(i,j);
	else if(j<i) defect+=rna->GetPairProbability(j,i);
      }
    }
    else if(rna->GetPair(i)>i)defect+=2*(1-rna->GetPairProbability(i,rna->GetPair(i)));
  }
  return defect/rna->GetSequenceLength();
  
}


double Multifind_object::entropy(vector<char> column){
  int A=0;
  int G=0;
  int U=0;
  int C=0;
  int gap=0;
  for(int i=0;i<column.size();++i){
    if(column[i]=='A')A++;
    if(column[i]=='G')G++;
    if(column[i]=='C')C++;
    if(column[i]=='U'||column[i]=='T')U++;
    if(column[i]=='-')gap++;
  }
  double pA=double(A)/(A+G+C+U+gap);
  // cout<<A<<"\n";
  //  cout<<A+G+C+U+gap<<"\n";
  double pU=double(U)/(A+G+C+U+gap);
  double pG=double(G)/(A+G+C+U+gap);
  double pC=double(C)/(A+G+C+U+gap);
  double pgap=double(gap)/(A+G+C+U+gap);
  
  double entropy=0;

  if(pA!=0)entropy+=pA*log(pA);
  if(pG!=0)entropy+=pG*log(pG);
  if(pC!=0)entropy+=pC*log(pC);
  if(pU!=0)entropy+=pU*log(pU);
  if(pgap!=0)entropy+=pgap*log(pgap);
 
  return entropy;
}


double Multifind_object::average_entropy(){
  vector<double> entropys;
  for(int i=0;i<input_alignment[0].size();++i){
    vector<char> column;
    for(int j=0;j<input_alignment.size();++j){
      column.push_back(input_alignment[j][i]);
    }
    //   cout<<entropy(column)<<"\n";
    entropys.push_back(entropy(column));
  }

  return average_multifind(entropys);
}
  

vector<double> Multifind_object::single_z_predict(string sequence,svm_model* model_folding_average,svm_model* model_folding_std,svm_model* model_ensemble_average,svm_model* model_ensemble_std){
    
  vector<double> results;
  int A=0;
  int G=0;
  int U=0;
  int C=0;
  for(int i=0;i<sequence.size();++i){
    if(sequence[i]=='A')A++;
    if(sequence[i]=='G')G++;
    if(sequence[i]=='C')C++;
    if(sequence[i]=='U'||sequence[i]=='T')U++;
  }
  


  double GC_content=double(G+C)/sequence.size();
  double G_content=double(G)/(G+C);
  double A_content=double(A)/(A+U);
  double length=double(sequence.size());
  
  double scaled_GC_content=(GC_content-0.25)/0.5*2-1;
  double scaled_G_content=(G_content-0.25)/0.5*2-1;
  double scaled_A_content=(A_content-0.25)/0.5*2-1;
  double scaled_length=(length-30)/120*2-1;


  struct svm_node* x_svm=new svm_node[5];
  x_svm[0].index=1;
  x_svm[0].value=scaled_GC_content;
  x_svm[1].index=2;
  x_svm[1].value=scaled_G_content;
  x_svm[2].index=3;
  x_svm[2].value=scaled_A_content;
  x_svm[3].index=4;
  x_svm[3].value=scaled_length;
  x_svm[4].index=-1;
  double average=svm_predict(model_folding_average,x_svm);
  double std=svm_predict(model_folding_std,x_svm);
  delete[] x_svm;  

  average=(average+1)/2*(-0.1638+80.40180000000009)-80.40180000000009;
  std=std*(4.81792216250574-0.514970634847032)+0.514970634847032;


  results.push_back(average);
  results.push_back(std);

  // vector<double> results;
  scaled_GC_content=(GC_content-0.233333333333333)/(0.775-0.233333333333333)*2-1;
  scaled_G_content=(G_content-0.222222222222222)/(0.8-0.222222222222222)*2-1;
  scaled_A_content=(A_content-0.2)/0.675*2-1;
  scaled_length=(length-30)/120*2-1;

  struct svm_node* y_svm=new svm_node[5];
  y_svm[0].index=1;
  y_svm[0].value=scaled_GC_content;
  y_svm[1].index=2;
  y_svm[1].value=scaled_G_content;
  y_svm[2].index=3;
  y_svm[2].value=scaled_A_content;
  y_svm[3].index=4;
  y_svm[3].value=scaled_length;
  y_svm[4].index=-1;
  double ensemble_average=svm_predict(model_ensemble_average,y_svm);
  double ensemble_std=svm_predict(model_ensemble_std,y_svm);
  delete[] y_svm;  

  results.push_back(ensemble_average);
  results.push_back(ensemble_std);
  return results;
}


int Multifind_object::get_gap(string sequence){
  
  for(int i=0;i<sequence.size()&&sequence[i]=='-';++i){
    sequence[i]=' ';
  }
    for(int i=sequence.size()-1;i>=0&&sequence[i]=='-';--i){
    sequence[i]=' ';
  }
  int number=0;
  for(int i=0;i<sequence.size();++i){
    if(sequence[i]=='-')number++;
  }
  return number;
}

/*
string Multifind_object::compact(string& sequence){
  string temp_sequence="";
  for(int i=0;i<sequence.size();i++)
    if(sequence[i]!='-'){temp_sequence+=sequence[i];}
  return temp_sequence;
}
*/  

double Multifind_object::common_energies(){
  
  vector<float> gapIndex;
  vector<float> gaps;
  for(int i=0;i<pair_alignments.size();++i){
    //     cout<<pair_alignments[i][0]<<"\n";
    // cout<<pair_alignments[i][1]<<"\n";
    gapIndex.push_back(get_gap(pair_alignments[i][0]));
    gaps.push_back(get_gap(pair_alignments[i][1]));
  }
  
  for(int i=0;i<pair_alignments.size();++i){
    // cout<<gapIndex[i]<<"\n";
    // cout<<dGIndex[i]<<"\n";
    //cout<<gaps[i]<<"\n";
    //cout<<energies[i]<<"\n";
  }
 
  double sum_common_energies=average_multifind(dGIndex)+0.4*average_multifind(gapIndex)+sum_multifind(energies)+0.4*sum_multifind(gaps);
   return sum_common_energies;

 
}


vector<double> Multifind_object::predict_ncRNA_probabilities(double sci,double entropy,double single_z,double ensemble_defect_z){
  string path=string(getDataPath());
		     
  string multi_model_name=path+"/"+"data_assemble_training_Multifind_predict_ensemble_z_final_svmformat.model";
 
  struct svm_node* multi=new svm_node[5];
  struct svm_model* multi_model=svm_load_model(multi_model_name.c_str());
  int nr_class_multi=svm_get_nr_class(multi_model);
  double scaled_sci=(sci+4.58621)/5.58621*2-1;
  double scaled_entropy=(entropy+0.958323)/0.958323*2-1;
  double scaled_single_z=(single_z+17.0798)/(1.9982+17.0798)*2-1;
  double scaled_ensemble_defect_z=(ensemble_defect_z+2.49514)/(1.98099+2.49514)*2-1;
  multi[0].index=1;
  multi[0].value=scaled_sci;
  multi[1].index=2;
  multi[1].value=scaled_entropy;
  multi[2].index=3;
  multi[2].value=scaled_single_z;
  multi[3].index=4;
  multi[3].value=scaled_ensemble_defect_z;
  multi[4].index=-1;

  double* multi_prob_estimates=NULL;
  multi_prob_estimates=(double *)malloc(nr_class_multi*sizeof(double));
  svm_predict_probability(multi_model,multi,multi_prob_estimates);

  vector<double> multi_probabilities;
  for(int j=0;j<nr_class_multi;j++){
    multi_probabilities.push_back(multi_prob_estimates[j]);
  }
  
  svm_free_and_destroy_model(&multi_model);
  delete[] multi_prob_estimates;
  delete[] multi;
  return multi_probabilities;
}


Multifind_object::Multifind_object(const string &outputmultifind, const vector<string> &ctfiles, const vector<string> &inputalignment, const vector<string> &inputsequences, const int &processors,ProgressHandler *Progress): 
  Multilign_object(true,outputmultifind,ctfiles,Progress),num_processors(processors){

  input_alignment=inputalignment;
  input_sequences=inputsequences;
  
}

int Multifind_object::Multifind_Predict(){
  vector<double> single_energies;
  vector<double> ensemble_defect_zs;
  vector<double> single_zs;

  if(ErrorCode=ProgressiveMultilign(num_processors)){
    string error_message=GetErrorMessage(ErrorCode);
    cerr<<error_message;
    exit(EXIT_FAILURE);
  }

  if(ErrorCode=WriteAlignment(output_multifind)){
    string error_message=GetErrorMessage(ErrorCode);
    cerr<<error_message;
    exit(EXIT_FAILURE);
  }
  
  float sum_common_energies=common_energies();
  //  cerr<<"lala!\n";
  string path=string(getDataPath());

  string file_name=path+"/"+"new_training_z_ave.scale.model";
  struct svm_model* model_folding_average=svm_load_model(file_name.c_str());
  file_name=path+"/"+"new_training_z_std.scale.model";
  struct svm_model* model_folding_std=svm_load_model(file_name.c_str());
  file_name=path+"/"+"average_ensemble_defect.model";
  struct svm_model* model_ensemble_average=svm_load_model(file_name.c_str());
  file_name=path+"/"+"std_ensemble_defect.model";
  struct svm_model* model_ensemble_std=svm_load_model(file_name.c_str());
  
  for(int i=0;i<input_sequences.size();++i){
    vector<double> single_predict_results=single_z_predict(input_sequences[i],model_folding_average,model_folding_std,model_ensemble_average,model_ensemble_std);
    RNA* point_RNA=new RNA(input_sequences[i].c_str());
    point_RNA->FoldSingleStrand(0,1,0);
    double single_energy=point_RNA->GetFreeEnergy(1);
    single_energies.push_back(single_energy);
    double single_z=(single_energy-single_predict_results[0])/(single_predict_results[1]);
    single_zs.push_back(single_z);

    point_RNA->PartitionFunction();
    point_RNA->MaximizeExpectedAccuracy(0,1,0);
    double ensemble_defect=normalized_ensemble_defect(point_RNA);
//    cout <<ensemble_defect<<" "<<single_predict_results[2]<<" "<<single_predict_results[3]<<"\n";
    double ensemble_defect_z=(ensemble_defect-single_predict_results[2])/(single_predict_results[3]);
    ensemble_defect_zs.push_back(ensemble_defect_z);
  }
  
  svm_free_and_destroy_model(&model_folding_average);
  svm_free_and_destroy_model(&model_folding_std);
  svm_free_and_destroy_model(&model_ensemble_average);
  svm_free_and_destroy_model(&model_ensemble_std);

  double sci=sum_common_energies/(sum_multifind(single_energies));
  ofstream OUTPUT(output_multifind.c_str(),ios::app);
  
  OUTPUT<<"\nSCI:"<<sci<<"\n";
  double Average_Entropy=average_entropy();
  OUTPUT<<"Average Entropy:"<<-1*Average_Entropy<<"\n";
  double Average_Z_Score=average_multifind(single_zs);
  OUTPUT<<"Average Single Z Score:"<<Average_Z_Score<<"\n";
  double Average_Ensemble_Defect_Z=average_multifind(ensemble_defect_zs);
  OUTPUT<<"Average Ensemble Defect Z Score:"<<Average_Ensemble_Defect_Z<<"\n";
 
  vector<double> ncRNA_probabilities=predict_ncRNA_probabilities(sci,Average_Entropy,Average_Z_Score,Average_Ensemble_Defect_Z);

  OUTPUT<<"Multifind Model Probability:"<<ncRNA_probabilities[0]<<"\n";
  OUTPUT<<"Multifind Model Non Probability:"<<ncRNA_probabilities[1]<<"\n";
  
  return 1;
}
