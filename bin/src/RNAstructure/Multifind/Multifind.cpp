#include "Multifind.h"


Multifind_Interface::Multifind_Interface(){
  processors=1;
}

string Multifind_Interface::compact(string& sequence){
  string temp_sequence="";
  for(int i=0;i<sequence.size();i++)
    if(sequence[i]!='-'){temp_sequence+=sequence[i];}
  return temp_sequence;
}
   

void Multifind_Interface::read_fasta(){
   ifstream in_fasta;
  in_fasta.open(fasta_file.c_str());
  string temp_buffer;
  //  string::size_type alignment_length;
  vector<string> content;
  while(in_fasta){
    getline(in_fasta,temp_buffer);
    if(temp_buffer.find_last_of('\n')!=string::npos)temp_buffer=temp_buffer.substr(0,temp_buffer.size()-1);
    if(temp_buffer!="")content.push_back(temp_buffer);
  }

  /*  cerr<<content.size()<<"size\n";

  for(int i=0;i<content.size();i++){
    cerr<<i<<" "<<content[i]<<"\n";
  }
  */
  for(vector<string>::size_type i=0;i<content.size();){
    
    string temp_sequence;
    
   
    if(content[i][0]=='>'){
      i++;
      while(i<content.size()){
	temp_sequence+=content[i];
	i++;
	if(i>=content.size())break;
	if(content[i][0]=='>')break;
      }
      input_alignment.push_back(temp_sequence);
    }
   
  }
  /*
  if(input_alignment.size()==0){
    cerr<<"There is no sequences in the input file. Please check the format of the input file. The input must be in FASTA format.\n";
      exit(EXIT_FAILURE);
      }*/

  for(vector<string>::iterator it=input_alignment.begin();it!=input_alignment.end();++it){
    string temp_sequence=compact(*it);
     input_sequences.push_back(temp_sequence);
  }
  
  number_seq=input_sequences.size();
  /*  for(int i=0;i<input_alignment.size();i++){
    cerr<<i<<"\n"<<input_alignment[i]<<"\n";
    cerr<<input_sequences[i]<<"\n";
  }
  */
}

bool Multifind_Interface::parse(int argc,char** argv){
#ifndef COMPILE_SMP
  ParseCommandLine* parser=new ParseCommandLine("Multifind");
#else
  ParseCommandLine* parser=new ParseCommandLine("Multifind-smp");
#endif

  //  int number_seq=0;
  /*
  ParseCommandLine* preparser=new ParseCommandLine("Multifind");
  vector<string> structureoutputOptions;
  structureoutputOptions.push_back("-s");
  structureoutputOptions.push_back("-S");
  structureoutputOptions.push_back("--structures");
  preparser->addOptionFlagsWithParameters( structureoutputOptions, "Specify if the user wants to output structures" );
  if( !preparser->isError() ) {
    preparser->setOptionInteger( structureoutputOptions, structure_output);
    if(structure_output!=0&&structure_output!=1){preparser->setError("the structures option should be set as 1 or 0");}
  }
  delete preparser;
  */
  /*  if(!structure_output){
    cerr<<"lala\n";
    parser->addParameterDescription( "fasta file", "The name of a file containing an multiple sequence alignment." );
    parser->addParameterDescription( "multifind file", "The name of a file to which the multifind output is written." );

#ifdef COMPILE_SMP
    vector<string> processorsOptions;
    processorsOptions.push_back("-p");
    processorsOptions.push_back("-P");
    processorsOptions.push_back("--processors");
    parser->addOptionFlagsWithParameters( processorsOptions, "Specify the number of processors in smp calculation" );
#endif

    parser->parseLine( argc, argv );
  
    if( !parser->isError() ) {
      fasta_file=parser->getParameter(1);
      multifind_file=parser->getParameter(2);
    }

#ifdef COMPILE_SMP
      if( !parser->isError() ) {
	parser->setOptionInteger( processorsOptions, processors);
	if(processors<=0){parser->setError("the number of processors can't be negative or zero");}
      }
#endif      
    read_fasta();
    bool noError = ( parser->isError() == false );
    delete parser;
    return noError;
  }
  else{
  */

    parser->addParameterDescription( "fasta file", "The name of a file containing an multiple sequence alignment." );
    parser->addParameterDescription( "multifind file", "The name of a file to which the multifind output is written." );

    if(argc>=2){
      fasta_file=string(argv[1]);
      read_fasta();
    }
    else{
      cerr<<"Please provide a FASTA alignment file as the first input\n";
      exit(EXIT_FAILURE);
    }
    for(int i=1;i<=number_seq;i++){
      stringstream ss;
      ss<<i;
      string str=ss.str();
      parser->addParameterDescription( "ct file "+str, "The name of a file to which the structure output is written." );
    }
 
#ifdef COMPILE_SMP
     vector<string> processorsOptions;
    processorsOptions.push_back("-p");
    processorsOptions.push_back("-P");
    processorsOptions.push_back("--processors");
    parser->addOptionFlagsWithParameters( processorsOptions, "Specify the number of processors in smp calculation" );
#endif      
   
    parser->parseLine( argc, argv );
    if( !parser->isError() ) {
      fasta_file=parser->getParameter(1);
      multifind_file=parser->getParameter(2);
      for(int i=1;i<=number_seq;i++){
	string ct_file=parser->getParameter(i+2);
	ct_files.push_back(ct_file);
      }
    }
#ifdef COMPILE_SMP
  if( !parser->isError() ) {
	parser->setOptionInteger( processorsOptions, processors);
	if(processors<=0){parser->setError("the number of processors can't be negative or zero");}
      }
#endif      
    bool noError = ( parser->isError() == false );
  delete parser;
  return noError;
  //  }

  


}

void Multifind_Interface::run(){
  TProgressDialog* progress = new TProgressDialog();
  Multifind_object* instance=new Multifind_object(multifind_file,ct_files,input_alignment,input_sequences,processors,progress);
  
  instance->Multifind_Predict();
  
  instance->StopProgress();
  delete progress;
  instance->CleanupIntermediateFiles();
  delete instance;	
  
  // return 1;
}

int main(int argc,char* argv[]){
 
  Multifind_Interface* runner=new Multifind_Interface();
  
  bool parseable = runner->parse(argc,argv);
  if(parseable==true){runner->run();}
  delete runner;
  return 0;
}

