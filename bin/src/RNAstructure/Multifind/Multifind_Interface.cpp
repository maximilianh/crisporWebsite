#include "Multifind_Interface.h"


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


  parser->addParameterDescription( "configuration file", "The name of a file containing configuration data." );


  // Parse the command line into pieces.
  parser->parseLine( argc, argv );


  // Open the config file.
  // If the file isn't valid, delete the parser and return false.
  ConfigFile file( parser->getParameter( 1 ) );
  if( file.isValid() == false ) {
      delete parser;
      return false;
  }
 
  
	// Check if the config file contains the Multifind flag, which is always necessary for reading data.
	// If it exists, read it in.
	// If the flag isn't present, set an error, delete the parser, and return false.
	bool isReadable = 
            file.contains( "Multifind" ) &&
            file.contains( "Fasta" );

	if( isReadable ) { 
            multifind_file = file.getOption<string>( "Multifind" ); 
            fasta_file = file.getOption<string>( "Fasta" );
            read_fasta();
        }
	else {
		parser->setErrorSpecialized( "Configuration file must contain a Multifind and a Fasta flag." );
		delete parser;
		return false;
	}


	// If SMP calculations are being done, make sure the processors flag has been specified.
	// If it has, set the number of processors if possible.
	// If it hasn't, show an error, delete the parser, and return false.
#ifdef COMPILE_SMP
	if( file.contains( "Processors" ) ) {
		processors = file.getOption<int>( "Processors" );
		if( processors < 0 ) {
			parser->setError( "number of processors" );
			delete parser;
			return false;
		}
	} else {
		parser->setErrorSpecialized( "Configuration file must contain the Processors flag." );
		delete parser;
		return false;
	}
#else
#endif

     isReadable =
         file.contains( "OutCT" );
     if ( isReadable ) {
         string ctData = file.getOption<string>( "OutCT" );
         unsigned int ctLast = ctData.length() - 1;
         if( ctData[0] == '{' && ctData[ctLast] == '}' ) {
             ctData = ctData.erase( 0, 1 );
             ctData = ctData.erase( ctLast - 1, 1 );
             ctLast = ctData.length() - 1;
             if( ctData[ctLast] == ';' ) { ctData = ctData.erase( ctLast, 1 ); }
             stringstream ctStr( ctData );
             string ctFile;
             while( ctStr.good() ) {
                 getline( ctStr, ctFile, ';' );
                 if( ctFile != "" ) { ct_files.push_back( ctFile ); }
             }
         } else {
             if( ctData[0] != '{' ) { parser->setErrorSpecialized( "OutCT group has no start bracket." ); }
             else { parser->setErrorSpecialized( "OutCT group has no end bracket." ); }
             delete parser;
             return false;
         }
         
         if ( ct_files.size() != number_seq ) {
             parser->setErrorSpecialized( "OutCT group does not have equal number of sequences as Fasta." ); 
             delete parser;
             return false;
         }
     }     
     
     

     // Delete the parser and return whether the parser encountered an error.
     bool noError = ( parser->isError() == false );
     delete parser;
     return noError;

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

