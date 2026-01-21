#include "../src/pfunction.h"
#include "ProbScan.h"
#include "../src/structure.h"
#include <algorithm>
#include <functional>
#include <numeric>
#include <iomanip>
#include <assert.h>
#include <cmath>
using std::vector;
using std::string;


const static int max_internal_loop=30;

ProbScan::ProbScan(std::string sequence, const char* const alphabetName):RNA(sequence.c_str(), alphabetName)
{
    PartitionFunction();
}

ProbScan::ProbScan(const char filename[], bool from_sequence_file, const char* const alphabetName):RNA(filename,from_sequence_file?FILE_SEQ:FILE_PFS, alphabetName)
{
  if (from_sequence_file)
    PartitionFunction();//calculate the partition function if it hasn't been
                        //done yet
}

double ProbScan::probability_of_hairpin(int i,int j)
{
  return TO_LINEAR( DIV(PROD(v->f(j,i+GetSequenceLength()) //V'(i,j)
         , erg3(i,j,GetStructure(),pfdata,0)) //K for hairpin
         , PFSCALE(w5[GetSequenceLength()],pfdata->scaling,2))); //Q
//divide by scaling^2 so the closing nucs aren't double counted
}

double ProbScan::probability_of_stemloop(int i, int j,int ip, int jp)
{
	//First calculate the equilibrium constant for the helix
	double helix_product = ONE;
	for (int position = 0; position < ip - i; ++position) {
		helix_product = PROD(helix_product, erg1(i+position,j-position,i+1+position,j-1-position,GetStructure(),pfdata));
	}

	return TO_LINEAR(DIV(PROD(v->f(j, i + GetSequenceLength()) //V'(ip,jp)
		, erg3(ip, jp, GetStructure(), pfdata, 0) //K for hairpin
		, helix_product)
		, PFSCALE(w5[GetSequenceLength()], pfdata->scaling, 2))); //Q
//divide by scaling^2 so the closing nucs aren't double counted
}

vector<hairpin_t> ProbScan::probability_of_all_hairpins(int min_size, int max_size,double threshold)
{
  vector<hairpin_t> hairpins;
  structure* st = GetStructure();
//#pragma omp parallel for
//parallelization causes datarace because of push_back method
  for(int i=1;i<GetSequenceLength()-min_size-1;i++){//search over all 0<i<j<n
    for(int j=i+min_size+1;j<std::max(i+max_size,GetSequenceLength());j++){
      if(GetDatatable()->can_pair(i,j,GetStructure()->numseq)){//if i and j can pair
        //get probability
        double probability = probability_of_hairpin(i,j);
        if (probability>threshold){ //add to the list if p>threshold
          hairpins.push_back(hairpin(probability,i,j));
        }
      }
    }
  }
  //hairpins now contains every hairpin where p>threshold
  return hairpins;
}

double ProbScan::probability_of_internal_loop(int i,int j, int k, int l)
{
  return TO_LINEAR(DIV(PROD(v->f(k,l) //V(k,l)
         , v->f(j,i+GetSequenceLength()) //V'(i,j)
         , erg2(i,j,k,l,GetStructure(),pfdata,0,0))//K for iloop
         , PFSCALE(w5[GetSequenceLength()],pfdata->scaling,2)));//Q
//we divide by scaling^2 so the closing nucs aren't double counted
 //previously, performed a correction for slipping by isoenergetic single nucleotide bulges
//but it was wrong. calculated probs have to be explicitly summed for this to be correct.
}


basestack_t basestack(double p,int i, int j, int k, int l){
    basestack_t stk;
    stk.i = i;
    stk.j = j;
    stk.k = k;
    stk.l = l;
    stk.probability = p;
    return stk;
}
double ProbScan::probability_of_stack(int i,int j)
{
  int k = i+1;
  int l = j-1;
  return TO_LINEAR(DIV(PROD(v->f(k,l) //V(k,l)
         , v->f(j,i+GetSequenceLength()) //V'(i,j)
         , erg1(i,j,k,l,GetStructure(),pfdata))//K for stack
         , PFSCALE(w5[GetSequenceLength()],pfdata->scaling,2)));//Q
}

double ProbScan::probability_of_helix(const int i, const int j, const int how_many_stacks){
    assert(how_many_stacks > 0);
    structure* st = GetStructure();
    if((j-i < how_many_stacks * 2 + minloop + 1) || (!GetDatatable()->can_pair(i, j, GetStructure()->numseq))){
        return ZERO;//no stack if i and j are too close or can't pair
    }
    PFPRECISION kStack = ONE;
    const int k = i+1;
    const int l = j-1;
    for(int s = 0; s < how_many_stacks; s++){
        if(!GetDatatable()->can_pair(k + s, l - s, GetStructure()->numseq)){//no stack if nucs can't pair
            return ZERO;
        }
        kStack = PROD(kStack, erg1(i+s,j-s,k+s,l-s,st,pfdata));
    }
    return TO_LINEAR(DIV(PROD(v->f(i+how_many_stacks,j-how_many_stacks) //V(k,l)
         , v->f(j,i+GetSequenceLength()) //V'(i,j)
         , kStack)//K for helix
         , PFSCALE(w5[GetSequenceLength()],pfdata->scaling,2)));//Q
}

std::vector<basestack_t> ProbScan::probability_of_all_helices(double threshold,int length)
{
    std::vector<basestack_t> stacks = std::vector<basestack_t>();
    for(int i=1;i<GetSequenceLength();i++){
        for(int j=i+minloop+1+length*2;j<GetSequenceLength();j++){
            double probability = probability_of_helix(i,j,length);
//            printf("%0.4f\n",probability);
            if(probability>threshold){
                stacks.push_back(basestack(probability,i,j,i+length,j-length));
            }
        }
    }
    return stacks;
}

double ProbScan::probability_of_multibranch_loop(const multibranch_loop_t& mb)
{
  assert(mb.branches.size()>=3);
  //holds the values from v array
  vector<PFPRECISION> vs;
  //V(j,i+numberofbases) for closing pair
  vs.push_back(PROD(v->f(mb.branches[0].second,mb.branches[0].first+GetSequenceLength())
                   ,penalty(mb.branches[0].second,mb.branches[0].first,GetStructure(),pfdata)));
  //V(i,j) for each branch, with AU/GU end penalty
  for(vector<std::pair<int,int> >::const_iterator it=mb.branches.begin()+1;it!=mb.branches.end();++it){
    vs.push_back(PROD(v->f(it->first,it->second)
                   ,penalty(it->first,it->second,GetStructure(),pfdata)));
  }
  //calculate equilibrium constant
  PFPRECISION Kmb = equilibrium_constant_for_multibranch_loop(mb);
  //take product of values from V array
#ifdef PF_LOG_CALC  
  PFPRECISION product_of_vs = std::accumulate(vs.begin(),vs.end(),(PFPRECISION) ONE);
#else
  PFPRECISION product_of_vs = std::accumulate(vs.begin(),vs.end(),(PFPRECISION) ONE
                                     ,std::multiplies<PFPRECISION>());
#endif  
  //return probability
  return TO_LINEAR(DIV(PROD(Kmb , product_of_vs) , w5[GetSequenceLength()]));
}

vector<internal_loop_t> ProbScan::probability_of_all_internal_loops(double threshold,std::string mode/*="both"*/)
{
  vector<internal_loop_t> iloops;//holds internal loops that we find
  int n = GetSequenceLength();
  structure* st = GetStructure();
  bool bulges_allowed = (mode == string("bulge")) || (mode == string("both"));
  bool iloops_allowed = (mode == string("internal")) || (mode == string("both"));
//#pragma omp parallel for //can't parallelize because of push_back
//search over all i,j,k,l with < max_internal_loop unpaired nucs
  for(int i=1;i<n-3;i++){
    //if we are considering iloops only(bulges_allowed is false), k-i>1
    int kmin = i + (bulges_allowed?1:2);
    //loops with up to 30 unpaired nucs
    int kmax = std::min(i+max_internal_loop,n-2);
    for(int k=kmin;k<kmax+1;k++){
      for(int l=k+minloop+1;l<n-1;l++){
        int jmin = 0;
        if(k-i == 1){
            //if k is next to i then the bulge is on the 3' side
            jmin = l+2;
        }
        else {
            //if we are considering iloops only (bulges_allowed is false), j-l>1
            jmin = l+(bulges_allowed?1:2);
        }
        //limit to 30 unpaired nucs, which depends on k-i
        int iloopjmax = std::min(l+(max_internal_loop-(k-i+1)),n);
        int jmax;
        if (iloops_allowed){
            jmax = iloopjmax;
        }
        //if we are restricted to bulge loops, then j-l must be 1 if k-i>1
        else {
            jmax = (k-i==1)?iloopjmax:l+1;
        }
        for(int j=jmin;j<jmax+1;j++){
          if(!GetDatatable()->can_pair(k, l, GetStructure()->numseq)) break;//if i can't pair to k we're done
          //if i can pair to j and k can pair to l
          if(GetDatatable()->can_pair(i, j, GetStructure()->numseq)){
            assert(i<k && l<j && k<l);//input validation
            assert(k-i>1 || j-l>1);
            //get probability of the internal loop
            double probability=probability_of_internal_loop(i,j,k,l);
            if (probability>threshold) {//add to list if prob>threshold
              iloops.push_back(internal_loop(probability,i,j,k,l));
            }
          }
        }
      }
    }
  }
  return iloops;//vector now holds all possible iloops with p>threshold
}


vector<mb_element> ProbScan::construct_mb_element_array(const multibranch_loop_t& mb)
{
//construct the mb_element array containing hairpins and nucleotides
//first make the closing hairpin,
//swapping its indices so its the same orientation as the others
  vector<mb_element> mb_element_array;
  mb_element closing_hairpin = mb_element(std::make_pair(mb.branches[0].second,mb.branches[0].first));
  mb_element_array.push_back(closing_hairpin);
  bool first=true;
  bool last=false;
//mb is a vector of pairs
//for each pair in the multibranch loop after the first one, add any unpaired
//nucleotides between the last pair and this pair, then add the new pair
  for(vector<std::pair<int,int> >::const_iterator it=mb.branches.begin()+1;it!=mb.branches.end();++it){
    mb_element last_pair = mb_element(*(it-1));//get last pair
    mb_element next_pair = mb_element(*it);//and next pair
    assert(last_pair.is_a_pair && next_pair.is_a_pair);//input validation
    for(int x = first?last_pair.i+1:last_pair.j+1;x<next_pair.i;x++){
      mb_element_array.push_back(mb_element(x));//insert unpaired nucs between
                                          //last_pair and next_pair
    }
    first=false;
    mb_element_array.push_back(next_pair); //add next_pair, which will be the new
                                        //last_pair
  }
//add the last few unpaired nucs
  for(int x = mb_element_array.back().j+1;x<mb_element_array.front().i;x++){
    mb_element_array.push_back(mb_element(x));
  }

  for(int i=0;i<4;i++){//duplicate the first 4 mb_elements of the array on the end
    mb_element_array.push_back(mb_element_array[i]);
  }
  return mb_element_array;
}

PFPRECISION prev_val(int index,int offset, vector<vector<PFPRECISION> >& arr)
{
  const PFPRECISION initial_value = ONE;
  if(index>=offset) return arr[index][offset];
  else return initial_value;
}

PFPRECISION ProbScan::equilibrium_constant_for_multibranch_loop(const multibranch_loop_t& mb)
{
  vector<mb_element> mb_elements = construct_mb_element_array(mb);
  if(mb_elements.size()<=8) return ZERO;//Rahul did this.. need to think more about justification
//N+4 by 4 array for accumulating the partition function
  vector<vector<PFPRECISION> > arr(mb_elements.size(),vector<PFPRECISION>(4,ONE));
  short* s = GetStructure()->numseq;//the nucleotide sequence
  const int n = mb_elements.size()-4;
  int stems=0,unpaired_nucs=0;
//first calculate the partition function around the circle starting at 4
//different start points "offsets". we'll handle interactions across the //starting points later
  for(int offset=0;offset<4;offset++){
    for(int x=offset;x<n+offset;x++){
      PFPRECISION result = ZERO;
      if(mb_elements[x].is_a_pair){
        if(offset==0) stems += 1;
  //case where new helix has no coaxial stack or danging end
        result = SUM(result, prev_val(x-1,offset,arr));

  //case where the previous nuc is making a 5' dangle on new helix
        if(!mb_elements[x-1].is_a_pair && x>offset){
          result = SUM(result, PROD(prev_val(x-2,offset,arr), erg4(mb_elements[x].j,mb_elements[x].i,mb_elements[x-1].i,2,
                                  GetStructure(),pfdata,0)));
        }
#ifndef disablecoax
  //case where new helix is coaxially stacking flush with helix at mb_element x-1
        if(mb_elements[x-1].is_a_pair && x>offset){
          result = SUM(result, PROD(prev_val(x-2,offset,arr) ,
                     ergcoaxflushbases(mb_elements[x-1].i,mb_elements[x-1].j,
                                       mb_elements[x].i,mb_elements[x].j,
                                       GetStructure(),pfdata)));
        }
  //case where new helix is coaxially stacking on mb_element x-2 with an intervening mismatch
  //coax 1: , || , ||
        if(x>2+offset && mb_elements[x-2].is_a_pair &&
            !mb_elements[x-1].is_a_pair && !mb_elements[x-3].is_a_pair){
          result = SUM(result, PROD(prev_val(x-4,offset,arr) ,
                     ergcoaxinterbases1(mb_elements[x-2].i,mb_elements[x-2].j,
                                        mb_elements[x].i,mb_elements[x].j,
                                        GetStructure(),pfdata)));
        }
#endif //disablecoax
      }
      else if(!mb_elements[x].is_a_pair){
        if(offset==0) unpaired_nucs += 1;
  //case where nuc is not dangling or participating in a mismatch
  //we have to scale for the unpaired nuc
        result = SUM(result, PFSCALE(prev_val(x-1,offset,arr),pfdata->scaling,1));
  //case where nuc is dangling 3' on helix at mb_element x-1
        if(x>offset && mb_elements[x-1].is_a_pair){
          result = SUM(result, PROD(prev_val(x-2,offset,arr), erg4(mb_elements[x-1].j,mb_elements[x-1].i,mb_elements[x].i,1,
                                  GetStructure(),pfdata,0)));
        }
  //case where mb_element x and x-2 are nucs forming a mismatch on helix at x-1
  // , || ,
        if(x>1+offset && mb_elements[x-1].is_a_pair && !mb_elements[x-2].is_a_pair){
          result = SUM(result, PROD(prev_val(x-3,offset,arr) ,
                    pfdata->tstkm[s[mb_elements[x-1].j]][s[mb_elements[x-1].i]]
                               [s[mb_elements[x].i]][s[mb_elements[x-2].i]]));

        }
  //case where new nuc is participating in a mismatch coaxial stack between helices at x-1 and x-3
  //coax 2: || , || ,
#ifndef disablecoax
        if(x>2+offset && !mb_elements[x-2].is_a_pair &&
            mb_elements[x-1].is_a_pair && mb_elements[x-3].is_a_pair){
          result = SUM(result, PROD(prev_val(x-4,offset,arr) ,
                     ergcoaxinterbases2(mb_elements[x-3].i,mb_elements[x-3].j,
                                        mb_elements[x-1].i,mb_elements[x-1].j,
                                        GetStructure(),pfdata)));
        }
#endif //disablecoax
      }
      arr[x][offset] = result;
    }
  }
//now let's paste the ends of the sequence together
//there are only a few cases because mb_element 0 is always
//a helix in current implementation
  PFPRECISION initiation = PROD(POWER(pfdata->eparam[10],stems) //per stem penalty "c"
                  , POWER(pfdata->eparam[6],unpaired_nucs) //per nuc penalty "b"
                  , pfdata->eparam[5]);//closing multibranch loop penality "a"
  PFPRECISION pfunc = arr[n-1][0];
#ifndef disablecoax
  if(mb_elements[0].is_a_pair && mb_elements[n-1].is_a_pair){
    pfunc = SUM(pfunc, PROD(arr[n-2][1], ergcoaxflushbases(mb_elements[n-1].i,mb_elements[n-1].j,
                                       mb_elements[0].i,mb_elements[0].j,
                                       GetStructure(),pfdata)));
  }
#endif//disablecoax
  if(mb_elements[0].is_a_pair && !mb_elements[n-1].is_a_pair){//5' dangle
    pfunc = SUM(pfunc, PROD(arr[n-2][1],erg4(mb_elements[0].j,mb_elements[0].i,mb_elements[n-1].i,2,
                                  GetStructure(),pfdata,0)));
  }
//terminal mismatch
  if(mb_elements[0].is_a_pair && !mb_elements[1].is_a_pair && !mb_elements[n-1].is_a_pair){
    pfunc = SUM(pfunc, PROD(arr[n-2][2],pfdata->tstkm[s[mb_elements[n].j]][s[mb_elements[n].i]]
                               [s[mb_elements[n+1].i]][s[mb_elements[n-1].i]]));
   }
//3 mismatch coax possibilities
// case like this: 5' || , ...  || ,  3'
#ifndef disablecoax
  if(mb_elements[0].is_a_pair && mb_elements[n-2].is_a_pair &&
     !mb_elements[1].is_a_pair && !mb_elements[n-1].is_a_pair){
    pfunc = SUM(pfunc, PROD(arr[n-3][2] , ergcoaxinterbases2(mb_elements[n-2].i,mb_elements[n-2].j,
                                        mb_elements[n].i,mb_elements[n].j,
                                        GetStructure(),pfdata)));
  }
// case like this: 5' || ... , || ,  3'
  if(!mb_elements[n-3].is_a_pair && !mb_elements[n-1].is_a_pair &&
     mb_elements[n-2].is_a_pair && mb_elements[n].is_a_pair){
    pfunc = SUM(pfunc, PROD(arr[n-4][1] , ergcoaxinterbases1(mb_elements[n-2].i,mb_elements[n-2].j,
                                      mb_elements[n].i,mb_elements[n].j,
                                      GetStructure(),pfdata)));
    }
// case like this: 5' || , ||   ...   ,  3'
  if(!mb_elements[n-1].is_a_pair && !mb_elements[n+1].is_a_pair &&
     mb_elements[n].is_a_pair && mb_elements[n+2].is_a_pair){
    pfunc = SUM(pfunc, PROD(arr[n-2][3] , ergcoaxinterbases1(mb_elements[n].i,mb_elements[n].j,
                                      mb_elements[n+2].i,mb_elements[n+2].j,
                                      GetStructure(),pfdata)));
    }
#endif//disablecoax
  return PROD(pfunc,initiation);
}


//functions for dealing with structures
hairpin_t hairpin(double p,int i, int j)
{
  hairpin_t h;
  h.probability = p;
  h.i = i;
  h.j = j;
  return h;
}

internal_loop_t internal_loop(double p,int i, int j, int k, int l)
{
  internal_loop_t il;
  il.i=i;
  il.j=j;
  il.k=k;
  il.l=l;
  il.probability=p;
  return il;
}

//to build a multibranch loop, call this with the closing base pair
//then call add_branch with the branches, from 5' to 3'
multibranch_loop_t multibranch_loop(int i, int j)
{
  multibranch_loop_t mb;
  mb.branches.push_back(std::make_pair(i,j));
  return mb;
}

void add_branch(multibranch_loop_t& mb,int k, int l)
{
  mb.branches.push_back(std::make_pair(k,l));
}

void show_hairpins(vector<hairpin_t> hairpins, bool fixed)//print hairpin output
{
  cout <<"--hairpins--"<<endl;
  cout << "prob i j" <<endl;
  for(vector<hairpin_t>::const_reverse_iterator it=hairpins.rbegin();it!=hairpins.rend();++it)
    if (fixed) cout << std::fixed<<std::setprecision(3)<<it->probability << " " << it->i << " " << it->j <<endl;
	else cout << std::scientific << std::setprecision(3) << it->probability << " " << it->i << " " << it->j << endl;
  cout<< "--hairpins end--"<<endl <<endl;
}

void show_internal_loops(vector<internal_loop_t> internals,bool fixed)//print iloop output
{
  cout << "--internal loops--"<<endl;
  cout << "prob i j k l"<<endl;
  for(vector<internal_loop_t>::const_reverse_iterator it=internals.rbegin();it!=internals.rend();++it)
    if (fixed) cout << std::fixed<<std::setprecision(3)<< it->probability << " " << it->i << " " << it->j <<" " << it->k << " " << it->l <<endl;
	else cout << std::scientific << std::setprecision(3) << it->probability << " " << it->i << " " << it->j << " " << it->k << " " << it->l << endl;
  cout<< "--internal loops end--"<<endl <<endl;
}

void show_bulge_loops(vector<internal_loop_t> internals, bool fixed)//print iloop output
{
  cout << "--bulge loops--"<<endl;
  cout << "prob i j k l"<<endl;
  for(vector<internal_loop_t>::const_reverse_iterator it=internals.rbegin();it!=internals.rend();++it)
    if (fixed) cout << std::fixed<<std::setprecision(3)<< it->probability << " " << it->i << " " << it->j <<" " << it->k << " " << it->l <<endl;
	else cout << std::scientific << std::setprecision(3) << it->probability << " " << it->i << " " << it->j << " " << it->k << " " << it->l << endl;
  cout<< "--bulge loops end--"<<endl <<endl;
}

void show_stacks(vector<basestack_t> stacks, bool fixed)//print stack output
{
  cout <<"--stacks--"<<endl;
  cout << "prob i j k l" <<endl;
  for(vector<basestack_t>::const_reverse_iterator it=stacks.rbegin();it!=stacks.rend();++it)
    if (fixed) cout << std::fixed<<std::setprecision(3)<< it->probability << " " << it->i << " " << it->j <<" " << it->k << " " << it->l <<endl;
	else cout << std::scientific << std::setprecision(3) << it->probability << " " << it->i << " " << it->j << " " << it->k << " " << it->l << endl;
  cout<< "--stacks end--"<<endl <<endl;
}
void show_mb_element_array(vector<mb_element> e)//just for debugging
{
  int count = 0;
  for(vector<mb_element>::iterator it=e.begin();it!=e.end();++it){
    cout<<count++<<" ";
    cout << (it->is_a_pair?"Pair: ":"Nuc ") << it->i<<" ";
    if(it->is_a_pair) cout <<it->j;
    cout <<std::endl;
  }
}

void show_mbl(multibranch_loop_t mb,bool fixed){//print multibranch loop output
	if (fixed) cout << std::fixed << std::setprecision(3) <<mb.probability;
	else cout << std::scientific << std::setprecision(3) << mb.probability;
  for(vector<std::pair<int,int> >::iterator it = mb.branches.begin();it!=mb.branches.end();++it){
    cout<<'\t'<<it->first<<'-'<<it->second;
  }
  cout<<'\n';
}
