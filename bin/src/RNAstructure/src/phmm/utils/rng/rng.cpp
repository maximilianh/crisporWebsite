#include <stdio.h>
#include <stdlib.h>
#include "rng.h"

#include <math.h>
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

t_rng::t_rng(int _prn_gen_seed)
{
	this->iff = 0; // not be modifed; see Knuth.
	this->prn_gen_seed = _prn_gen_seed;
}

t_rng::~t_rng()
{};

double t_rng::random_double_ran3()
{
	long idum = this->prn_gen_seed;
	double random_num = (double)this->ran3(&idum);
	//printf("Generated %f\n",random_num);
	return(random_num);
}

vector<int>* t_rng::permute_indices(int i_max, int i_to_exclude_from_first_position)
{
	vector<int>* permuted_indices = new vector<int>();
	vector<int>* index_bag = new vector<int>();

	// Fill the remaining indices.
	for(int cur_i = 0; cur_i < i_max; cur_i++)
	{
		// do not add the index to exclude to the bag.
		if(cur_i != i_to_exclude_from_first_position)
		{
			index_bag->push_back(cur_i);
		}
	}

	// generate the permutation.
	//for(int cur_i = 0; cur_i < index_bag->size(); cur_i++)
	int cur_i = 0;
	while(!index_bag->empty())
	{
		// Choose an index randomly from remaining indices.
		//int rand_i_bag = rand() % index_bag->size();

		// Following remainder operation is necessary for making sure that the 
		// generated random index is in range of index_bag. Note that the casting
		// takes ceiling of the floating point value rand_num, which is in range
		// [0, index_bag->size].
		double rand_num = index_bag->size() * this->random_double_ran3();
		int rand_i_bag = (int)rand_num % index_bag->size();

		// Copy the selected index to permuted indices.
		permuted_indices->push_back(index_bag->at(rand_i_bag));

		// Remove the element at selected position from remaining indices.
		index_bag->erase(index_bag->begin() + rand_i_bag);

		// After choosing first indices from index_bag, add the index to exclude from first
		// position to the bag so that it can be chosen again.
		if(cur_i == 0 && i_to_exclude_from_first_position >= 0 && i_to_exclude_from_first_position < i_max)
		{
			index_bag->push_back(i_to_exclude_from_first_position);
		}

		cur_i++;
	}

	if(index_bag->size() != 0)
	{
		printf("There are indices in the bag..\n");
		exit(0);
	}

	delete(index_bag);

	// Dump the permuted indices.
	for(int cur_i = 0; cur_i < i_max; cur_i++)
	{
		printf("%d ", permuted_indices->at(cur_i));
	}

	return(permuted_indices);
}


// idum number is the seed that is used to "warm up" the generator by 
// setting ma array. idum is set to 1 after that, i.e., it is not used any more.
// ma array is important and should not be changed while generating a sequence of random numbers.
// After initialization of ma array is utilized to generate random numbers.
float t_rng::ran3(long* idum)
{
	long mj,mk;
	int i,ii,k;

	// Initialization.
	if (*idum < 0 || iff == 0) 
	{ 
		iff=1;
		mj=labs(MSEED-labs(*idum)); // Initialize ma[55] using the seed idum and the
		mj %= MBIG; // large number MSEED.
		ma[55]=mj;
		mk=1;

		// Now initialize the rest of the table,
		for (i=1;i<=54;i++) 
		{ 
			ii=(21*i) % 55; // in a slightly random order,
			ma[ii]=mk; // with numbers that are not especially random.
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}

		// We randomize them by warming up the generator
		for (k=1;k<=4;k++) 
		{
			for(i=1;i<=55;i++) 
			{
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) 
				{
					ma[i] += MBIG;
				}
			}
		}

		inext=0; 
		inextp=31; // The constant 31 is special; see Knuth.
		*idum=1;
	}

	if (++inext == 56) 
	{
		inext=1; // Increment inext and inextp, wrapping around
	}

	if (++inextp == 56) 
	{
		inextp=1; 
	}

	// Generate a new random number subtractively.
	mj=ma[inext]-ma[inextp]; 

	// Be sure that it is in range.
	if (mj < MZ) 
	{
		mj += MBIG; 
	}

	ma[inext]=mj; // Store it,
	return mj*FAC; // and output the derived uniform deviate.
}


