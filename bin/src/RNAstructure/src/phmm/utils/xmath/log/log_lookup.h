#pragma once

#include <cmath>
#include <iostream>
#include <iomanip>
#include "../../../../defines.h"
//#define ENABLE_DEBUG_LOGS
#include "../../../../debug_logging.h"
#define LUPRECISION double
// #define LU_INTERP_LINEAR // enable linear interpolation in the lookup table.
#define LUINTERP_CUBIC // enable quartic interpolation.

// Provides a fast lookup of the function `log1p(exp(x))`, which provides a more 
// accurate result for  `log(1+exp(x))` 
// This is used to optimize the calculation of the sum of two log-scale values.
//      log(exp(a)+exp(b)) = log( exp(a) * (1+exp(b-a)) ) = a + log1p(exp(b-a))
class log_lookup_sum{
    private:
        double delta, delta2, delta3, low;
        int steps;
        LUPRECISION *lookup_table;

#if defined LUINTERP_CUBIC
        LUPRECISION *a, *b, *c;
#endif
    public:
        // low_in:   The lowest negative number to be evaluated using lookup tables
        // steps_in: The number of table entries
        log_lookup_sum(LUPRECISION low_in, int steps_in){
            low = low_in;
            steps = steps_in;

            lookup_table = new LUPRECISION[steps+3];
            delta = (-low)/steps;
            delta2 = pow((double) delta,2);
            delta3 = pow((double) delta,3);

            for(int i=0;i<steps+3;i++){
		        double x = (1-i)*delta;
                lookup_table[i] = log1p(exp(x));
            }

#if defined LUINTERP_CUBIC
            a = new LUPRECISION[steps+3];
            b = new LUPRECISION[steps+3];
            c = new LUPRECISION[steps+3];
            
            for(int i=0;i<steps+3;i++){
                if (i>0 && i<steps+1){
                    a[i] = delta2 * (2*lookup_table[i-1] + 3*lookup_table[i] - 6*lookup_table[i+1] + lookup_table[i+2])/(6*delta3);
                    b[i] = 3*delta/(6*delta3)*(lookup_table[i-1] - 2*lookup_table[i] + lookup_table[i+1]);
                    c[i] = (lookup_table[i-1] - 3*lookup_table[i] + 3*lookup_table[i+1] - lookup_table[i+2])/(6*delta3);
                }
                else{
                    a[i]=0;
                    b[i]=0;
                    c[i]=0;
                }
            }
#endif

        }
        ~log_lookup_sum(){
            delete[] lookup_table;
#if defined LUINTERP_CUBIC
			delete[] a;
			delete[] b;
			delete[] c;
#endif
        }

        LUPRECISION exp_lu(LUPRECISION x){
            if (x < low) return 0;
            else {
                
#if defined LU_INTERP_LINEAR
            // Linear interpolation
            LUPRECISION f = -x/delta;
            int index = (int)f;
            f -= index; // get fractional part [0 to 1)
            return lookup_table[index+2]*f+lookup_table[index+1]*(1-f);
#elif defined LUINTERP_CUBIC
            x=-x;
            //quartic interpolation
            int index = (int) (x/delta) + 1;
            LUPRECISION mu = x - (index-1)*delta;

//            LUPRECISION a0,a1,a2,a3,mu2;
			LUPRECISION mu2;

            mu2 = mu*mu;

            return  (lookup_table[index]
                    - a[index]*mu 
                    + b[index]*mu2 
                    - c[index]*mu2*mu);
#else 
            // No interpolation
            // LOG_DEBUG4("[log_lookup]\tindex: " << index)
            return lookup_table[1+(int)std::round(-x/delta)];
#endif
            }
        //    }
        }

        void test(){
	    LUPRECISION x;
        int sampling = 100;
            for(int i=0;i<-low/delta*sampling;i++){
        		x = exp_lu(-delta*(i-1)/sampling);
                // std::cout << std::setprecision(16) << x << "\t" << exp_lu(x) << "\t" << log1p(exp(x)) << std::endl;
            }
        }
};

// Provides a fast lookup of the function `log1p(-exp(-x))`, which provides a more 
// accurate result for  `log(1-exp(-x))` 
// This is used to optimize the calculation of the difference between two log-scale values.
//      log(exp(a)-exp(b)) = log( exp(a) * (1-exp(b-a)) ) = a + log1p(-exp(-(b-a)))
class log_lookup_sub{
    private:
        LUPRECISION delta, delta2, delta3, high;
        int steps;
        LUPRECISION *lookup_table;

    public:
        log_lookup_sub(LUPRECISION high_in, int steps_in){
            high= high_in;
            steps = steps_in;

            lookup_table = new LUPRECISION[steps+2];
            delta = high/steps;
            delta2 = pow(delta,2);
            delta3 = pow(delta,3);
            for(int i=1;i<steps+2;i++){
                double x = i*delta;
                lookup_table[i-1] = log1p(-exp(-x));
            }
        }
        ~log_lookup_sub(){
            delete[] lookup_table;
        }

        LUPRECISION exp_lu(LUPRECISION x){
            //std::cout << x << std::endl;
            if (x < 3*delta) return log1p(-exp(-x));
            else{
                if (x > high-2*delta) return 0;
                else{
                    int index = (int) (x/delta-0.5);
                    return lookup_table[index];

                    // int index = (int) (x/delta) + 1;
                    // LUPRECISION x0 = (index-1)*delta;
                    // return lookup_table[index]+(x-x0)*(lookup_table[index+1]-lookup_table[index])/delta;
                    
                            //     LUPRECISION x0 = (index+1)*delta;
            
                //     return -((-6*lookup_table[index]*delta3 + delta2 * (2*lookup_table[index-1] + 3*lookup_table[index] - 6*lookup_table[index+1] + lookup_table[index+2])*(x - x0) 
                //   - 3*delta*(lookup_table[index-1] - 2*lookup_table[index] + lookup_table[index+1])*(x - x0)*(x - x0) 
                //   + (lookup_table[index-1] - 3*lookup_table[index] + 3*lookup_table[index+1] - lookup_table[index+2])*(x - x0)*(x - x0)*(x - x0))/(6*delta3));
                }
            }
        }

        void test(){
            int steps = (int) high/delta;
	        LUPRECISION x;
            for(int i=1;i<steps;i++){
		        x = -delta*i;
                std::cout << std::setprecision(16) << x << "\t" << exp_lu(-x) << "\t" << log1p(-exp(x)) << std::endl;
            }
        }
        
};//*/