#ifndef _EXTRINSIC_H_
#define _EXTRINSIC_H_
#include "alignment.h"
#include "arrays.h"
//#include "logdouble.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"
#include <vector>

typedef std::vector<std::vector<double>> table_t;

double proclivity(const alignment& aln, const table_t& probs);
double similarity(const sequence& s, const sequence& t);

template<typename T>
class extrinsic {
    public:
	    extrinsic();
        extrinsic(const sequence& s,
                  const std::vector<const sequence*> seqs,
                  const std::vector<const table_t*> probs,
                  const std::vector<const alignment*> alignments,
				  const T gamma=T(0.3));
        T bonus(const int i, const int j) const;
		void show() const;

    private:
		const bool initialized;
		const T gamma;
        std::vector<std::vector<T> > bonuses;
		T max_value() const;
		void elementwise_divide(const T denom);
		void elementwise_pow(const T denom);
		void normalize();
		void apply_gamma();
};

#endif
