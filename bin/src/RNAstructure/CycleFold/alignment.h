#include <vector>
#include "sequence.h"
#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

class alignment {
    public:
        alignment(const sequence& a, const sequence& b);
        double coincidence(const int i, const int k) const;
        double similarity() const;
        void show() const;

    private:
        std::vector<std::vector<double> > probs;
        double _similarity;
};
#endif
