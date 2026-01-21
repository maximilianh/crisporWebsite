#ifndef LOOP_UTILS
#define LOOP_UTILS

#include "../RNA_class/RNA.h"
#include "../RNA_class/ProbScan.h"
#include <vector>

namespace loop{
    class basepair{
        public:
        int i;
        int j;
        basepair(int i, int j);
        bool operator==(const basepair& other) const;
        basepair & operator=(const basepair & other);
    };

    class loop{
        public:
        loop(int i, int j);
        basepair outer;
        virtual std::vector<int> nucs() const = 0;
        virtual double getProbability(ProbScan& p) const = 0;
    };

    class hairpin : public loop
    {
        public:
        hairpin(int i, int j);
        std::vector<int> nucs() const;
        double getProbability(ProbScan& p) const;
    };

    class internal : public loop
    {
        public:
        internal(int i, int j, int k, int l);
        basepair inner;
        std::vector<int> nucs() const;
        double getProbability(ProbScan& p) const;
    };

    class multibranch : public loop
    {
        public:
        multibranch(std::vector<basepair> pairs);
        std::vector<basepair> pairs;
        std::vector<int> nucs() const;
        double getProbability(ProbScan& p) const;
    };

    class stem : public loop
    {
        public:
        stem(int i, int j, int k, int l);
        basepair inner;
        std::vector<int> nucs() const;
        double getProbability(ProbScan& p) const;
    };
    std::ostream& operator<<(std::ostream& output, const hairpin& h);
    std::ostream& operator<<(std::ostream& output, const internal& h);
    std::ostream& operator<<(std::ostream& output, const multibranch& h);
    std::ostream& operator<<(std::ostream& output, const stem& h);

    std::vector<hairpin> find_hairpins(RNA& r,int structurenumber=1);
    std::vector<internal> find_internals(RNA& r,int structurenumber=1);
    std::vector<multibranch> find_multibranch(RNA& r,int structurenumber=1);
    std::vector<stem> find_stems(RNA& r,int structurenumber=1);

}
#endif
