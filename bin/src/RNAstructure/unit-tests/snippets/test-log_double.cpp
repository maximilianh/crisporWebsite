#include <cmath>
#include <cstdlib>

#define USE_XLOG_ZERO
#define STDIO

#ifdef STDIO
    #include <iostream>
    #include <iomanip>
    #define IF_STDIO(...) __VA_ARGS__
#else
    #define IF_STDIO(...) 
#endif

inline const bool logeq(const double& a, const double& b) {
    return abs((a-b)/a) < 1E-4;
}

#ifdef D1
    #define PF_LOG_CLASS
    #define PF_LOG_CLASS_STRICT
    #define INNER(X) X.logval()
#else
    #define INNER(X) X
#endif
#include "../../src/pfunction_math.h"
#include "../../src/log_double.h"

#define DTYPE PFPRECISION
#define DISPLAY(X) "{ log=" << INNER(X) << ", linear=" << TO_LINEAR(X) << " }"

#ifndef COUNT
    #define COUNT -1
#endif

int main() {
    IF_STDIO(std::cout << "COUNT=" << COUNT << "\n");
    IF_STDIO(std::cout << std::setprecision(16));

    #define MATH_OP(P,Q) PROD(P,Q,P,Q,P,Q)
    #define MATH_RES(D) pow(D, 9)

    int errors = 0;

    //#if COUNT == 0
    double d1, r1, s1; DTYPE p1, q1, m1;
    d1 = rand()/100.0;
    p1 = TO_XLOG(d1);
    q1 = PROD(p1, p1);
    m1 = MATH_OP(p1,q1);
    r1 = TO_LINEAR(p1);
    s1 = TO_LINEAR(m1);
    double expected = MATH_RES(d1);
    IF_STDIO(std::cout << -1 << "." << 
        "\td=" << d1 << "\tr=" << r1 << "\ts=" << s1 <<
        "\n\tp=" << DISPLAY(p1) << 
        "\n\tq=" << DISPLAY(q1) <<
        "\n\tm=" << DISPLAY(m1) << 
        "\n");
    if (!logeq(d1,r1) || !logeq(s1,expected)) {
        IF_STDIO(std::cout << "!!!BAD " << -1 << "\n");
        errors++;
    }    //#endif //COUNT==0

    #if COUNT > 0
        double d[COUNT], r[COUNT], s[COUNT];
        DTYPE  p[COUNT], q[COUNT], m[COUNT];
        const int count = COUNT;
    #elif COUNT < 0
        const int count = 1+10*rand()/RAND_MAX;
        IF_STDIO(std::cout << "dyn-count=" << count << "\n");
        double *d, *r, *s; 
        DTYPE *p, *q, *m;
        d = new double[count];
        r = new double[count];
        s = new double[count];
        p = new DTYPE[count];
        q = new DTYPE[count];
        m = new DTYPE[count];
    #endif

    #if COUNT != 0
    for(int i=0; i<count;i++) {
        d[i] = rand()/100.0;
        p[i] = TO_XLOG(d[i]);
        q[i] = PROD(p[i], p[i]);
        m[i] = MATH_OP(p[i], q[i]);
        r[i] = TO_LINEAR(p[i]);
        s[i] = TO_LINEAR(m[i]);
        double expected = MATH_RES(d[i]);
        IF_STDIO(std::cout << i << "." << 
        "\td=" << d[i] << "\tr=" << r[i] << "\ts=" << s[i] <<
        "\n\tp=" << DISPLAY(p[i]) << 
        "\n\tq=" << DISPLAY(q[i]) <<
        "\n\tm=" << DISPLAY(m[i]) << 
        "\n");
        if (!logeq(d[i],r[i]) || !logeq(s[i],expected)) {
            IF_STDIO(std::cout << "!!!BAD " << i << "\n");
            errors++;
        }
    }
    #endif

    #if count < 0
        delete[] d, r, p, q;
    #endif

    // for(int i=0; i<count;i++)
    //     d[i] = rand()/100.0;

    // for(int i=0; i<count;i++) {
    //     p[i] = TO_XLOG(d[i]);
    //     q[i] = TO_XLOG(2*d[i]);
    // }

    // for(int i=0; i<count;i++)
    //     r[i] = TO_LINEAR(p[i]);

    // for(int i=0; i<count;i++) {
    //     IF_STDIO(std::cout << i << ".  " << d[i] << ", " << DISPLAY(p[i]) << ", " << DISPLAY(q[i]) << ", " << r[i] << "\n");
    //     if (!logeq(d[i],r[i])) {
    //         IF_STDIO(std::cout << "!!!BAD[" << i << "] " << d[i] << ", " << DISPLAY(p[i]) << ", " << DISPLAY(q[i]) << ", " << r[i] << "\n");
    //         errors++;
    //     }
    // }


    // PFPRECISION rsum = SUM(p[0], p[1], p[2]);
    // double dsum = TO_LINEAR(rsum);

    // PFPRECISION rprod = PROD(p[0], p[1], p[2]);
    // double dprod = TO_LINEAR(rprod);

    //std::cout << "PROD: " << dprod << " SUM:" << dsum << "\n";

    IF_STDIO(std::cout.flush());

    return errors;
}

/*  //other template test
#include <iostream>

#ifdef D1
inline double scale(const double& A, const double& B) {
    #ifdef MSG
    std::cout << "scaling " << A << " and " << B << std::endl;
    #endif
    return A;
}
#elif defined D2
inline double scale(const double& A) {
    #ifdef MSG
    std::cout << "scaling " << A << std::endl;
    #endif
    return A;
}
#else
#define scale(X) X
#endif

int main() {
    double X, Y, Z;
    std::cin >> X;
    std::cin >> Y;
    std::cin >> Z;

    #ifdef D1
    std::cout << scale(X, Y) << std::endl;
    std::cout << scale(Z, Y) << std::endl;
    #else
    std::cout << scale(X) << std::endl;
    std::cout << scale(Z) << std::endl;
    #endif

    return 0;
}
*/