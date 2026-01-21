#include <iostream>
#include <type_traits>
#include <cstring>
#include <string>
#include <vector>

// binary_fold provides a means to aggregate/reduce a list of arguments over a binary
// operation/function. The fold can be a left- or right-associative fold.
// `Operation` must be a callable type, e.g. a functor type -- a class/struct with operator() defined.
//
// The advantage of this approach is that the operation can be templated, so the arguments and 
// return value can be all of different types.
//
// Examples of `Operation` include op_add and op_multiply.
template<typename Operation>
struct binary_fold {
    // Fold operations are divided into two parts (for clarity):
    //   1) Return Type: result_type<bool,Args...>::type (aka result_t) determines the return type of the entire fold operation.
    //   2) Functions: lfold<Args...>() and  rfold<Args...>() perform the actual fold calculation.

    // op_t<A,B> is the return type of calling `Operation` with argument types A and B.
    template<typename A, typename B> using op_t = typename std::result_of<Operation(A,B)>::type;  // alternative: decltype(std::declval<Operation>()(std::declval<A>(), std::declval<B>()));
    // result_type: determines the return-type of a fold operation.
    // This is the default (invalid) case with 0 type args.
    template<bool LeftFold, typename ... Args> struct result_type {};
    // Alias for result_type<...>::type;
    template<bool LeftFold, typename ... Args> using result_t = typename result_type<LeftFold, Args...>::type;
    // terminal case with 1 arg -- just return the type of the arg.
    template<bool LeftFold, typename T> struct result_type<LeftFold,T> { using type=T; };
    // general case -- handles left fold (due to subsequent specialization)
    template<typename T, typename U, typename...Args> 
    struct result_type<true,T,U,Args...> { using type = result_t<true,op_t<T,U>, Args...>; };  // for left fold:  (A,B,C,D) --> (((A,B),C),D) type is op_t<op_t<op_t<A,B>,C>,D>
    // partial specialization handles right fold
    template<typename T, typename U, typename...Args> 
    struct result_type<false,T,U,Args...> { using type = op_t<T, result_t<false,U,Args...>>; }; // for right fold: (A,B,C,D) --> (A,(B,(C,D))) type is op_t<A, op_t<B, op_t<C,D>>>

    // left-fold terminal case (1 arg)
    template<typename T> static result_t<true,T> lfold(const T& t) { return t; }
    // left-fold general case
    template<typename T, typename U, typename...Args>
    static result_t<true,T,U,Args...> lfold(const T& t, const U& u, const Args&... args) { return lfold<op_t<T,U>,Args...>(Operation()(t, u), args...); }// left fold

    // right-fold terminal case (1 arg)
    template<typename T> static result_t<false,T> rfold(const T& t) { return t; }
    // right-fold general case
    template<typename T, typename...Args>
    static result_t<false,T,Args...> rfold(const T& t, const Args&... args) { return Operation()(t, rfold<Args...>(args...)); }
};
// Alias for binary_fold<Operation>::result_type<bool,...>::type
template<typename Operation, bool LeftFold, typename... Args>
using fold_result_t = typename binary_fold<Operation>::template result_t<LeftFold,Args...>;

// Left Fold
template<typename Operation, typename... Args> 
fold_result_t<Operation,true,Args...> fold(Args...args) { return binary_fold<Operation>::lfold(args...); }
// Right Fold
template<typename Operation, typename... Args> 
fold_result_t<Operation,false,Args...> rfold(Args...args) { return binary_fold<Operation>::rfold(args...); }

// Provides the addition operation.
struct op_add {
    template<typename T, typename U>
    inline auto operator()(const T& lhs, const U& rhs) -> decltype(lhs+rhs) {
        return lhs+rhs;
    }
};

// Provides the multiplication operation.
struct op_multiply {
    template<typename T, typename U>
    inline auto operator()(const T& lhs, const U& rhs) -> decltype(lhs*rhs) {
        return lhs*rhs;
    }
};

// reduces a list of objects by addition 
template<typename...Ts>
fold_result_t<op_add,true,Ts...> fold_add(const Ts&...ts) { return fold<op_add,Ts...>(ts...); }
// reduces a list of objects by multiplication 
template<typename...Ts>
fold_result_t<op_multiply,true,Ts...> fold_mul(const Ts&...ts) { return fold<op_multiply,Ts...>(ts...); }


/* This section contains fold_left and fold_right templates that
   accept a binary function of the form `T func(const T& a, const T& b)`
   This is a simpler approach, but the downside is that the type T cannot
   be inferred (because binary_op depends on it and is not itself a type, 
   but a non-type parameter), so T needs to be specified explicitly.
*/
template <typename T> using binary_op = T(*)(const T& a, const T& b);
// The terminal (single-argument) cases of the variadic functions defined later. 
template<typename T, binary_op<T> Operation> inline T bin_fold_right(const T& t) { return t; }
template<typename T, binary_op<T> Operation> inline T bin_fold_left(const T& t) { return t; }

// Combines arguments using right-associative operation
// i.e. bin_fold_right<T,op>(A,B,C) --> op(A, op(B,C))
template<typename T, binary_op<T> Operation, typename ... Rest> 
inline T bin_fold_right(const T& t, Rest... rest) {
	return Operation(t, bin_fold_right<T, Operation>(rest...));
}
// Combines arguments using left-associative operation
// i.e. bin_fold_left<T,op>(A,B,C) --> op(op(A,B), C)
template<typename T, binary_op<T> Operation, typename ... Rest>
inline T bin_fold_left(const T& a, const T& b, Rest... rest) {
    return bin_fold_left<T, Operation>(Operation(a,b), rest...);
}
template<typename T> inline T add_op(const T& a, const T& b) { return a+b; }
template<typename T> inline T mul_op(const T& a, const T& b) { return a*b; }
template<typename T> inline T div_op(const T& a, const T& b) { return a/b; }



// Classes for tests of asymmetric reductions
struct Base { const virtual char* const out() const = 0; };
struct A:Base { const char* const out() const { return "A"; } };
struct B:Base { const char* const out() const { return "B"; } };
struct C:Base { const char* const out() const { return "C"; } };
C operator+(A a, B b) { return C(); } // A+B=C
C operator+(B a, A b) { return C(); } // B+A=C
B operator+(A a, C b) { return B(); } // A+C=B
B operator+(C a, A b) { return B(); } // C+A=B
A operator+(B a, C b) { return A(); } // B+C=A
A operator+(C a, B b) { return A(); } // C+B=A
A operator+(A a, A b) { return A(); } // A+A=A
B operator+(B a, B b) { return B(); } // B+B=B
C operator+(C a, C b) { return C(); } // C+C=C

#define CONCAT_(x,y) x##y
#define CONCAT(x,y) CONCAT_(x,y)
#define TOSTRING_(x) #x
#define TOSTRING(x) TOSTRING_(x)
#define TEST(T,E) TEST_(T,CONCAT(T,__LINE__),E)
#define TEST_(T,N,E) T N = (E); std::cout<< TOSTRING(N) << "=" << N.out() << std::endl;
#define TEST_FL(T,E) TEST_F_(T,CONCAT(T,__LINE__),E,fold) //fold_left
#define TEST_FR(T,E) TEST_F_(T,CONCAT(T,__LINE__),E,rfold) //fold_right
#define TEST_F_(T,N,E,FOLD) T N = FOLD<op_add>E; std::cout << TOSTRING(N) " " TOSTRING(FOLD) TOSTRING(E) " = " << N << " ex " TOSTRING(T) << std::endl;
#define DBG(S) std::cout << #S <<" = "<< S << "\n"

std::ostream& operator<<(std::ostream &out, const Base& obj) { 
    return out << obj.out();
}

int main() {
    A a; B b; C c;
    TEST(A,a+a);
    TEST(B,a+(a+b));
    TEST(C,a+((a+b)+a));
    TEST(B,b+((a+b)+a));
    TEST(B,b+((a+b)+a));
    TEST(C,a+b);
    TEST(B,(a+b)+a);
    TEST(A,(a+b)+b);
    TEST(C,(a+b)+c);

    TEST_FL(A,(a,a));
    TEST_FL(B,(b,b));
    TEST_FL(C,(c,c));

    TEST_FL(A,(a,a,a));
    TEST_FL(B,(b,b,b));
    TEST_FL(C,(c,c,c));

    TEST_FL(B,(c,c,a));
    TEST_FL(C,(a,b));
    TEST_FL(B,(a,b,a));
    TEST_FL(C,(a,a,b));
    TEST_FL(A,(a,b,b));
    TEST_FL(C,(a,b,c));
    TEST_FL(B,(c,b,c));
    TEST_FL(B,(a,a,b,a));
    TEST_FL(A,(b,a,b,a));
    TEST_FL(B,(b,c,b,a));
    TEST_FL(C,(b,c,b,c));
    
    TEST_FR(A,(a,a));
    TEST_FR(C,(a,b));
    TEST_FR(B,(a,b,a));
    TEST_FR(B,(a,a,b));
    TEST_FR(C,(a,b,b));
    TEST_FR(A,(a,b,c));
    TEST_FR(B,(c,b,c));
    TEST_FR(C,(a,a,b,a));
    TEST_FR(B,(b,a,b,a));
    TEST_FR(A,(b,c,b,a));
    TEST_FR(B,(b,c,b,c));

    fold_result_t<op_add, false, A> zR1 = A();
    fold_result_t<op_add, false, A, B> zR2 = C();
    fold_result_t<op_add, false, A, B, C> zR3 = A();
    fold_result_t<op_add, false, A, C, C, B> zR4 = C();
    fold_result_t<op_add, false, A, A, C, C, B> zR5 = B();

    fold_result_t<op_add, true, A> zL1 = A();
    fold_result_t<op_add, true, A, B> zL2 = C();
    fold_result_t<op_add, true, A, B, C> zL3 = C();
    fold_result_t<op_add, true, A, C, C, B> zL4 = C();
    fold_result_t<op_add, true, A, A, C, C, B> zL5 = C();

    std::cout << fold<op_add>(1,2,3,4) << std::endl;
    std::cout << fold<op_multiply>(1,2,3,4) << std::endl;

    std::cout << "fold_add="<<fold_add(1,2,3,4) << std::endl;
    std::cout << "fold_mul="<<fold_mul(1,2,3,4) << std::endl;

    std::cout << "fl_ta:" << fold<op_add>(a) << std::endl;
    std::cout << "fl_ab:" << fold<op_add>(a,b) << std::endl;
    std::cout << "fl_abc:" << fold_add(a,b,c) << std::endl;
    std::cout << "fr_abc:" << rfold<op_add>(a,b, c) << std::endl;
 
    return 0;
}