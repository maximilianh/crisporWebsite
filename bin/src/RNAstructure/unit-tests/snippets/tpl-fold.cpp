#include <iostream>

// Binary function that accepts two arguments of the same type (T) and returns the same type.
template <typename T> using binary_op = T(*)(const T& a, const T& b);

template<typename T, binary_op<T> Operation> inline T fold_right_assoc(T t) { return t; }
template<typename T, binary_op<T> Operation> inline T fold_left_assoc(T t) { return t; }

template<typename T, binary_op<T> Operation, typename ... Rest> 
inline T fold_right_assoc(T t, Rest... rest) {
	return Operation(t, fold_right_assoc<T, Operation>(rest...));
}
template<typename T, binary_op<T> Operation, typename ... Rest> 
inline T fold_left_assoc(Rest... rest, T t) {
	return Operation(fold_left_assoc<T, Operation>(rest...), t);
}

inline int add(const int& A, const int& B) { return A+B; }

#if defined D1
    //#define add_all(...) fold_left_assoc<int,add>(__VA_ARGS__)
    #define add_all(A,B,C,D)  add(add(add(A,B),C),D)
#elif defined D2
    #define add_all(A,B,C,D)  add(add(add(A,B),C),D)
#elif defined D3
    #define add_all(A,B,C,D)  A+B+C+D
#elif defined D4
    #define add_all(...) fold_right_assoc<int,add>(__VA_ARGS__)
#elif defined D5
    #define add_all(A,B,C,D)  add(A, add(B, add(C,D))) // add(add(add(A,B),C),D)
#elif defined D6
    #define add_all(A,B,C,D)  A+(B+(C+D))
#else 
    #error No test defined!
#endif

int main() {
    int i=rand();
    int j=rand();
    int k=rand();
    int l=rand();
    int ans = add_all(i,j,k,l);
    std::cout << "Ans=" << ans << std::endl;
    return 0;
}