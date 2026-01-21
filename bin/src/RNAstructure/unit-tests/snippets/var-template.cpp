#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include <typeinfo>

template<typename T> 
void f(const T t) { std::cout << "T0: " << t << std::endl; }

template<typename T> 
void f(const T t, std::nullptr_t p) { f(t); }

template<typename T, typename ...Args> 
void f(const T t, const Args... args) { 
    std::cout << "T" << (sizeof...(args)) << " (" <<  typeid(typeof(T)).name() <<"): " <<t<< ", "; 
    f(args...);
}

template<typename T, typename ...Args> 
void f(const T t, std::nullptr_t p, const Args... args) { 
    std::cout << "T" << (sizeof...(args)) << " (" <<  typeid(typeof(T)).name() <<"): " <<t<< ", "; 
    f(args...);
}

template<typename ...Args> 
void f(const std::nullptr_t p, const Args... args) { 
    f(args...);
}

// template<typename T> 
// void f(T t, const unsigned int i) { std::cout << "T: " << t << " int i==" << i << std::endl; }

// template<typename T> 
// void f(T t, const char* const c) { std::cout << "T: " << t << " char c==" << c << std::endl; }

int main() {
    #ifdef D0
    f("banana", 1.025);
    f("acorn",45);
    const std::string john = "John Wallace";
    f(1,"frogs",2.023,true,john,"pizza","face");
    f(nullptr,2,3);
    f(1,nullptr,3);
    f(1,7,8,nullptr,2,3);
    f(1,2,nullptr);
    #else
        #ifdef D1
            #define RESULT ,nullptr
        #else
            #define RESULT 
        #endif
        f(1 RESULT);
    #endif

    return 0;
}