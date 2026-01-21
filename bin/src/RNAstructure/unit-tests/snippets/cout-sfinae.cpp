#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <cxxabi.h>
#include <boost/type_index.hpp> // optional. if missing, use typeid in show_info

/*
This demonstrates how to use SFINAE to write an ostream insertion operator (<<) 
that will only be used for some types, depending on whether:
    A). The type has specific member (e.g. a function or Functor, etc) defined
        for itself or in a base class.  (e.g. `obj.serialize()` in this example)
 or
    B). A specific function has been defined and/or overloaded 
        for that type (e.g. `to_string` in this example).

This technique can be extended to enable/disable any template 
function for specific classes.
*/

// class A has a to_string overload defined for it, but no A.serialize()
struct A {};
std::string to_string(const A&) { return "I am to_string(A)"; }

// class B has a serialize method.
struct B {
    std::string serialize() const { return "I am B.serialize()!"; }
};

// class C has a serialize that returns an int (which we allow in this example)
struct C {
    int serialize() const {
        return 42;
    }
};

// type D inherits from A, so it has a to_string. But it also implements serialize()
struct D : A {
    std::string serialize() const {
        return "I am D.serialize()!";
    }
};

// type E has a field called serialize, which is a Functor, so it is callable via E().serialize()
struct E {
    struct Functor {
        std::string operator()() const{ 
            return "I am E's Functor!";
        }
    };
    Functor serialize;
};

// Type F has a "bad" serialize member --i.e. it is not callable. 
// But F does have a to_string overload.
struct F {
    std::string serialize;
};
std::string to_string(const F&) { return "I am to_string(F)"; }

// type G is just a double, which does not have serialize or to_string
// but it already has an ostream<< overload.
typedef double G;

struct H{}; // no serialize or to_string or ostream<< overload -- cannot be inserted with <<

template <class T> struct has_serialize {
    // We test if the type has serialize using decltype and declval.
    template <typename C> static constexpr decltype(std::declval<C>().serialize(), bool()) test(int) {
        // We can return values, thanks to constexpr instead of playing with sizeof.
        return true;
    }
    template <typename C> static constexpr bool test(...) {
        return false;
    }
    // int is used to give the precedence!
    static constexpr bool value = test<T>(int());
};

template <class T> struct has_to_string {
    // We test if the type has serialize using decltype and declval.
    template <typename C> static constexpr decltype(to_string(std::declval<C>()), bool()) test(int) {
        // We can return values, thanks to constexpr instead of playing with sizeof.
        return true;
    }
    template <typename C> static constexpr bool test(...) {
        return false;
    }
    // int is used to give the precedence!
    static constexpr bool value = test<T>(int());
};

template<typename T, typename TrueType=T>
using can_serialize_t = typename std::enable_if<has_serialize<T>::value||has_to_string<T>::value, TrueType>::type;

template<typename T>
can_serialize_t<T,std::ostream&> operator<<(std::ostream& out, const T& obj) { return out << serialize(obj); }


template<typename T>
auto serialize(T& obj) -> typename std::enable_if<has_serialize<T>::value, decltype(obj.serialize())>::type {
    return obj.serialize();
}

template<typename T> 
auto serialize(T& obj) -> typename std::enable_if<has_to_string<T>::value&&!has_serialize<T>::value, decltype(to_string(obj))>::type {
    return to_string(obj);
}

#define DBG(X) std::cout<< #X << " = " << X << "\n"

template<typename T>
void show_info(const T& t) {
    std::ostream& out = std::cout;
    auto name = boost::typeindex::type_id<T>().pretty_name(); // alternative: typeid(T).name()
    out << "Type: ";
    out.width(6); out.setf(std::ios_base::left);  // pad-right to 6 chars for type name.
    out << name << "\tValue: ";
    out.width(20); out.setf(std::ios_base::left); // pad-right to 20 chars for serialized value.
    
    out << t;     // uses serialize or to_string or ostream<< overload
    
    out << "\tInfo: ";
    if (has_to_string<T>::value && has_serialize<T>::value)
        out << "Has both to_string and serialize.";
    else if (has_to_string<T>::value)
        out << "Has to_string.";
    else if (has_serialize<T>::value)
        out << "Has serialize.";
    else
        out << "Has neither.";
    out << std::endl;
}


int main() {
    A a;
    B b;
    C c;
    D d;
    E e;
    F f;
    G g; // double 
    H h; // empty struct

    show_info(a);
    show_info(b);
    show_info(c);
    show_info(d);
    show_info(e);
    show_info(f);
    show_info(g); // double -- does not use serialize or to_string
    // show_info(h); // error. no serialize or to_string or ostream<< overload

    return 0;
}