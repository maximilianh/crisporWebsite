#ifndef LOG_DOUBLE_H
#define LOG_DOUBLE_H

#include "phmm/utils/xmath/log/xlog_math.h"
// #include <cmath>
// #include <cstdio>
// #include <cstdlib>
#include <iostream> // for insertion operator ( >> )

#ifdef PF_LOG_CLASS_STRICT
#define DBL_CONVERSION explicit
#else
#define DBL_CONVERSION
#endif

// This class is used to store values in the log scale while treating them as if they were in linear scale.
class log_double {
    public:
        double _value; // the number stored in log-scale 

        // Constructor and destructor
        log_double() {}
        // Create log_double from linear-scale int. Currently required for DynProgArray, but may be fixed (removed) in future.
        explicit log_double(const int& linear_value): _value(xlog(linear_value)){}
        // Copy constructor
        //log_double(const log_double& copy): _value(copy._value){}
        // An overload constructor that directly sets the inner value to the first argument. 
        // The first argument is not used, but should be set to true.
        log_double(const bool direct_init, const double& log_value): _value(log_value){}
        // Destructor not necessary
        //~log_double() {}

        // Convert between double and log_double when called
        DBL_CONVERSION operator double() const { return xexp(_value); }
        DBL_CONVERSION operator float() const { return (float)xexp(_value); }
        DBL_CONVERSION operator long double() const { return (long double)xexp(_value); }

        // Assignment from an existing log_double. (not needed. default assignment is fine)
        // log_double& operator=(const log_double& rhs) { _value = rhs._value; return *this; }

        // Assignment Operators for math operations between two log_double values.
        log_double& operator +=(const log_double& other) { _value = xlog_sum(_value, other._value);  return *this; }
        log_double& operator -=(const log_double& other) { _value = xlog_sub(_value, other._value);  return *this; }
        log_double& operator *=(const log_double& other) { _value = xlog_mul(_value, other._value);  return *this; }
        //log_double& operator *=(const log_double& other) { _value = _value + other._value;  return *this; } // faster, but skips LOG_OF_ZERO-checks
		log_double& operator /=(const log_double& other) { _value = xlog_div(_value, other._value);  return *this; }
        //log_double & operator /=(const log_double& other) { _value = _value - other._value;  return *this; } // faster, but skips LOG_OF_ZERO-checks

        #ifndef PF_LOG_CLASS_STRICT // disallow implicit conversions and operations with doubles in strict mode.
        // Constructor from linear-scale double, which sets _value to log(arg).
        DBL_CONVERSION log_double(const double& linear_value): _value(xlog(linear_value)){}

        // Assignment from a linear-scale double.
        log_double& operator=(const double& linear_value) { _value = xlog(linear_value); return *this; }

        // Binary operators with linear-scale doubles.
        log_double& operator +=(const double& linear_rhs) { _value = xlog_sum(_value, xlog(linear_rhs));  return *this; }
        log_double& operator -=(const double& linear_rhs) { _value = xlog_sub(_value, xlog(linear_rhs));  return *this; }
        log_double& operator *=(const double& linear_rhs) { _value = xlog_mul(_value, xlog(linear_rhs));  return *this; }
        log_double& operator /=(const double& linear_rhs) { _value = xlog_div(_value, xlog(linear_rhs));  return *this; }

        bool operator >(const double& other) const { return _value > xlog(other); }
        bool operator >=(const double& other) const { return _value >= xlog(other); }
        bool operator <(const double& other) const { return _value < xlog(other); }
        bool operator <=(const double& other) const { return _value <= xlog(other); }
        bool operator ==(const double& other) const { return _value == xlog(other); }
        #endif // !PF_LOG_CLASS_STRICT
		
        // Unary negation (-X)
    	//log_double operator-() const {return log_double(xexp(xlog_pow(_value, -1.0))); }

        bool operator ==(const log_double& other) const { return _value == other._value; }
        bool operator >(const log_double& other) const { return _value > other._value; }
        bool operator >=(const log_double& other) const { return _value >= other._value; }
        bool operator <(const log_double& other) const { return _value < other._value; }
        bool operator <=(const log_double& other) const { return _value <= other._value; }

        // log_double operator^(const int& rhs) const;

        // Initialize new log_double with a log-scale value.
        static inline log_double from_log(const double& log_value) {
            return log_double(true, log_value);
        }

        // Same as copy constructor, but can be used for polymorphism, e.g.:
        // 
        // a = log_double::from_log(X)  //  X might be `#define`d as log_double(3) ..or.. xlog(3)
        static inline log_double from_log(const log_double& copy) {
            return log_double(copy);
        }

        // static log_double cast(const double& linear_value) { 
        //     log_double out;
        //     out._value = xlog(linear_value);
        //     return out;
        // }
        // static log_double* cast(double *d) { 
        //     return reinterpret_cast<log_double*>(d);
        // }

        // Get the internal log-scale value.
        const double& log() const { return _value; }
        
        // The following static const members must be implemented in a linked object file (i.e. a cpp)
        static log_double zero() { return log_double(true, LOG_OF_ZERO); }
        static log_double one() { return log_double(true, LOG_OF_ONE); }
};

inline log_double ld_from_xlog(const double& xlog_value) {
    return log_double::from_log(xlog_value);
}
inline log_double ld_from_linear(const double& linear_value) { 
    #ifdef PF_LOG_CLASS_STRICT
    return log_double(true, xlog(linear_value));
    #else
    return log_double(linear_value);
    #endif
}

#ifndef PF_LOG_CLASS_STRICT // disallow implicit conversions and operations with doubles in strict mode.
// Binary operations with lhs log_double and rhs (linear scale) double.
inline log_double operator+(const double& linear_lhs, const log_double& rhs) { return ld_from_xlog(xlog_sum(xlog(linear_lhs), rhs._value));}
inline log_double operator-(const double& linear_lhs, const log_double& rhs) { return ld_from_xlog(xlog_sub(xlog(linear_lhs), rhs._value));}
inline log_double operator*(const double& linear_lhs, const log_double& rhs) { return ld_from_xlog(xlog_mul(xlog(linear_lhs), rhs._value));}
inline log_double operator/(const double& linear_lhs, const log_double& rhs) { return ld_from_xlog(xlog_div(xlog(linear_lhs), rhs._value));}
// Binary operations with lhs (linear scale) double and rhs log_double.
inline log_double operator+(const log_double& lhs, const double& linear_rhs) { return ld_from_xlog(xlog_sum(lhs._value,xlog(linear_rhs))); }
inline log_double operator-(const log_double& lhs, const double& linear_rhs) { return ld_from_xlog(xlog_sub(lhs._value,xlog(linear_rhs))); }
inline log_double operator*(const log_double& lhs, const double& linear_rhs) { return ld_from_xlog(xlog_mul(lhs._value,xlog(linear_rhs))); }
inline log_double operator/(const log_double& lhs, const double& linear_rhs) { return ld_from_xlog(xlog_div(lhs._value,xlog(linear_rhs))); }
#endif // !PF_LOG_CLASS_STRICT

// Binary operations between two log_double instances
inline log_double operator+(const log_double& lhs, const log_double& rhs) { return ld_from_xlog(xlog_sum(lhs._value,rhs._value)); }
inline log_double operator-(const log_double& lhs, const log_double& rhs) { return ld_from_xlog(xlog_sub(lhs._value,rhs._value)); }
inline log_double operator*(const log_double& lhs, const log_double& rhs) { return ld_from_xlog(xlog_mul(lhs._value,rhs._value)); }
inline log_double operator/(const log_double& lhs, const log_double& rhs) { return ld_from_xlog(xlog_div(lhs._value,rhs._value)); }

// Take the log of a log_double to return its inner (log-scale) value.
inline const double& log(const log_double& val){ return val._value; }
inline const double& xlog(const log_double& val){ return val._value; }

// ostream insertion operator
inline std::ostream& operator<<(std::ostream& os, const log_double& out)  { return os << xexp(out._value); }

inline log_double pow(const log_double& val, const double& exponent)     { return ld_from_xlog(xlog_pow(val._value, exponent)); }
inline log_double pow(const log_double& val, const log_double& exponent) { return ld_from_xlog(xlog_pow(val._value, xexp(exponent._value))); }
inline log_double pow(const log_double& val, const int exponent)        { return ld_from_xlog(xlog_pow(val._value, exponent)); }
inline log_double exp(const log_double& val) { return ld_from_xlog(xexp(val._value)); }

#ifndef __CUDA_CALC_
//typedef log_double real_t;
#endif

#endif // LOG_DOUBLE_H