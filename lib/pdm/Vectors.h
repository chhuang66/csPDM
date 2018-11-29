///////////////////////////////////////////////////////////////////////////////
// Vectors.h
// =========
// 2D
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2007-02-14
// UPDATED: 2013-01-20
//
// Copyright (C) 2007-2013 Song Ho Ahn
///////////////////////////////////////////////////////////////////////////////


#ifndef VECTORS_H_DEF
#define VECTORS_H_DEF

#include <cmath>
#include <iostream>
#include <stdexcept>

///////////////////////////////////////////////////////////////////////////////
// 2D vector
///////////////////////////////////////////////////////////////////////////////
struct Vector2
{
    double x;
    double y;

    // ctors
    Vector2() : x(0), y(0) {};
    Vector2(double x, double y) : x(x), y(y) {};

    // utils functions
    void        set(double x, double y);
    double       length() const;                         //
    double       distance(const Vector2& vec) const;     // distance between two vectors
    Vector2&    normalize();                            //
    double       dot(const Vector2& vec) const;          // dot product
    bool        equal(const Vector2& vec, double e) const; // compare with epsilon
    double ave() const;
    // operators
    Vector2     operator-() const;                      // unary operator (negate)
    Vector2     operator+(const Vector2& rhs) const;    // add rhs
    Vector2     operator-(const Vector2& rhs) const;    // subtract rhs
    Vector2&    operator+=(const Vector2& rhs);         // add rhs and update this object
    Vector2&    operator-=(const Vector2& rhs);         // subtract rhs and update this object
    Vector2     operator*(const double scale) const;     // scale
    Vector2     operator*(const Vector2& rhs) const;    // multiply each element
    Vector2&    operator*=(const double scale);          // scale and update this object
    Vector2&    operator*=(const Vector2& rhs);         // multiply each element and update this object
    Vector2     operator/(const double scale) const;     // inverse scale
    Vector2&    operator/=(const double scale);          // scale and update this object
    bool        operator==(const Vector2& rhs) const;   // exact compare, no epsilon
    bool        operator!=(const Vector2& rhs) const;   // exact compare, no epsilon
    bool        operator<(const Vector2& rhs) const;    // comparison for sort
    double       operator[](int index) const;            // subscript operator v[0], v[1]
    double&      operator[](int index);                  // subscript operator v[0], v[1]

    friend Vector2 operator*(const double a, const Vector2 vec);
    friend std::ostream& operator<<(std::ostream& os, const Vector2& vec);
};

inline double Vector2::ave() const
{
 return (x+y)/2.0;
}

// fast math routines from Doom3 SDK
inline double invSqrt(double x)
{
    double xhalf = 0.5f * x;
    int i = *(int*)&x;          // get bits for doubleing value
    i = 0x5f3759df - (i>>1);    // gives initial guess
    x = *(double*)&i;            // convert bits back to double
    x = x * (1.5f - xhalf*x*x); // Newton step
    return x;
}



///////////////////////////////////////////////////////////////////////////////
// inline functions for Vector2
///////////////////////////////////////////////////////////////////////////////
inline Vector2 Vector2::operator-() const {
    return Vector2(-x, -y);
}

inline Vector2 Vector2::operator+(const Vector2& rhs) const {
    return Vector2(x+rhs.x, y+rhs.y);
}

inline Vector2 Vector2::operator-(const Vector2& rhs) const {
    return Vector2(x-rhs.x, y-rhs.y);
}

inline Vector2& Vector2::operator+=(const Vector2& rhs) {
    x += rhs.x; y += rhs.y; return *this;
}

inline Vector2& Vector2::operator-=(const Vector2& rhs) {
    x -= rhs.x; y -= rhs.y; return *this;
}

inline Vector2 Vector2::operator*(const double a) const {
    return Vector2(x*a, y*a);
}

inline Vector2 Vector2::operator*(const Vector2& rhs) const {
    return Vector2(x*rhs.x, y*rhs.y);
}

inline Vector2& Vector2::operator*=(const double a) {
    x *= a; y *= a; return *this;
}

inline Vector2& Vector2::operator*=(const Vector2& rhs) {
    x *= rhs.x; y *= rhs.y; return *this;
}

inline Vector2 Vector2::operator/(const double a) const {
    return Vector2(x/a, y/a);
}

inline Vector2& Vector2::operator/=(const double a) {
    x /= a; y /= a; return *this;
}

inline bool Vector2::operator==(const Vector2& rhs) const {
    return (x == rhs.x) && (y == rhs.y);
}

inline bool Vector2::operator!=(const Vector2& rhs) const {
    return (x != rhs.x) || (y != rhs.y);
}

inline bool Vector2::operator<(const Vector2& rhs) const {
    if(x < rhs.x) return true;
    if(x > rhs.x) return false;
    if(y < rhs.y) return true;
    if(y > rhs.y) return false;
    return false;
}

inline double Vector2::operator[](int index) const {
	if(index < 0 || index >= 2)
    	throw std::out_of_range("blah");
    return (&x)[index];
}

inline double& Vector2::operator[](int index) {
    return (&x)[index];
}

inline void Vector2::set(double x, double y) {
    this->x = x; this->y = y;
}

inline double Vector2::length() const {
    return sqrtf(x*x + y*y);
}

inline double Vector2::distance(const Vector2& vec) const {
    return sqrtf((vec.x-x)*(vec.x-x) + (vec.y-y)*(vec.y-y));
}

inline Vector2& Vector2::normalize() {
    //@@const double EPSILON = 0.000001f;
    double xxyy = x*x + y*y;
    //@@if(xxyy < EPSILON)
    //@@    return *this;

    //double invLength = invSqrt(xxyy);
    double invLength = 1.0f / sqrtf(xxyy);
    x *= invLength;
    y *= invLength;
    return *this;
}

inline double Vector2::dot(const Vector2& rhs) const {
    return (x*rhs.x + y*rhs.y);
}

inline bool Vector2::equal(const Vector2& rhs, double epsilon) const {
    return fabs(x - rhs.x) < epsilon && fabs(y - rhs.y) < epsilon;
}

inline Vector2 operator*(const double a, const Vector2 vec) {
    return Vector2(a*vec.x, a*vec.y);
}

inline std::ostream& operator<<(std::ostream& os, const Vector2& vec) {
    os << "(" << vec.x << ", " << vec.y << ")";
    return os;
}
// END OF VECTOR2 /////////////////////////////////////////////////////////////

#endif
