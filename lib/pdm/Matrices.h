///////////////////////////////////////////////////////////////////////////////
// Matrice.h
// =========
// NxN Matrix Math classes
//
// The elements of the matrix are stored as column major order.
// | 0 2 |    | 0 3 6 |    |  0  4  8 12 |
// | 1 3 |    | 1 4 7 |    |  1  5  9 13 |
//            | 2 5 8 |    |  2  6 10 14 |
//                         |  3  7 11 15 |
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2005-06-24
// UPDATED: 2013-09-30
//
// Copyright (C) 2005 Song Ho Ahn
///////////////////////////////////////////////////////////////////////////////

#ifndef MATH_MATRICES_H
#define MATH_MATRICES_H

#include <iostream>
#include <iomanip>
#include "Vectors.h"

///////////////////////////////////////////////////////////////////////////
// 2x2 matrix
///////////////////////////////////////////////////////////////////////////

class Matrix2
{
public:
    // constructors
    Matrix2();  // init with identity
    Matrix2(const double src[4]);
    Matrix2(double m0, double m1, double m2, double m3);

    void        set(const double src[4]);
    void        set(double m0, double m1, double m2, double m3);
    void        setRow(int index, const double row[2]);
    void        setRow(int index, const Vector2& v);
    void        setColumn(int index, const double col[2]);
    void        setColumn(int index, const Vector2& v);
    double ave() const;
    const double* get() const;
    double       getDeterminant();
    
    Matrix2&    identity();
    Matrix2&    transpose();                            // transpose itself and return reference
    Matrix2&    invert();

    // operators
    Matrix2     operator+(const Matrix2& rhs) const;    // add rhs
    Matrix2     operator-(const Matrix2& rhs) const;    // subtract rhs
    Matrix2&    operator+=(const Matrix2& rhs);         // add rhs and update this object
    Matrix2&    operator-=(const Matrix2& rhs);         // subtract rhs and update this object
    Vector2     operator*(const Vector2& rhs) const;    // multiplication: v' = M * v
    Matrix2     operator*(const Matrix2& rhs) const;    // multiplication: M3 = M1 * M2
    Matrix2&    operator*=(const Matrix2& rhs);         // multiplication: M1' = M1 * M2
    bool        operator==(const Matrix2& rhs) const;   // exact compare, no epsilon
    bool        operator!=(const Matrix2& rhs) const;   // exact compare, no epsilon
    double       operator[](int index) const;            // subscript operator v[0], v[1]
    double&      operator[](int index);                  // subscript operator v[0], v[1]

    friend Matrix2 operator-(const Matrix2& m);                     // unary operator (-)
    friend Matrix2 operator*(double scalar, const Matrix2& m);       // pre-multiplication
    friend Vector2 operator*(const Vector2& vec, const Matrix2& m); // pre-multiplication
    friend std::ostream& operator<<(std::ostream& os, const Matrix2& m);

protected:

private:
    double m[4];

};


///////////////////////////////////////////////////////////////////////////
// inline functions for Matrix2
///////////////////////////////////////////////////////////////////////////
inline double Matrix2::ave() const
{
 return (m[0]+m[1]+m[2]+m[3])/4.0;
}

inline Matrix2::Matrix2()
{
    // initially identity matrix
    identity();
}



inline Matrix2::Matrix2(const double src[4])
{
    set(src);
}



inline Matrix2::Matrix2(double m0, double m1, double m2, double m3)
{
    set(m0, m1, m2, m3);
}



inline void Matrix2::set(const double src[4])
{
    m[0] = src[0];  m[1] = src[1];  m[2] = src[2];  m[3] = src[3];
}



inline void Matrix2::set(double m0, double m1, double m2, double m3)
{
    m[0]= m0;  m[1] = m1;  m[2] = m2;  m[3]= m3;
}



inline void Matrix2::setRow(int index, const double row[2])
{
    m[index] = row[0];  m[index + 2] = row[1];
}



inline void Matrix2::setRow(int index, const Vector2& v)
{
    m[index] = v.x;  m[index + 2] = v.y;
}



inline void Matrix2::setColumn(int index, const double col[2])
{
    m[index*2] = col[0];  m[index*2 + 1] = col[1];
}



inline void Matrix2::setColumn(int index, const Vector2& v)
{
    m[index*2] = v.x;  m[index*2 + 1] = v.y;
}



inline const double* Matrix2::get() const
{
    return m;
}



inline Matrix2& Matrix2::identity()
{
    m[0] = m[3] = 1.0f;
    m[1] = m[2] = 0.0f;
    return *this;
}



inline Matrix2 Matrix2::operator+(const Matrix2& rhs) const
{
    return Matrix2(m[0]+rhs[0], m[1]+rhs[1], m[2]+rhs[2], m[3]+rhs[3]);
}



inline Matrix2 Matrix2::operator-(const Matrix2& rhs) const
{
    return Matrix2(m[0]-rhs[0], m[1]-rhs[1], m[2]-rhs[2], m[3]-rhs[3]);
}



inline Matrix2& Matrix2::operator+=(const Matrix2& rhs)
{
    m[0] += rhs[0];  m[1] += rhs[1];  m[2] += rhs[2];  m[3] += rhs[3];
    return *this;
}



inline Matrix2& Matrix2::operator-=(const Matrix2& rhs)
{
    m[0] -= rhs[0];  m[1] -= rhs[1];  m[2] -= rhs[2];  m[3] -= rhs[3];
    return *this;
}



inline Vector2 Matrix2::operator*(const Vector2& rhs) const
{
    return Vector2(m[0]*rhs.x + m[2]*rhs.y,  m[1]*rhs.x + m[3]*rhs.y);
}



inline Matrix2 Matrix2::operator*(const Matrix2& rhs) const
{
    return Matrix2(m[0]*rhs[0] + m[2]*rhs[1],  m[1]*rhs[0] + m[3]*rhs[1],
                   m[0]*rhs[2] + m[2]*rhs[3],  m[1]*rhs[2] + m[3]*rhs[3]);
}



inline Matrix2& Matrix2::operator*=(const Matrix2& rhs)
{
    *this = *this * rhs;
    return *this;
}



inline bool Matrix2::operator==(const Matrix2& rhs) const
{
    return (m[0] == rhs[0]) && (m[1] == rhs[1]) && (m[2] == rhs[2]) && (m[3] == rhs[3]);
}



inline bool Matrix2::operator!=(const Matrix2& rhs) const
{
    return (m[0] != rhs[0]) || (m[1] != rhs[1]) || (m[2] != rhs[2]) || (m[3] != rhs[3]);
}



inline double Matrix2::operator[](int index) const
{
    return m[index];
}



inline double& Matrix2::operator[](int index)
{
    return m[index];
}



inline Matrix2 operator-(const Matrix2& rhs)
{
    return Matrix2(-rhs[0], -rhs[1], -rhs[2], -rhs[3]);
}



inline Matrix2 operator*(double s, const Matrix2& rhs)
{
    return Matrix2(s*rhs[0], s*rhs[1], s*rhs[2], s*rhs[3]);
}



inline Vector2 operator*(const Vector2& v, const Matrix2& rhs)
{
    return Vector2(v.x*rhs[0] + v.y*rhs[1],  v.x*rhs[2] + v.y*rhs[3]);
}



inline std::ostream& operator<<(std::ostream& os, const Matrix2& m)
{
    os << std::fixed << std::setprecision(5);
    os << "[" << std::setw(10) << m[0] << " " << std::setw(10) << m[2] << "]\n"
       << "[" << std::setw(10) << m[1] << " " << std::setw(10) << m[3] << "]\n";
    os << std::resetiosflags(std::ios_base::fixed | std::ios_base::floatfield);
    return os;
}
// END OF MATRIX2 INLINE //////////////////////////////////////////////////////


#endif
