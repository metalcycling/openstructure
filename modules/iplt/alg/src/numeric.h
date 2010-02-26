// Copyright (c) 2007, 2008 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// Matrix and vector classes, based on Eigen2.
//
// Avoid using Eigen2 classes directly; instead typedef them here.

#ifndef LIBMV_NUMERIC_NUMERIC_H
#define LIBMV_NUMERIC_NUMERIC_H

#include <Eigen/Array>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/SVD>

#include <ost/base.hh>

#include <iostream>

#if _WIN32 || __APPLE__
  void static sincos (Real x, Real *sinx, Real *cosx) {
    *sinx = sin(x);
    *cosx = cos(x);
  }
#endif //_WIN32 || __APPLE__

// #if _WIN32
  // inline long lround(double d) {
    // return (long)(d>0 ? d+0.5 : ceil(d-0.5));
  // }
  // inline int round(double d) {
    // return (d>0) ? int(d+0.5) : int(d-0.5);
  // }
  // typedef unsigned int uint;
// #endif //_WIN32

namespace ost { namespace iplt { namespace alg {

typedef Eigen::MatrixXf Mat;
typedef Eigen::VectorXf Vec;

typedef Eigen::MatrixXf Matf;
typedef Eigen::VectorXf Vecf;

typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> Matu;
typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> Vecu;

typedef Eigen::Matrix<Real, 2, 2> Mat2;
typedef Eigen::Matrix<Real, 2, 3> Mat23;
typedef Eigen::Matrix<Real, 3, 3> Mat3;
typedef Eigen::Matrix<Real, 3, 4> Mat34;
typedef Eigen::Matrix<Real, 3, 5> Mat35;
typedef Eigen::Matrix<Real, 4, 3> Mat43;
typedef Eigen::Matrix<Real, 4, 4> Mat4;
typedef Eigen::Matrix<Real, 4, 6> Mat46;
typedef Eigen::Matrix<float, 2, 2> Mat2f;
typedef Eigen::Matrix<float, 2, 3> Mat23f;
typedef Eigen::Matrix<float, 3, 3> Mat3f;
typedef Eigen::Matrix<float, 3, 4> Mat34f;
typedef Eigen::Matrix<float, 3, 5> Mat35f;
typedef Eigen::Matrix<float, 4, 3> Mat43f;
typedef Eigen::Matrix<float, 4, 4> Mat4f;
typedef Eigen::Matrix<float, 4, 6> Mat46f;

typedef Eigen::Matrix<Real, 3, 3, Eigen::RowMajor> RMat3;

typedef Eigen::Matrix<Real, 2, Eigen::Dynamic> Mat2X;
typedef Eigen::Matrix<Real, 3, Eigen::Dynamic> Mat3X;
typedef Eigen::Matrix<Real, 4, Eigen::Dynamic> Mat4X;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 2> MatX2;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 3> MatX3;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 4> MatX4;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 9> MatX9;

typedef Eigen::Vector2f Vec2;
typedef Eigen::Vector3f Vec3;
typedef Eigen::Vector4f Vec4;
typedef Eigen::Matrix<Real, 5, 1> Vec5;
typedef Eigen::Matrix<Real, 6, 1> Vec6;
typedef Eigen::Matrix<Real, 7, 1> Vec7;
typedef Eigen::Matrix<Real, 8, 1> Vec8;
typedef Eigen::Matrix<Real, 9, 1> Vec9;
typedef Eigen::Matrix<Real, 10, 1> Vec10;
typedef Eigen::Matrix<Real, 13, 1> Vec13;

typedef Eigen::Vector2f Vec2f;
typedef Eigen::Vector3f Vec3f;
typedef Eigen::Vector4f Vec4f;

typedef Eigen::VectorXi VecXi;

typedef Eigen::Vector2i Vec2i;
typedef Eigen::Vector3i Vec3i;
typedef Eigen::Vector4i Vec4i;

typedef Eigen::Matrix<float,
                      Eigen::Dynamic,
                      Eigen::Dynamic,
                      Eigen::RowMajor> RMatf;

using Eigen::Map;
using Eigen::Dynamic;
using Eigen::Matrix;

// Find U, s, and VT such that
//
//   A = U * diag(s) * VT
//
template <typename TMat, typename TVec>
inline void SVD(TMat *A, Vec *s, Mat *U, Mat *VT) {
  assert(0);
}

// Solve the linear system Ax = 0 via SVD. Store the solution in x, such that
// ||x|| = 1.0. Return the singluar value corresponding to the solution.
// Destroys A and resizes x if necessary.
// TODO(maclean): Take the SVD of the transpose instead of this zero padding.
template <typename TMat, typename TVec>
Real Nullspace(TMat *A, TVec *nullspace) {
  if (A->rows() >= A->cols()) {
    Eigen::SVD<TMat> svd(*A);
    (*nullspace) = svd.matrixV().col(A->cols()-1);
    return svd.singularValues()(A->cols()-1);
  }
  // Extend A with rows of zeros to make it square. It's a hack, but is
  // necessary until Eigen supports SVD with more columns than rows.
  Mat A_extended(A->cols(), A->cols());
  A_extended.block(A->rows(), 0, A->cols() - A->rows(), A->cols()).setZero();
  A_extended.block(0,0, A->rows(), A->cols()) = (*A);
  return Nullspace(&A_extended, nullspace);
}

// Solve the linear system Ax = 0 via SVD. Finds two solutions, x1 and x2, such
// that x1 is the best solution and x2 is the next best solution (in the L2
// norm sense). Store the solution in x1 and x2, such that ||x|| = 1.0. Return
// the singluar value corresponding to the solution x1.  Destroys A and resizes
// x if necessary.
template <typename TMat, typename TVec1, typename TVec2>
Real Nullspace2(TMat *A, TVec1 *x1, TVec2 *x2) {
  if (A->rows() >= A->cols()) {
    Eigen::SVD<TMat> svd(*A);
    Mat V = svd.matrixV();
    *x1 = V.col(A->cols() - 1);
    *x2 = V.col(A->cols() - 2);
    return svd.singularValues()(A->cols()-1);
  }
  // Extend A with rows of zeros to make it square. It's a hack, but is
  // necessary until Eigen supports SVD with more columns than rows.
  Mat A_extended(A->cols(), A->cols());
  A_extended.block(A->rows(), 0, A->cols() - A->rows(), A->cols()).setZero();
  A_extended.block(0,0, A->rows(), A->cols()) = (*A);
  return Nullspace2(&A_extended, x1, x2);
}

// In place transpose for square matrices.
template<class TA>
inline void TransposeInPlace(TA *A) {
  *A = A->transpose().eval();
}

template<typename TVec>
inline Real NormL1(const TVec &x) {
  return x.cwise().abs().sum();
}

template<typename TVec>
inline Real NormL2(const TVec &x) {
  return x.norm();
}

template<typename TVec>
inline Real NormLInfinity(const TVec &x) {
  return x.cwise().abs().maxCoeff();
}

template<typename TVec>
inline Real DistanceL1(const TVec &x, const TVec &y) {
  return (x - y).cwise().abs().sum();
}

template<typename TVec>
inline Real DistanceL2(const TVec &x, const TVec &y) {
  return (x - y).norm();
}
template<typename TVec>
inline Real DistanceLInfinity(const TVec &x, const TVec &y) {
  return (x - y).cwise().abs().maxCoeff();
}

// Normalize a vector with the L1 norm, and return the norm before it was
// normalized.
template<typename TVec>
inline Real NormalizeL1(TVec *x) {
  Real norm = NormL1(*x);
  *x /= norm;
  return norm;
}

// Normalize a vector with the L2 norm, and return the norm before it was
// normalized.
template<typename TVec>
inline Real NormalizeL2(TVec *x) {
  Real norm = NormL2(*x);
  *x /= norm;
  return norm;
}

// Normalize a vector with the L^Infinity norm, and return the norm before it
// was normalized.
template<typename TVec>
inline Real NormalizeLInfinity(TVec *x) {
  Real norm = NormLInfinity(*x);
  *x /= norm;
  return norm;
}

// Return the square of a number.
template<typename T>
inline T Square(T x) {
  return x * x;
}

Mat3 RotationAroundX(Real angle);
Mat3 RotationAroundY(Real angle);
Mat3 RotationAroundZ(Real angle);

// Returns the rotation matrix of a rotation of angle |axis| around axis.
// This is computed using the Rodrigues formula, see:
//  http://mathworld.wolfram.com/RodriguesRotationFormula.html
Mat3 RotationRodrigues(const Vec3 &axis);

// Make a rotation matrix such that center becomes the direction of the
// positive z-axis, and y is oriented close to up.
Mat3 LookAt(Vec3 center);

// Return a diagonal matrix from a vector containg the diagonal values.
template <typename TVec>
inline Mat Diag(const TVec &x) {
  return x.asDiagonal();
}

template<typename TMat>
inline Real FrobeniusNorm(const TMat &A) {
  return sqrt(A.cwise().abs2().sum());
}

template<typename TMat>
inline Real FrobeniusDistance(const TMat &A, const TMat &B) {
  return FrobeniusNorm(A - B);
}

inline Vec3 CrossProduct(const Vec3 &x, const Vec3 &y) {
  return x.cross(y);
}

Mat3 CrossProductMatrix(const Vec3 &x);

void MeanAndVarianceAlongRows(const Mat &A,
                              Vec *mean_pointer,
                              Vec *variance_pointer);

#if _WIN32
  // TODO(bomboze): un-#if this for both platforms once tested under Windows
  /* This solution was extensively discussed here http://forum.kde.org/viewtopic.php?f=74&t=61940 */
  #define SUM_OR_DYNAMIC(x,y) (x==Eigen::Dynamic||y==Eigen::Dynamic)?Eigen::Dynamic:(x+y)

  template<typename Derived1, typename Derived2>
  struct hstack_return {
    typedef typename Derived1::Scalar Scalar;
    enum {
         RowsAtCompileTime = Derived1::RowsAtCompileTime,
         ColsAtCompileTime = SUM_OR_DYNAMIC(Derived1::ColsAtCompileTime, Derived2::ColsAtCompileTime),
         Options = Derived1::Flags&Eigen::RowMajorBit ? Eigen::RowMajor : 0,
         MaxRowsAtCompileTime = Derived1::MaxRowsAtCompileTime,
         MaxColsAtCompileTime = SUM_OR_DYNAMIC(Derived1::MaxColsAtCompileTime, Derived2::MaxColsAtCompileTime)
    };
    typedef Eigen::Matrix<Scalar,
                RowsAtCompileTime,
                ColsAtCompileTime,
                Options,
                MaxRowsAtCompileTime,
                MaxColsAtCompileTime> type;
  };

  template<typename Derived1, typename Derived2>
  typename hstack_return<Derived1,Derived2>::type
  HStack (const Eigen::MatrixBase<Derived1>& lhs, const Eigen::MatrixBase<Derived2>& rhs) {
    typename hstack_return<Derived1,Derived2>::type res;
    res.resize(lhs.rows(), lhs.cols()+rhs.cols());
    res << lhs, rhs;
    return res;
  };


  template<typename Derived1, typename Derived2>
  struct vstack_return {
    typedef typename Derived1::Scalar Scalar;
    enum {
         RowsAtCompileTime = SUM_OR_DYNAMIC(Derived1::RowsAtCompileTime, Derived2::RowsAtCompileTime),
         ColsAtCompileTime = Derived1::ColsAtCompileTime,
         Options = Derived1::Flags&Eigen::RowMajorBit ? Eigen::RowMajor : 0,
         MaxRowsAtCompileTime = SUM_OR_DYNAMIC(Derived1::MaxRowsAtCompileTime, Derived2::MaxRowsAtCompileTime),
         MaxColsAtCompileTime = Derived1::MaxColsAtCompileTime
    };
    typedef Eigen::Matrix<Scalar,
                RowsAtCompileTime,
                ColsAtCompileTime,
                Options,
                MaxRowsAtCompileTime,
                MaxColsAtCompileTime> type;
  };

  template<typename Derived1, typename Derived2>
  typename vstack_return<Derived1,Derived2>::type
  VStack (const Eigen::MatrixBase<Derived1>& lhs, const Eigen::MatrixBase<Derived2>& rhs) {
    typename vstack_return<Derived1,Derived2>::type res;
    res.resize(lhs.rows()+rhs.rows(), lhs.cols());
    res << lhs, rhs;
    return res;
  };


#else //_WIN32

  // Since it is not possible to typedef privately here, use a macro.
  // Always take dynamic columns if either side is dynamic.
  #define COLS \
    ((ColsLeft == Eigen::Dynamic || ColsRight == Eigen::Dynamic) \
     ? Eigen::Dynamic : (ColsLeft + ColsRight))

  // Same as above, except that prefer fixed size if either is fixed.
  #define ROWS \
    ((RowsLeft == Eigen::Dynamic && RowsRight == Eigen::Dynamic) \
     ? Eigen::Dynamic \
     : ((RowsLeft == Eigen::Dynamic) \
        ? RowsRight \
        : RowsLeft \
       ) \
    )

  // TODO(keir): Add a static assert if both rows are at compiletime.
  template<typename T, int RowsLeft, int RowsRight, int ColsLeft, int ColsRight>
  Eigen::Matrix<T, ROWS, COLS>
  HStack(const Eigen::Matrix<T, RowsLeft,  ColsLeft>  &left,
         const Eigen::Matrix<T, RowsRight, ColsRight> &right) {
    assert(left.rows() == right.rows());
    int n = left.rows();
    int m1 = left.cols();
    int m2 = right.cols();

    Eigen::Matrix<T, ROWS, COLS> stacked(n, m1 + m2);
    stacked.block(0, 0,  n, m1) = left;
    stacked.block(0, m1, n, m2) = right;
    return stacked;
  }

  // Reuse the above macros by swapping the order of Rows and Cols. Nasty, but
  // the duplication is worse.
  // TODO(keir): Add a static assert if both rows are at compiletime.
  // TODO(keir): Mail eigen list about making this work for general expressions
  // rather than only matrix types.
  template<typename T, int RowsLeft, int RowsRight, int ColsLeft, int ColsRight>
  Eigen::Matrix<T, COLS, ROWS>
  VStack(const Eigen::Matrix<T, ColsLeft,  RowsLeft>  &top,
         const Eigen::Matrix<T, ColsRight, RowsRight> &bottom) {
    assert(top.cols() == bottom.cols());
     int n1 = top.rows();
    int n2 = bottom.rows();
    int m = top.cols();

    Eigen::Matrix<T, COLS, ROWS> stacked(n1 + n2, m);
    stacked.block(0,  0, n1, m) = top;
    stacked.block(n1, 0, n2, m) = bottom;
    return stacked;
  }
  #undef COLS
  #undef ROWS
#endif //_WIN32



void HorizontalStack(const Mat &left, const Mat &right, Mat *stacked);

template<typename TTop, typename TBot, typename TStacked>
void VerticalStack(const TTop &top, const TBot &bottom, TStacked *stacked) {
  assert(top.cols() == bottom.cols());
  int n1 = top.rows();
  int n2 = bottom.rows();
  int m = top.cols();

  stacked->resize(n1 + n2, m);
  stacked->block(0,  0, n1, m) = top;
  stacked->block(n1, 0, n2, m) = bottom;
}

void MatrixColumn(const Mat &A, int i, Vec2 *v);
void MatrixColumn(const Mat &A, int i, Vec3 *v);
void MatrixColumn(const Mat &A, int i, Vec4 *v);

template <typename TMat, typename TCols>
TMat ExtractColumns(const TMat &A, const TCols &columns) {
  TMat compressed(A.rows(), columns.size());
  for (int i = 0; i < columns.size(); ++i) {
    compressed.col(i) = A.col(columns[i]);
  }
  return compressed;
}

template <typename TMat, typename TDest>
void reshape(const TMat &a, int rows, int cols, TDest *b) {
  assert(a.rows()*a.cols() == rows*cols);
  b->resize(rows, cols);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      (*b)(i, j) = a[cols*i + j];
    }
  }
}

}}}  // namespace mv

#endif  // LIBMV_NUMERIC_NUMERIC_H
