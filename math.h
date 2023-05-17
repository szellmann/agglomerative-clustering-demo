// ======================================================================== //
// Copyright 2022-2023 Stefan Zellmann                                      //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

#include <cmath>
#include <cstddef>
#include <ostream>

namespace math {
struct vec2f
{
  vec2f() = default;
  vec2f(float s) : x(s), y(s) {}
  vec2f(float x, float y) : x(x), y(y) {}
  float &operator[](int i) { return ((float*)this)[i]; }
  const float &operator[](int i) const { return ((float*)this)[i]; }
  float x, y;
};

inline
vec2f operator+(vec2f u, vec2f v) {
  return {u.x+v.x,u.y+v.y};
}

inline
vec2f operator-(vec2f u, vec2f v) {
  return {u.x-v.x,u.y-v.y};
}

inline
vec2f operator*(vec2f u, vec2f v) {
  return {u.x*v.x,u.y*v.y};
}

inline
vec2f operator/(vec2f u, vec2f v) {
  return {u.x/v.x,u.y/v.y};
}

inline
vec2f min(vec2f u, vec2f v) {
  return {fminf(u.x,v.x),fminf(u.y,v.y)};
}

inline
vec2f max(vec2f u, vec2f v) {
  return {fmaxf(u.x,v.x),fmaxf(u.y,v.y)};
}

inline
float dot(vec2f u, vec2f v) {
  return u.x*v.x+u.y*v.y;
}

inline
float norm2(vec2f u) {
  return dot(u,u);
}

inline
float length(vec2f u) {
  return sqrtf(dot(u,u));
}

inline
std::ostream& operator<<(std::ostream &out, vec2f v) {
  out << '(' << v.x << ',' << v.y << ')';
  return out;
}

struct vec3f
{
  vec3f() = default;
  vec3f(float s) : x(s), y(s), z(s) {}
  vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
  float &operator[](int i) { return ((float*)this)[i]; }
  const float &operator[](int i) const { return ((float*)this)[i]; }
  float x, y, z;
};

inline
vec3f operator+(vec3f v, float a) {
  return {v.x+a,v.y+a,v.z+a};
}

inline
vec3f operator+(vec3f u, vec3f v) {
  return {u.x+v.x,u.y+v.y,u.z+v.z};
}

inline
vec3f operator-(vec3f u, vec3f v) {
  return {u.x-v.x,u.y-v.y,u.z-v.z};
}

inline
vec3f operator*(vec3f u, vec3f v) {
  return {u.x*v.x,u.y*v.y,u.z*v.z};
}

inline
vec3f operator/(vec3f u, vec3f v) {
  return {u.x/v.x,u.y/v.y,u.z/v.z};
}

inline
vec3f min(vec3f u, vec3f v) {
  return {fminf(u.x,v.x),fminf(u.y,v.y),fminf(u.z,v.z)}; }

inline
vec3f max(vec3f u, vec3f v) {
  return {fmaxf(u.x,v.x),fmaxf(u.y,v.y),fmaxf(u.z,v.z)};
}

inline
float reduce_min(vec3f u) {
  return fminf(fminf(u.x,u.y),u.z);
}

inline
float reduce_max(vec3f u) {
  return fmaxf(fmaxf(u.x,u.y),u.z);
}

inline
float dot(vec3f u, vec3f v) {
  return u.x*v.x+u.y*v.y+u.z*v.z;
}

inline
vec3f normalize(vec3f u) {
  return u / sqrtf(dot(u,u));
}

inline
std::ostream& operator<<(std::ostream &out, vec3f v) {
  out << '(' << v.x << ',' << v.y <<',' << v.z << ')';
  return out;
}

struct vec2i
{
  vec2i() = default;
  vec2i(int s) : x(s), y(s) {}
  vec2i(int x, int y) : x(x), y(y) {}
  int &operator[](int i) { return ((int*)this)[i]; }
  const int &operator[](int i) const { return ((int*)this)[i]; }
  int x, y;
};

inline
vec2i operator+(vec2i u, vec2i v) {
  return {u.x+v.x,u.y+v.y};
}

inline
vec2i operator-(vec2i u, vec2i v) {
  return {u.x-v.x,u.y-v.y};
}

inline
vec2i operator*(vec2i u, vec2i v) {
  return {u.x*v.x,u.y*v.y};
}

inline
vec2i operator/(vec2i u, vec2i v) {
  return {u.x/v.x,u.y/v.y};
}

inline
std::ostream& operator<<(std::ostream &out, vec2i v) {
  out << '(' << v.x << ',' << v.y << ')';
  return out;
}

struct vec3i
{
  vec3i() = default;
  vec3i(int s) : x(s), y(s), z(s) {}
  vec3i(int x, int y, int z) : x(x), y(y), z(z) {}
  int &operator[](int i) { return ((int*)this)[i]; }
  const int &operator[](int i) const { return ((int*)this)[i]; }
  int x, y, z;
};

inline
std::ostream& operator<<(std::ostream &out, vec3i v) {
  out << '(' << v.x << ',' << v.y <<',' << v.z << ')';
  return out;
}

struct vec2ui
{
  vec2ui() = default;
  vec2ui(int s) : x(s), y(s) {}
  vec2ui(unsigned x, unsigned y) : x(x), y(y) {}
  unsigned &operator[](int i) { return ((unsigned*)this)[i]; }
  const unsigned &operator[](int i) const { return ((unsigned*)this)[i]; }
  unsigned x, y;
};

inline
vec2ui operator-(vec2ui u, vec2ui v) {
  return {u.x-v.x,u.y-v.y};
}

inline
std::ostream& operator<<(std::ostream &out, vec2ui v) {
  out << '(' << v.x << ',' << v.y << ')';
  return out;
}

inline
bool operator==(vec2ui u, vec2ui v)
{
  return u.x==v.x && u.y==v.y;
}

struct box1f
{
  box1f() = default;
  box1f(float lo, float up) : lower(lo), upper(up) {}

  inline
  bool empty() const {
    return upper <= lower;
  }

  inline
  float center() const {
    return (lower+upper)/2;
  }

  inline
  void extend(float v) {
    lower = fminf(lower,v);
    upper = fmaxf(upper,v);
  }

  float lower, upper;
};

inline
std::ostream& operator<<(std::ostream &out, box1f b) {
  out << '(' << b.lower << ',' << b.upper << ')';
  return out;
}

struct box2f
{
  box2f() = default;
  box2f(vec2f lo, vec2f up) : lower(lo), upper(up) {}

  inline
  bool empty() const {
    return upper.x <= lower.x || upper.y <= lower.y;
  }

  inline
  vec2f center() const {
    return (lower+upper)/2;
  }

  inline
  vec2f size() const {
    return upper-lower;
  }

  inline
  bool contains(vec2f p) const {
    return lower.x<=p.x && p.x<=upper.x
        && lower.y<=p.y && p.y<=upper.y;
  }

  inline
  void extend(vec2f v) {
    lower = min(lower,v);
    upper = max(upper,v);
  }

  inline
  void extend(box2f other) {
    extend(other.lower);
    extend(other.upper);
  }

  vec2f lower, upper;
};

inline
float area(box2f b) {
  vec2f v = b.upper-b.lower;
  return v.x*v.y;
}

inline
std::ostream& operator<<(std::ostream &out, box2f b) {
  out << '(' << b.lower << ',' << b.upper << ')';
  return out;
}

struct  box3f
{
  box3f() = default;
  box3f(vec3f lo, vec3f up) : lower(lo), upper(up) {}

  inline
  bool empty() const {
    return upper.x <= lower.x || upper.y <= lower.y || upper.z <= lower.z;
  }

  inline
  vec3f center() const {
    return (lower+upper)/2;
  }

  vec3f lower, upper;
};

inline
std::ostream& operator<<(std::ostream &out, box3f b) {
  out << '(' << b.lower << ',' << b.upper << ')';
  return out;
}

struct box3i
{
  box3i() = default;
  box3i(vec3i lo, vec3i up) : lower(lo), upper(up) {}

  vec3i lower, upper;
};

inline
std::ostream& operator<<(std::ostream &out, box3i b) {
  out << '(' << b.lower << ',' << b.upper << ')';
  return out;
}


// ==================================================================
// misc
// ==================================================================

inline
float lerp(float a, float b, float x) {
  return x*a + (1.f-x)*b;
}

inline
vec3f lerp(vec3f a, vec3f b, float x) {
  return x*a + (1.f-x)*b;
}

inline
float clamp(float x, float a, float b) {
  return fmaxf(a,fminf(x,b));
}

inline
vec3f clamp(vec3f x, vec3f a, vec3f b) {
  return max(a,min(x,b));
}

inline
size_t linearIndex(int x, int y, int z, int dims[3]) {
  return z*dims[0]*dims[1] + y*size_t(dims[0]) + x;
}


// ==================================================================
// sliceT
// ==================================================================

template <typename VecT>
struct sliceT
{
  typedef typename VecT::value_type value_type;

  size_t lower, upper;
  VecT &vec;

  inline
  size_t size() const;

  template <typename VecT2>
  inline
  sliceT &operator=(const sliceT<VecT2> &other);

  inline
  value_type &operator[](size_t i);

  inline
  const value_type &operator[](size_t i) const;
};

template <typename VecT>
size_t sliceT<VecT>::size() const {
  return upper-lower;
}

template <typename VecT>
template <typename VecT2>
sliceT<VecT> &sliceT<VecT>::operator=(const sliceT<VecT2> &other) {
  assert(size()==other.size());

  if (&other != (const sliceT<VecT2> *)this) {
    for (size_t i=0; i<size(); ++i) {
      (*this)[i] = other[i];
    }
  }

  return *this;
}

template <typename VecT>
typename sliceT<VecT>::value_type &sliceT<VecT>::operator[](size_t i) {
  return vec[lower+i];
}

template <typename VecT>
const typename sliceT<VecT>::value_type &sliceT<VecT>::operator[](size_t i) const {
  return vec[lower+i];
}


// ==================================================================
// blockT
// ==================================================================

template <typename MatT>
struct blockT
{
  typedef typename MatT::value_type value_type;

  vec2ui lower, upper;
  MatT &mat;

  inline
  vec2ui size() const;

  template <typename MatT2>
  inline
  blockT &operator=(const blockT<MatT2> &other);

  inline
  value_type &operator()(unsigned x, unsigned y);

  inline
  const value_type &operator()(unsigned x, unsigned y) const;
};

template <typename MatT>
vec2ui blockT<MatT>::size() const {
  return upper-lower;
}

template <typename MatT>
template <typename MatT2>
blockT<MatT> &blockT<MatT>::operator=(const blockT<MatT2> &other) {
  assert(size()==other.size());

  if (&other != (const blockT<MatT2> *)this) {
    for (unsigned y=0; y<size().y; ++y) {
      for (unsigned x=0; x<size().x; ++x) {
        (*this)(x,y) = other(x,y);
      }
    }
  }

  return *this;
}

template <typename MatT>
typename blockT<MatT>::value_type &blockT<MatT>::operator()(unsigned x, unsigned y) {
  return mat(lower.x+x,lower.y+y);
}

template <typename MatT>
const typename blockT<MatT>::value_type &blockT<MatT>::operator()(unsigned x, unsigned y) const {
  return mat(lower.x+x,lower.y+y);
}


// ==================================================================
// variable-size vector type
// ==================================================================

template <typename T, typename Allocator>
struct vectorN
{
  typedef T value_type;

  vectorN() = default;
  vectorN(size_t N);
  vectorN(const vectorN &other);
  vectorN &operator=(const vectorN &other);
 ~vectorN();
  size_t N;
  T *data=nullptr;
  Allocator alloc;
  static_assert(std::is_same<T,typename Allocator::value_type>::value,"Type mismatch");

  inline
  size_t size() const;

  inline
  value_type &operator[](size_t i);

  inline
  const value_type &operator[](size_t i) const;

  inline
  sliceT<vectorN> slice(size_t lower, size_t upper);

  inline
  sliceT<const vectorN> slice(size_t lower, size_t upper) const;
};

template <typename T, typename Allocator>
vectorN<T,Allocator>::vectorN(size_t N)
  : N(N)
{
  data = alloc.allocate(N);
}

template <typename T, typename Allocator>
vectorN<T,Allocator>::vectorN(const vectorN &other)
  : N(other.N)
{
  data = alloc.allocate(N);
  for (size_t i=0; i<N; ++i) {
    data[i] = other.data[i];
  }
}

template <typename T, typename Allocator>
vectorN<T,Allocator> &vectorN<T,Allocator>::operator=(const vectorN &other) {
  if (&other != this) {
    alloc.deallocate(data,N);
    N = other.N;
    data = alloc.allocate(N);
    for (size_t i=0; i<N; ++i) {
      data[i] = other.data[i];
    }
  }
  return *this;
}

template <typename T, typename Allocator>
vectorN<T,Allocator>::~vectorN() {
  alloc.deallocate(data,N);
}

template <typename T, typename Allocator>
size_t vectorN<T,Allocator>::size() const {
  return N;
}

template <typename T, typename Allocator>
typename vectorN<T,Allocator>::value_type &vectorN<T,Allocator>::operator[](size_t i) {
  assert(i<N);
  return data[i];
}

template <typename T, typename Allocator>
const typename vectorN<T,Allocator>::value_type &vectorN<T,Allocator>::operator[](size_t i) const {
  assert(i<N);
  return data[i];
}

template <typename T, typename Allocator>
sliceT<vectorN<T,Allocator>> vectorN<T,Allocator>::slice(size_t lower, size_t upper) {
  return {lower,upper,*this};
}

template <typename T, typename Allocator>
sliceT<const vectorN<T,Allocator>> vectorN<T,Allocator>::slice(size_t lower, size_t upper) const {
  return {lower,upper,*this};
}

template <typename T, typename Allocator>
vectorN<T,Allocator> operator-(const vectorN<T,Allocator> &u) {
  vectorN<T,Allocator> result(u.size());
  for (size_t i=0; i<result.size(); ++i) {
    result[i] = -u[i];
  }
  return result;
}

template <typename T, typename Allocator>
vectorN<T,Allocator> operator+(const vectorN<T,Allocator> &u, const vectorN<T,Allocator> &v) {
  assert(u.size()==v.size());
  vectorN<T,Allocator> result(u.size());
  for (size_t i=0; i<result.size(); ++i) {
    result[i] = u[i]+v[i];
  }
  return result;
}

template <typename T, typename Allocator>
vectorN<T,Allocator> operator-(const vectorN<T,Allocator> &u, const vectorN<T,Allocator> &v) {
  assert(u.size()==v.size());
  vectorN<T,Allocator> result(u.size());
  for (size_t i=0; i<result.size(); ++i) {
    result[i] = u[i]-v[i];
  }
  return result;
}

template <typename T, typename Allocator>
vectorN<T,Allocator> operator*(const vectorN<T,Allocator> &u, const vectorN<T,Allocator> &v) {
  assert(u.size()==v.size());
  vectorN<T,Allocator> result(u.size());
  for (size_t i=0; i<result.size(); ++i) {
    result[i] = u[i]*v[i];
  }
  return result;
}

template <typename T, typename Allocator>
vectorN<T,Allocator> operator/(const vectorN<T,Allocator> &u, const vectorN<T,Allocator> &v) {
  assert(u.size()==v.size());
  vectorN<T,Allocator> result(u.size());
  for (size_t i=0; i<result.size(); ++i) {
    result[i] = u[i]/v[i];
  }
  return result;
}

template <typename T, typename Allocator>
vectorN<T,Allocator> operator*(const vectorN<T,Allocator> &u, const T &a) {
  vectorN<T,Allocator> result(u.size());
  for (size_t i=0; i<result.size(); ++i) {
    result[i] = u[i]*a;
  }
  return result;
}

template <typename T, typename Allocator>
vectorN<T,Allocator> operator/(const vectorN<T,Allocator> &u, const T &a) {
  vectorN<T,Allocator> result(u.size());
  for (size_t i=0; i<result.size(); ++i) {
    result[i] = u[i]/a;
  }
  return result;
}

template <typename T, typename Allocator>
vectorN<T,Allocator> operator+=(vectorN<T,Allocator> &u, const vectorN<T,Allocator> &v) {
  u = u + v;
  return u;
}

template <typename T, typename Allocator>
T dot(const vectorN<T,Allocator> &u, const vectorN<T,Allocator> &v) {
  assert(u.size()==v.size());
  T result(0.0);
  for (size_t i=0; i<u.size(); ++i) {
    result += u[i]*v[i];
  }
  return result;
}

template <typename T, typename Allocator>
T length(const vectorN<T,Allocator> &u) {
  return sqrtf(dot(u,u));
}

template <typename T, typename Allocator>
vectorN<T,Allocator> normalize(const vectorN<T,Allocator> &u) {
  return u / sqrtf(dot(u,u));
}

template <typename T, typename Allocator>
size_t arg_min(const vectorN<T,Allocator> &u) {
  size_t biggestDim = 0;
  for (size_t i=1; i<u.size(); ++i)
    if (u[i] < u[biggestDim]) biggestDim = i;
  return biggestDim;
}

template <typename T, typename Allocator>
size_t arg_max(const vectorN<T,Allocator> &u) {
  size_t biggestDim = 0;
  for (size_t i=1; i<u.size(); ++i)
    if (u[i] > u[biggestDim]) biggestDim = i;
  return biggestDim;
}

template <typename T, typename Allocator>
vectorN<T,Allocator> clamp(const vectorN<T,Allocator> &u, const T &a, const T &b) {
  vectorN<T,Allocator> result(u.size());
  for (size_t i=0; i<u.size(); ++i) {
    result[i] = u[i];
    result[i] = result[i] < a ? a : result[i];
    result[i] = result[i] > b ? b : result[i];
  }
  return result;
}

template <typename T, typename Allocator>
std::ostream &operator<<(std::ostream &out, const vectorN<T,Allocator> &v) {
  out << '(';
  for (size_t i=0; i<v.N; ++i) {
    out << v[i];
    if (i < v.N-1) out << ',';
  }
  out << ')';
  return out;
}


// ==================================================================
// variable-size matrix type
// ==================================================================

template <typename T, typename Allocator>
struct matrixN
{
  typedef T value_type;

  matrixN() = default;
  matrixN(unsigned numRows, unsigned numCols);
  matrixN(const matrixN &other);
  matrixN &operator=(const matrixN &other);
 ~matrixN();
  unsigned numRows=0, numCols=0;
  T *data=nullptr;
  Allocator alloc;
  static_assert(std::is_same<T,typename Allocator::value_type>::value,"Type mismatch");

  inline
  vec2ui size() const;

  inline
  value_type &operator[](const vec2ui &i);

  inline
  const value_type &operator[](const vec2ui &i) const;

  inline
  value_type &operator()(unsigned x, unsigned y);

  inline
  const value_type &operator()(unsigned x, unsigned y) const;

  inline
  blockT<matrixN> block(vec2ui lower, vec2ui upper);

  inline
  blockT<const matrixN> block(vec2ui lower, vec2ui upper) const;
};

template <typename T, typename Allocator>
matrixN<T,Allocator>::matrixN(unsigned numRows, unsigned numCols)
  : numRows(numRows), numCols(numCols)
{
  data = alloc.allocate(numRows*size_t(numCols));
}

template <typename T, typename Allocator>
matrixN<T,Allocator>::matrixN(const matrixN &other)
  : numRows(other.numRows), numCols(other.numCols)
{
  data = alloc.allocate(numRows*size_t(numCols));
  for (unsigned y=0; y<numRows; ++y) {
    for (unsigned x=0; x<numCols; ++x) {
      (*this)(x,y) = other(x,y);
    }
  }
}

template <typename T, typename Allocator>
matrixN<T,Allocator> &matrixN<T,Allocator>::operator=(const matrixN &other)
{
  if (&other != this) {
    alloc.deallocate(data,numRows*size_t(numCols));
    numRows = other.numRows;
    numCols = other.numCols;
    data = alloc.allocate(numRows*size_t(numCols));
    for (unsigned y=0; y<numRows; ++y) {
      for (unsigned x=0; x<numCols; ++x) {
        (*this)(x,y) = other(x,y);
      }
    }
  }
  return *this;
}

template <typename T, typename Allocator>
matrixN<T,Allocator>::~matrixN() {
  static_assert(std::is_same<T,typename Allocator::value_type>::value,"Type mismatch");
  alloc.deallocate(data,numRows*size_t(numCols));
}

template <typename T, typename Allocator>
vec2ui matrixN<T,Allocator>::size() const {
  return {numRows,numCols};
}

template <typename T, typename Allocator>
typename matrixN<T,Allocator>::value_type &matrixN<T,Allocator>::operator[](const vec2ui &i) {
  assert(i.x<numCols && i.y<numRows);
  return data[i.y*size_t(numCols)+i.x];
}

template <typename T, typename Allocator>
const typename matrixN<T,Allocator>::value_type &matrixN<T,Allocator>::operator[](const vec2ui &i) const {
  assert(i.x<numCols && i.y<numRows);
  return data[i.y*size_t(numCols)+i.x];
}

template <typename T, typename Allocator>
typename matrixN<T,Allocator>::value_type &matrixN<T,Allocator>::operator()(unsigned x, unsigned y) {
  assert(x<numCols && y<numRows);
  return data[y*size_t(numCols)+x];
}

template <typename T, typename Allocator>
const typename matrixN<T,Allocator>::value_type &matrixN<T,Allocator>::operator()(unsigned x, unsigned y) const {
  assert(x<numCols && y<numRows);
  return data[y*size_t(numCols)+x];
}

template <typename T, typename Allocator>
blockT<matrixN<T,Allocator>> matrixN<T,Allocator>::block(vec2ui lower, vec2ui upper) {
  return {lower,upper,*this};
}

template <typename T, typename Allocator>
blockT<const matrixN<T,Allocator>> matrixN<T,Allocator>::block(vec2ui lower, vec2ui upper) const {
  return {lower,upper,*this};
}

template <typename T, typename Allocator>
matrixN<T,Allocator> transpose(const matrixN<T,Allocator> &m) {
  matrixN<T,Allocator> result(m.numCols,m.numRows);
  for (unsigned y=0; y<m.numRows; ++y) {
    for (unsigned x=0; x<m.numCols; ++x) {
      result(y,x) = m(x,y);
    }
  }
  return result;
}

template <typename T, typename Allocator>
vectorN<T,Allocator> operator*(const vectorN<T,Allocator> &v, const matrixN<T,Allocator> &m) {
  assert(v.N==m.numRows);
  vectorN<T,Allocator> result(m.numCols);
  for (size_t i=0; i<result.size(); ++i) {
    result[i] = T(0.0);
    for (size_t j=0; j<v.size(); ++j) {
      result[i] += v[j]*m(i,j);
    }
  }
  return result;
}

template <typename T, typename Allocator>
vec2ui arg_min(const matrixN<T,Allocator> &m) {
  vec2ui biggestDim(0,0);
  for (unsigned y=0; y<m.numRows; ++y) {
    for (unsigned x=0; x<m.numCols; ++x) {
      if (m(x,y) < m[biggestDim])
        biggestDim = vec2ui(x,y);
    }
  }
  return biggestDim;
}

template <typename T, typename Allocator>
vec2ui arg_max(const matrixN<T,Allocator> &m) {
  vec2ui biggestDim(0,0);
  for (unsigned y=0; y<m.numRows; ++y) {
    for (unsigned x=0; x<m.numCols; ++x) {
      if (m(x,y) > m[biggestDim])
        biggestDim = vec2ui(x,y);
    }
  }
  return biggestDim;
}

template <typename T, typename Allocator>
std::ostream &operator<<(std::ostream &out, const matrixN<T,Allocator> &m) {
  out << '(';
  for (unsigned y=0; y<m.numRows; ++y) {
    out << "col[" << y << "]:(";
    for (unsigned x=0; x<m.numCols; ++x) {
      out << m(x,y);
      if (x < m.numCols-1) out << ',';
    }
    out << ')';
    if (y < m.numRows-1) out << ',';
  }
  out << ')';
  return out;
}


// ==================================================================
// ray tracing
// ==================================================================

struct Ray
{
  vec3f org;
  float tmin;
  vec3f dir;
  float tmax;
};

inline
bool boxTest(const Ray &ray, const box3f &box, float &t0, float &t1) {
  const vec3f t_lo = (box.lower - ray.org) / ray.dir;
  const vec3f t_hi = (box.upper - ray.org) / ray.dir;

  const vec3f t_nr = min(t_lo,t_hi);
  const vec3f t_fr = max(t_lo,t_hi);

  t0 = fmaxf(ray.tmin,reduce_max(t_nr));
  t1 = fminf(ray.tmax,reduce_min(t_fr));
  return t0 < t1;
}

} // namespace math


