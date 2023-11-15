// =========================================================================
// File:       tledHelper.h
// Purpose:    Class providing general purpose routines
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledHelper_H
#define tledHelper_H

#ifdef _Visualisation_
#include <vtkVersion.h>
#endif

#include <cmath>
#include <map>
#include <functional>
#include <algorithm>
#include <vector>
#include <limits>
#include <sstream>
#include <iostream>
#ifndef __CUDACC__
#if defined __GNUC__ && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 2))
#include <ext/hash_map>
#elif defined  __GNUC__ && !defined __llvm__ && !defined __GXX_EXPERIMENTAL_CXX0X__
#include <tr1/unordered_map>
namespace tledCXX0XNamespace = std::tr1;
#else
/* If on VC++ 2008, download "feature pack" from http://www.microsoft.com/download/en/details.aspx?displaylang=en&id=6922 */
#include <unordered_map>
namespace tledCXX0XNamespace = std;
#endif
#endif /* #ifndef __CUDACC__ */

#ifdef __GNUG__
#define __tled_inline extern inline
#else
#define __tled_inline inline
#endif

/**
 * \defgroup helper Helper Routines and Data Structures
 *
 * General purpose routines.
 */

#define tledPi 3.141592653589793238462643383279f

#ifdef _Visualisation_
#if VTK_MAJOR_VERSION >= 6
#define tledVTK6CompatSetInput(P_VTKOBJ, P_INPUT) P_VTKOBJ->SetInputData(P_INPUT);
#else
#define tledVTK6CompatSetInput(P_VTKOBJ, P_INPUT) P_VTKOBJ->SetInput(P_INPUT);			    
#endif 
#endif

/**
 * \brief Misc. helper routines
 * \ingroup helper 
 */
namespace tledHelper {
  /**
   * \brief Returns a vector with sorted, unique entries for any input vector with user-defined equality and ordering operators
   */
  template <typename TValue, class TEqualPredicate, class TOrdering>
  std::vector<TValue> MakeSortedUnique(const std::vector<TValue> &vec, const TEqualPredicate &equal, const TOrdering &order) {
    std::vector<TValue> uniqueVec(vec);
    typename std::vector<TValue>::iterator i_end;

    std::sort(uniqueVec.begin(), uniqueVec.end(), order);
    i_end = std::unique(uniqueVec.begin(), uniqueVec.end(), equal);    
    uniqueVec.erase(i_end, uniqueVec.end());

    return uniqueVec;
  } /* MakeSortedUnique */

  /**
   * \brief Returns a vector with sorted, unique entries for any input vector (req.: binary ==, < operators for element type)
   */
  template <typename TValue>
  std::vector<TValue> MakeSortedUnique(const std::vector<TValue> &vec) {
    return MakeSortedUnique(vec, std::equal_to<TValue>(), std::less<TValue>());
  } /* MakeSortedUnique */

  /** \brief Returns true if v is numerically insignificant compared to ref. */
  template <typename TFloat>
  bool IsNumericallyZero(const TFloat &v, const TFloat &ref) {
#ifndef __USE_FAST_ZERO_TEST
    return v == 0 || (ref != 0 && std::fabs(v/ref) < std::numeric_limits<TFloat>::epsilon());
#else
    /* seems to lead to problems with a number of compilers (always false?) */
    return ref + v == ref;
#endif
  }

  class Error {
  private:
    std::ostringstream m_MsgStream;

  protected:
    virtual bool IsFatal(void) const = 0;
    virtual std::ostringstream& GetMsgStream(void) { return m_MsgStream; }
    virtual std::ostream& GetOutStream(void) { return std::cerr; }
    virtual const char* GetPrefix(void) const = 0;

    /**
     * \name std::ostream-like API
     * @{
     */
  public:
    Error& operator<<(const char c);
    Error& operator<<(const std::string &str);
    Error& operator<<(const char *str) { return this->operator<<(std::string(str)); }
    Error& operator<<(const int num);
    Error& operator<<(const double num);
    Error& operator<<(const void *ptr);    
    /** @} */

  public:
    static void Log(Error &r_error, const char *whereFile, const int whereLine, const char *whereFnct);

  public:
    Error(void) {}
    virtual ~Error(void) {}
  };

  class FatalError : public Error {
  public:
    virtual bool IsFatal(void) const { return true; }
    virtual const char* GetPrefix(void) const { return "Fatal error"; }
  };

  class NonFatalError : public Error {
  public:
    virtual bool IsFatal(void) const { return false; }
    virtual const char* GetPrefix(void) const { return "Error"; }
  };

  class Info : public Error {
  protected:
    virtual bool IsFatal(void) const { return false; }
    virtual std::ostream& GetOutStream(void) { return std::cout; }
    virtual const char* GetPrefix(void) const { return "Info"; }
  };

  class Warning : public NonFatalError {
  protected:
    virtual const char* GetPrefix(void) const { return "Warning"; }
  };
}

#define tledLogErrorStream(STREAM) tledHelper::Error::Log((STREAM), __FILE__, __LINE__, __FUNCTION__)
#define tledFatalError(MSG) tledLogErrorStream(tledHelper::FatalError() << MSG)
#define tledNonFatalError(MSG) tledLogErrorStream(tledHelper::NonFatalError() << MSG)
#define tledWarning(MSG) tledLogErrorStream(tledHelper::Warning() << MSG)
#define tledFatalNotYetImplementedError tledFatalError("Feature not yet implemented")
#define tledFatalFeatureNotEnabledError tledFatalError("Feature not enabled at compile time.")

#ifndef NDEBUG
/** Only works in debug builds */
#define tledLogDebugStream(STREAM) tledLogErrorStream(STREAM)
#else
/** Only works in debug builds */
#define tledLogDebugStream(STREAM)
#endif

/** 
 * \brief Generates a sequence of integers
 * \ingroup helper
 */
class tledSequenceGenerator {
private:
  int m_Counter;

  /**
   * \name Generator API
   * @{
   */
public:
  int operator()(void) { return (m_Counter += 1) - 1; }
  /** @} */

  /**
   * \name Sequence Creation
   * @{
   */
public:
  static std::vector<int> MakeSequence(const int start, const int end) {
    std::vector<int> seq(end - start);

    tledSequenceGenerator(seq.begin(), seq.end(), start);

    return seq;
  }

public:
  /** Generates an integer sequence {startValue, .., (rangeEnd - 1 - rangeBegin) + startValue} */
  template <typename TIterator>
  tledSequenceGenerator(TIterator rangeBegin, TIterator rangeEnd, const int startValue) : m_Counter(startValue) { std::generate(rangeBegin, rangeEnd, *this); }
  /** @} */
};

/**
 * \brief Simple fixed length array
 * \ingroup helper
 *
 *
 * Emulates C++0x std::array
 */
template <typename TValue, const int t_numElements>
class tledArray {
  /**
   * \name Element Access 
   * @{
   */
public:
  typedef const TValue* const_iterator;
  typedef TValue* iterator;

private:
  TValue m_Array[t_numElements];

public:
  const TValue& operator[](const int eInd) const { return m_Array[eInd]; }
  TValue& operator[](const int eInd) { return m_Array[eInd]; } 

  bool operator==(const tledArray &x) const { return std::equal(this->begin(), this->end(), x.begin()); }

  iterator begin(void) { return m_Array; }
  const_iterator begin(void) const { return m_Array; }

  iterator end(void) { return m_Array + t_numElements; }
  const_iterator end(void) const { return m_Array + t_numElements; }  
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledArray(void) {}
  tledArray(const TValue initVals[]) { std::copy(initVals, initVals + t_numElements, m_Array); }
  /** @} */
};

#ifndef __CUDACC__
/**
 * \brief General purpose hash maps
 * \ingroup helper
 *
 * Not available in CUDA modules. 
 */
template <typename TKeyType, typename TPayloadType, class THasherType>
#if defined __GNUC__ && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 2))
class tledMap : public __gnu_cxx::hash_map<TKeyType, TPayloadType, THasherType>
#else
class tledMap : public tledCXX0XNamespace::unordered_map<TKeyType, TPayloadType, THasherType, std::equal_to<TKeyType> >
#endif
{};
#endif /* #ifndef __CUDACC__ */

#ifdef _MSC_VER
namespace std {
  template <typename T>
  inline bool isnan(T x) {
    return _isnan(x) != 0;
  }

  /**
   * \deprecated use algorithm's max
   */
  template <typename T>
  inline T fmax(T x, T y) {
    return max(x, y);
  }

  /**
   * \deprecated use algorithm's min
   */
  template <typename T>
  inline T fmin(T x, T y) {
    return min(x, y);
  }
}
#endif

#endif // tledHelper_H
