// =========================================================================
// File:       tledStaticBVHCreator.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2015
// 
// Copyright (c) 2015, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledStaticBVHCreator_H
#define tledStaticBVHCreator_H

#include "tledBVHCreator.h"

/**
 * \brief BVH generator for rigid geometry
 * \ingroup contact
 *
 * Primarily intended for use with AABB BV geometry.
 */
template <class TBVH>
class tledStaticBVHCreator : public tledBVHCreatorImpl<TBVH, typename TBVH::ContactMesh> {
  /**
   * \name BVH Generation
   * @{
   */
private:
  void _SplitBinary(std::vector<int>::iterator ppi_childIndexBounds[][2], const std::vector<int>::iterator &pIndicesBegin, const std::vector<int>::iterator &pIndicesEnd);

protected:
  float FindCentroidAxisAvg(const int cInd, const std::vector<int>::const_iterator ic_pIndsBegin, const std::vector<int>::const_iterator ic_pIndsEnd) const;
  int FindCentroidMaxAxis(const std::vector<int>::const_iterator ic_pIndsBegin, const std::vector<int>::const_iterator ic_pIndsEnd) const;

  virtual void GenerateMain(void);
  virtual void SplitBV(std::vector<int>::iterator p_childIndexBounds[][2], const int BVIndex, const std::vector<int>::iterator &pIndicesBegin, const std::vector<int>::iterator &pIndicesEnd);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledStaticBVHCreator(void) {}
  /** @} */
}; 

#include "tledStaticBVHCreator.tpp"
#endif
