// =========================================================================
// File:       tledMembraneSelfCollisionBVHCreator.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledMembraneSelfCollisionBVHCreator_H
#define tledMembraneSelfCollisionBVHCreator_H

#include "tledHelper.h"
#include "tledSelfCollisionBVHCreator.h"

/**
 * \brief BVH creator adapted for use when membranes are involved in contact problem.
 */
template <class TBVH>
class tledMembraneSelfCollisionBVHCreator : public tledSelfCollisionBVHCreator<TBVH> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledSelfCollisionBVHCreator<TBVH> Superclass;
  typedef typename TBVH::ContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  typedef typename TBVH::BoundingVolume BoundingVolume;
  /** @} */

  /**
   * \name BVH Generation
   * @{
   */
private:
  std::vector<int> m_MembraneSeedFacetInds;

private:
  int _FindMembraneClusterRoot(const int primitiveInd) const;

protected:
  virtual void PreBottomUpHook(std::vector<int> &r_activeBVs);
  virtual void GenerateMain(void);
  /** @} */

public:
  virtual ~tledMembraneSelfCollisionBVHCreator(void) {}
};

#include "tledMembraneSelfCollisionBVHCreator.tpp"
#endif
