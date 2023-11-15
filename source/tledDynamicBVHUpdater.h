// =========================================================================
// File:       tledDynamicBVHUpdater.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDynamicBVHUpdater_H
#define tledDynamicBVHUpdater_H

class tledDynamicBVH;

/**
 * \brief Dynamic-BVH Updater API
 * \ingroup contact
 */
class tledDynamicBVHUpdater {
  /**
   * \name BVH
   * @{
   */
public:
  /** Access to BVH through reference of most general type */
  virtual tledDynamicBVH& GetUnspecifiedBVH(void) = 0;
  virtual const tledDynamicBVH& GetUnspecifiedBVH(void) const = 0;

  virtual void SetBVH(tledDynamicBVH &r_bvh) = 0;
  /** @} */

  /**
   * \name Updating
   * @{
   */
public:
  virtual void UpdateBVH(void) = 0;
  /** @} */

  /**
   * \name Initialisation
   * @{
   */
public:
  virtual void Init(void) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledDynamicBVHUpdater(void) {}

public:
  virtual ~tledDynamicBVHUpdater(void) {}
  /** @} */
};

/**
 * \brief Base class for dynamic-BVH updaters
 * \ingroup contact
 */
template <class TBVH>
class tledDynamicBVHUpdaterImpl : public tledDynamicBVHUpdater {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBVH BVH;
  typedef typename TBVH::BoundingVolume BoundingVolume;
  typedef typename TBVH::ContactMesh ContactMesh;
  typedef typename ContactMesh::Facet Facet;
  /** @} */

  /**
   * \name BVH
   * @{
   */
private:
  BVH *mp_BVH;
  
public:
  const BVH& GetBVH(void) const { return *mp_BVH; }
  BVH& GetBVH(void) { return *mp_BVH; }
  virtual void SetBVH(tledDynamicBVH &r_bvh) { mp_BVH = static_cast<BVH*>(&r_bvh); }

  virtual tledDynamicBVH& GetUnspecifiedBVH(void) { return GetBVH(); }
  virtual const tledDynamicBVH& GetUnspecifiedBVH(void) const { return GetBVH(); }
  /** @} */

  /**
   * \name Updating
   * @{
   */
protected:
#ifndef NDEBUG
  /** Test routine (only availble in debug builds) that tests that the nodes of a given subtree are properly contained by the BVH when translated by trans */
  void CheckContainmentRecursive(std::vector<int> &r_nodeInds, const int subtreeRoot, const float trans[]) const;  

  /** Test routine (only availble in debug builds) that tests that the nodes of a given subtree are properly contained by the BVH when translated by trans */
  void CheckContainmentRecursive(const int subtreeRoot, const float trans[]) const;

  /** Test routine (only availble in debug builds) that tests that the nodes of a given subtree are properly contained by the BVH */
  void CheckContainmentRecursive(const int subtreeRoot) const;
#endif
  /** @} */

  /**
   * \name Mesh
   * @{
   */
public:
  const ContactMesh& GetMesh(void) const { return this->GetBVH().GetMesh(); }
  ContactMesh& GetMesh(void) { return this->GetBVH().GetMesh(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledDynamicBVHUpdaterImpl(BVH &r_bvh) : mp_BVH(&r_bvh) {}
  tledDynamicBVHUpdaterImpl(void) : mp_BVH(NULL) {}
  virtual ~tledDynamicBVHUpdaterImpl(void) {}
  /** @} */
};

#include "tledDynamicBVHUpdater.tpp"
#endif
