// =========================================================================
// File:       tledDynamicBVH.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    September 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDynamicBVH_H
#define tledDynamicBVH_H

#include "tledBVH.h"
#include "tledMovingRigidContactSurface.h"

#include <vector>

class tledDynamicBVHUpdater;

/**
 * \brief Basic API for dynamic BVHs, i.e. such bounding moving geometry
 * \ingroup contact
 */
class tledDynamicBVH : public tledBVH {
  /**
   * \name Refitting
   * @{
   */
private:
  tledDynamicBVHUpdater *mp_Updater;

public:
  /**
   * \brief Adapts the bounds of the BVs to accomodate geometry changes.
   */
  virtual void Update(void) = 0;

  tledDynamicBVHUpdater& GetUpdater() { return *mp_Updater; }
  const tledDynamicBVHUpdater& GetUpdater() const { return *mp_Updater; }

  /** The updater is assumed to be an object residing in dynamically allocated memory, is destroyed with delete together with the BVH */
  void SetUpdater(tledDynamicBVHUpdater &r_updater);

  /** Standard top-down refit starting at BVH node bvInd. */
  virtual void UpdateTopDownRecursive(const int bvInd) = 0;
  /** @} */

  /**
   * \name Safety Margins and Related Constants
   * @{
   */
private:
  float m_BVMaxDisplacement;

public:
  virtual void SetMargin(const float bvMargin);

  /**
   * \brief Getter for maximum displacement allowed without some kind of updating of BVs/BVH subtree.
   */
  float GetBVMaxDisplacement(void) const { return m_BVMaxDisplacement; }

  virtual void SetBVMaxDisplacement(const float maxDisp);

  /** \brief BV safety margin that is maintained at all times. */
  float GetRealMargin(void) const { return GetMargin() - GetBVMaxDisplacement(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  static tledDynamicBVH* CreateBVH(const tledRigidContactSurface &mesh, const std::string &bvType, const float margin, const bool useGPU);

  /** XML import */
  static tledDynamicBVH* CreateBVH(const tledRigidContactSurface &mesh, const XMLNode root, const bool useGPU);

  tledDynamicBVH(void) : mp_Updater(NULL), m_BVMaxDisplacement(std::numeric_limits<float>::quiet_NaN()) {}
  virtual ~tledDynamicBVH(void);
  /** @} */
};

/**
 * \brief Base class for BVH implementations for dynamic geometry
 * \ingroup contact
 */
template <class TSurface, class TBV, class TAPI = tledDynamicBVH>
class tledDynamicBVHImpl : public tledBVHImpl<TSurface, TBV, TAPI> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBVHImpl<TSurface, TBV, TAPI> Superclass;
  typedef TBV BoundingVolume;
  typedef TSurface ContactMesh;
  /** @} */

  /**
   * \name Mesh
   * @{
   */
public:
  const ContactMesh& GetMesh(void) const { return Superclass::GetMesh(); }
  ContactMesh& GetMesh(void) { return const_cast<ContactMesh&>(Superclass::GetMesh()); }
  /** @} */

  /**
   * \name Update
   * @{
   */
private:
  int m_UpdateCounter;

public:
  virtual void Update(void);  

  /** Applies a translation to a subtree */
  virtual void TranslateSubTree(const int rootBVInd, const float t[]);

  /** Applies a rigid transform to a subtree */
  virtual void TransformSubTree(const int rootBVInd, const float m[][3], const float cor[], const float t[]);

  /** Standard top-down refit starting at BVH node bvInd. */
  virtual void UpdateTopDownRecursive(const int bvInd);

  virtual void ComputeBoundsFromNodes(BoundingVolume &r_bv, const int *nodeListStart, const int *nodeListEnd) const;

  int GetUpdateCounter(void) const { return m_UpdateCounter; }
  void ResetUpdateCounter(void) { m_UpdateCounter = 0; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void Init(tledBVHCreator &r_builder);

  tledDynamicBVHImpl(ContactMesh &r_mesh) : Superclass(r_mesh), m_UpdateCounter(0) {}
  virtual ~tledDynamicBVHImpl(void) {}
  /** @} */
};

#include "tledDynamicBVH.tpp"
#endif
