// =========================================================================
// File:       tledBVHTraverser.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVHTraverser_H
#define tledBVHTraverser_H

class tledUnstructuredContactManager;

/**
 * \brief Base API for BVH traversal algorithms.
 * \ingroup contact
 */
class tledBVHTraverser {
  /**
   * \name Collision Detection
   * @{
   */
protected:
  virtual void RunNarrowPhase(void) = 0;
  virtual void RunBroadPhase(void) = 0;

public:
  virtual void FindCollisions(void) = 0;
  /** @} */

  /**
   * \name Master/Slave Control
   * @{
   */
protected:
  bool m_DoMaster;

public:
  /**
   * \brief Master/slave switch, if set to true the surface bounded by the BVH will be the master in the next collision detection pass. 
   */
  void SetDoMaster(const bool doMaster) { m_DoMaster = doMaster; }

  /**
   * \return Master/slave flag, if true master.
   */
  bool DoMaster(void) const { return m_DoMaster; }
  /** @} */

  /**
   * \name Safety Margins and Related Constants
   * @{
   */
private:
  float m_ContactMaxDist;
  
public:
  /**
   * \brief Maximum geometry distance at which an event gets passed on to the contact solver (important for rate constraints).
   *
   * Default: 1/2*real BV margin
   */
  void SetNarrowPhaseMaxDistance(const float maxDist) { m_ContactMaxDist = maxDist; }

  /**
   * \brief Maximum geometry distance at which an event gets passed on to the contact solver (important for rate constraints).
   */
  float GetNarrowPhaseMaxDistance(void) const { return m_ContactMaxDist; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  /** Final initialisation after BVH has been constructed etc. */
  virtual void Init(tledUnstructuredContactManager &r_manager) = 0;

  tledBVHTraverser(void) : m_DoMaster(true) {}
  virtual ~tledBVHTraverser(void) {}
  /** @} */
};

/**
 * \brief Base class for BVH traversal algorithm implementations
 * \ingroup contact
 */
template <class TMasterBVH, class TSlaveBVH, class TAPI>
class tledBVHTraverserImpl : public TAPI {
  /**
   * \name Types
   * @{
   */
public:
  typedef TAPI Superclass;
  typedef TMasterBVH MasterBVH;
  typedef TSlaveBVH SlaveBVH;
  typedef typename TMasterBVH::ContactMesh MasterMesh;
  typedef typename TSlaveBVH::ContactMesh SlaveMesh;
  /** @} */

  /**
   * \name BVH Access
   * @{
   */
private:
  const MasterBVH *mpc_MasterBVH;
  SlaveBVH *mp_SlaveBVH;

public:
  const MasterBVH& GetMasterBVH(void) const { return *mpc_MasterBVH; }

  const SlaveBVH& GetSlaveBVH(void) const { return *mp_SlaveBVH; }
  SlaveBVH& GetSlaveBVH(void) { return *mp_SlaveBVH; }
  /** @} */

  /**
   * \name Detection
   * @{
   */
public:
  virtual void FindCollisions(void);
  /** @} */

  /**
   * \name Mesh Access
   * @{
   */
public:
  const MasterMesh& GetMasterMesh(void) const { return GetMasterBVH().GetMesh(); }
  
  const SlaveMesh& GetSlaveMesh(void) const { return GetSlaveBVH().GetMesh(); }
  SlaveMesh& GetSlaveMesh(void) { return GetSlaveBVH().GetMesh(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledBVHTraverserImpl(SlaveBVH &r_slaveBVH, const MasterBVH &masterBVH) : mpc_MasterBVH(&masterBVH), mp_SlaveBVH(&r_slaveBVH) {}
  virtual ~tledBVHTraverserImpl(void) {}
  /** @} */  
};

#include "tledBVHTraverser.tpp"
#endif
