// =========================================================================
// File:       tledDeformableMembraneContactSurfaceCPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDeformableMembraneContactSurfaceCPU_H
#define tledDeformableMembraneContactSurfaceCPU_H

#include "tledDeformableMembraneContactSurface.h"
#include "tledDeformableContactSurfaceCPU.h"
#include "tledVectorArithmetic.h"

/**
 * \brief Untemplated API for contact surfaces comprising a membrane component
 * \ingroup contact
 */
class tledDeformableMembraneContactSurfaceCPU : public tledDeformableContactSurfaceCPU, public tledDeformableMembraneContactSurface {
  /**
   * \name Membrane Specific API
   * @{
   */
private:
  const float *mpc_MembraneThicknesses;
  
protected:
  float GetMembraneElementThickness(const int meInd) const { return mpc_MembraneThicknesses[meInd]; }

public:
  virtual void SetMembraneFacetThickenesses(const float ts[]) { mpc_MembraneThicknesses = ts; }
  virtual void UpdateMembraneFacetNormal(const int membraneFacetIndex, const float n[]) = 0;
  /** @} */

  /**
   * \name Construction, Destruction, Instantiation
   * @{
   */
public:
  static tledDeformableMembraneContactSurfaceCPU* CreateSurface(const tledMesh &volumeMesh, const tledSurface &membraneMesh);
  static tledDeformableMembraneContactSurfaceCPU* CreateSurface(const XMLNode xmlRep);
  static tledDeformableMembraneContactSurfaceCPU* CreateSurface(const std::string &type);

  tledDeformableMembraneContactSurfaceCPU(void) : tledDeformableContactSurfaceCPU() {}
  virtual ~tledDeformableMembraneContactSurfaceCPU(void) {}
  /** @} */
};

/**
 * \brief Membrane contact surface implementation
 * \ingroup contact
 */
template <class TBaseSurface>
class tledDeformableMembraneContactSurfaceImplCPU : public tledDeformableContactSurfaceImplCPU<TBaseSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledDeformableContactSurfaceImplCPU<TBaseSurface> Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name Normals
   * @{
   */
public:
  virtual const float* GetNodeNormalCached(const int nodeIndex); 
  virtual const float* GetFacetNormalCached(const int facetIndex);
  /** @} */

  /**
   * \name Membrane Geometry Definition
   * @{
   */
public:
  float GetFacetThickness(const int facetIndex) const;
  /** @} */

  /**
   * \name Pre-Computation/Geometry Update Routines
   * @{
   */
protected:
  /** Computes the contact surface node positions from the corresponding membrane nodes, normal, and thickness information */
  void OffsetNodes(void);

public:
  virtual void UpdateMembraneFacetNormal(const int membraneFacetIndex, const float n[]);

  virtual void Init(void);
  virtual void Update(const float us[]); 
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledDeformableMembraneContactSurfaceImplCPU(const tledMesh &volMesh, const tledSurface &membrane) { this->ConstructFromSolidMeshAndMembrane(volMesh, membrane); }
  tledDeformableMembraneContactSurfaceImplCPU(void) : Superclass() {}
  virtual ~tledDeformableMembraneContactSurfaceImplCPU(void) {}
  /** @} */
};

/**
 * \brief Shorthand for deformable triangle surface meshes with a membrane component.
 * \ingroup contact
 */
typedef tledDeformableMembraneContactSurfaceImplCPU<tledDeformableMembraneContactSurfaceImpl<3, tledDeformableMembraneContactSurfaceCPU> > tledDeformableMembraneContactSurfaceT3CPU;

/**
 * \brief Shorthand for deformable quad surface meshes with a membrane component.
 * \ingroup contact
 */
typedef tledDeformableMembraneContactSurfaceImplCPU<tledDeformableMembraneContactSurfaceImpl<4, tledDeformableMembraneContactSurfaceCPU> > tledDeformableMembraneContactSurfaceQ4CPU;

#include "tledDeformableMembraneContactSurfaceCPU.tpp"
#endif
