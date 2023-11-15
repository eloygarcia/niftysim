// =========================================================================
// File:       tledSurfaceConstraint.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSurfaceConstraint_H
#define tledSurfaceConstraint_H

#include "tledSurface.h"
#include "tledHelper.h"
#include "tledSolver.h"

#include <vector>

/**
 * \brief Base class for surface mesh-based constraints (traction, pressure)
 * \ingroup constraints
 */
class tledSurfaceConstraint : public tledConstraint {
  /**
   * \name Geometric definition
   * @{
   */
protected:
  class ConstraintSurface;
  class QuadConstraintSurface;
  class TriangleConstraintSurface;

private:
  tledSolver *mp_CurrentConfSolver;
  ConstraintSurface *mp_Surface;
  int m_LastUpdateStep;

protected:
  /** Last step in which the surface was updated to current configuration */
  int GetLastUpdateStep(void) const { return m_LastUpdateStep; }

  /** Updates face geometry to current configuration */
  void UpdateGeometry(const int step);

  const ConstraintSurface& GetSurface(void) const { return *mp_Surface; }
  ConstraintSurface& GetSurface(void) { return *mp_Surface; }

  const tledSolver& GetSolver(void) const { return *mp_CurrentConfSolver; }  
  tledSolver& GetSolver(void) { return *mp_CurrentConfSolver; }  
  bool HaveSolver(void) const { return mp_CurrentConfSolver != NULL; }

  /** Factory member for surface object, called by SetFaces. */
  virtual ConstraintSurface* CreateSurface(const int type, const std::vector<int> &faces);

public:
  /** 
   * \brief Configures the constraint geometry.
   *
   * Overloading sub-classes <i>must</i> call this function. If the intent is to add new surface types, the CreateSurface factory member function
   * should be overloaded instead.
   * \param type 1 = triangle, 0 = quadrilateral 
   */
  virtual void SetFaces(const int type, const std::vector<int> &faces);

  int GetNumberOfFaceVertices(void) const;
  int GetNumberOfFaces(void) const;
  int GetFaceType(void) const;

  /** 
   * \brief Set solver object from which to retrieve current-configuration information. 
   *
   * Needs R/W access for GPU solvers!
   */
  void SetSolver(tledSolver &r_solver) { mp_CurrentConfSolver = &r_solver; }
  /** @} */

  /**
   * \name Constraint Manager Interface
   * @{
   */
private:
  std::vector<float> m_ConstraintForces;

protected:
  /** Can be used as a force holding buffer for GetForcesVal */
  std::vector<float>& GetConstraintForcesBuffer(void) { return m_ConstraintForces; }
  
  /** Extracts  for all nodes contained in the boundary the component corresponding to the given DOF index from an accumulation buffer (size 3x#mesh nodes). */ 
  std::vector<float>* ConvertToConstraintDOFForces(const std::vector<float>::const_iterator allForcesBegin, const std::vector<float>::const_iterator allForcesEnd, const int dof);

public:
  virtual std::vector<int>& GetInd(void);
  virtual const std::vector<int>& GetInd(void) const;

  virtual std::vector<int>* GetForceInd(int dof);
  virtual std::vector<float>* GetForceVal(int dof, int step, double dt, double T) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledSurfaceConstraint(void);
  virtual ~tledSurfaceConstraint(void);
  /** @} */
};

class tledSurfaceConstraint::ConstraintSurface : public tledSurface {
  /**
   * \name Surface Configuration
   * @{
   */
private:
  std::vector<int> m_NodeIndices;
  std::vector<float> m_NodeCoordinates;

protected:
  const float* GetNodeVector(void) const { return &m_NodeCoordinates.front(); }
  int GetNodeVectorSize(void) const { return int(m_NodeCoordinates.size()); }

public:
  /** R/W reference to node list that can be returned to constraint manager */
  std::vector<int>& GetNodeIndices(void) { return m_NodeIndices; }
  const std::vector<int>& GetNodeIndices(void) const { return m_NodeIndices; }

  virtual void Init(const std::vector<int> &facetDefs, const int numMeshNodes);
  void Update(const float X0[], const float U[]);
  /** @} */

  /**
   * \name Geometric Queries
   *
   * Mostly just wrappers for implementations in tledSurfaceImpl. May be more convenient, but if performance is critical use 
   * corresponding tledSurfaceImpl member functions.
   * @{
   */
public:
  /** Unnormalised normal, s.t. \f$||n|| = 2\cdot\mbox{area}\f$ */ 
  virtual float* ComputeConstraintFacetNormal(float *p_dst, const int facetIndex) const = 0;
  virtual float ComputeConstraintFacetArea(const int facetIndex) const = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  ConstraintSurface(void) {}
  virtual ~ConstraintSurface(void) {}
  /** @} */
};

class tledSurfaceConstraint::QuadConstraintSurface : public tledSurfaceImpl<tledBasicSurfaceFacet<4>, ConstraintSurface> {
  /**
   * \name Surface Configuration
   * @{
   */
public:
  virtual void Init(const std::vector<int> &facetDefs, const int numMeshNodes);
  /** @} */

  /**
   * \name Geometric Queries
   * @{
   */
public:
  virtual float* ComputeConstraintFacetNormal(float *p_dst, const int facetIndex) const { return this->ComputeFacetNormal(p_dst, facetIndex); }
  virtual float ComputeConstraintFacetArea(const int facetIndex) const { return this->ComputeFacetArea(facetIndex); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  QuadConstraintSurface(void) {}
  virtual ~QuadConstraintSurface(void) {}
  /** @} */  
};

class tledSurfaceConstraint::TriangleConstraintSurface : public tledSurfaceImpl<tledBasicSurfaceFacet<3>, ConstraintSurface> {
  /**
   * \name Surface Configuration
   * @{
   */
public:
  virtual void Init(const std::vector<int> &facetDefs, const int numMeshNodes);
  /** @} */

  /**
   * \name Geometric Queries
   * @{
   */
public:
  virtual float* ComputeConstraintFacetNormal(float *p_dst, const int facetIndex) const { return this->ComputeFacetNormal(p_dst, facetIndex); }
  virtual float ComputeConstraintFacetArea(const int facetIndex) const { return this->ComputeFacetArea(facetIndex); }  
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  TriangleConstraintSurface(void) {}
  virtual ~TriangleConstraintSurface(void) {}
  /** @} */  
};

#endif
