// =========================================================================
// File:       tledDeformableContactSurfaceCPU.h
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
#ifndef tledDeformableContactSurfaceCPU_H
#define tledDeformableContactSurfaceCPU_H

#include "tledContactSurfaceCPU.h"
#include "tledDeformableContactSurface.h"
#include "tledVectorArithmetic.h"

#include <algorithm>

/**
 * \brief Deformable contact surface (CPU) interface.
 * \ingroup contact
 */
class tledDeformableContactSurfaceCPU : public tledDeformableContactSurface {
  /**
   * \name Pre-Computation/Geometry Update Routines
   * @{
   */
public:
  /**
   * Updates the surface configuration using the given displacement vector.
   * Surface normals are not updated, and computed ad-hoc in computeContactForces.
   */
  virtual void Update(const float us[]) { tledDeformableContactSurface::IncUpdateCounter(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledDeformableContactSurfaceCPU(void) {}

public:
  static tledDeformableContactSurfaceCPU* CreateSurface(const tledMesh &mesh);
  static tledDeformableContactSurfaceCPU* CreateSurface(const XMLNode xmlRep);
  static tledDeformableContactSurfaceCPU* CreateSurface(const std::string &type);

  virtual ~tledDeformableContactSurfaceCPU(void) {}
  /** @} */
};

/**
 * \brief CPU deformable contact surface  implementation.
 * \ingroup contact
 *
 * TBaseSurface must extend tledDeformableContactSurfaceCPU, tledDeformableContactSurfaceImpl.
 */
template <class TBaseSurface>
class tledDeformableContactSurfaceImplCPU : public tledContactSurfaceImplCPU<TBaseSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSurfaceImplCPU<TBaseSurface> Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */
  
  /**
   * \name Geometry Update
   * @{
   */
protected:
  typedef typename Superclass::CachedQuantity CachedQuantity;
  typedef typename Superclass::CachedNormal CachedNormal;
  typedef typename Superclass::NodeNormalData NodeNormalData;

public:
  virtual void Update(const float us[]);
  virtual void Init(void);

  /**
   * \brief Stores the surface after it has been ensured there are no intersections.
   *
   * Stored node positions are accessible through GetOldNodeCoordinates.
   */
  virtual void Save(void);
  /** @} */

  /**
   * \name Nodes
   * @{
   */
private:
  float const *mpc_OldNodeCoordinates;
  std::vector<float> m_OldNodeCoordinates;

public:  
  /** @{ */
  /** Previous time-stpe node positions */
  const float* GetAllOldNodeCoordinates(void) const { return mpc_OldNodeCoordinates; }
  const float* GetOldNodeCoordinates(const int nodeIndex) const { return mpc_OldNodeCoordinates + 3*nodeIndex; }
  void SetAllOldNodeCoordinates(const float *nodeCds);
  /** @} */

  /**
   * \brief Computes the <i>normalised</i> normal of a surface node.
   *
   * Facet normals required in the computation of the node normal are updated automatically.
   */
  virtual const float* GetNodeNormalCached(const int nodeIndex); 

  /**
   * \brief Computes the <i>normalised</i> normal of a surface node.
   *
   * Facet normals required in the computation of the node normal are updated automatically.
   */
  const float* GetNodeNormal(const int nodeIndex) const;

  /**
   * \brief Computes a continuous normal from vertex normals.
   *
   * Updates node normals if required.
   */
  float* ComputeC0NormalCached(float *p_dstNormal, const int facet[], const float shapeValues[]);  
  /** @} */

  /**
   * \name Facets
   * @{
   */
public:
  struct CachedProjectionOperator : public CachedQuantity {
    float Operator[12];
  };

protected:
  struct FacetData {
    CachedNormal Normal, OldNormal;
    CachedProjectionOperator ProjectionOperator;
  };

private:
  std::vector<FacetData> m_FacetData;

protected:
  CachedNormal& GetFacetNormalCache(const int facetInd) { 
    FacetData &r_data = m_FacetData[facetInd]; /* Do not remove this reference, required due to NVCC compiler bug (error: __T11 undefined etc.)! */
    return r_data.Normal; 
  }

  CachedNormal& GetOldFacetNormalCache(const int facetInd) { 
    FacetData &r_data = m_FacetData[facetInd];
    return r_data.OldNormal; 
  }

  const CachedNormal& GetOldFacetNormalCache(const int facetInd) const { 
    const FacetData &data = m_FacetData[facetInd];
    return data.OldNormal; 
  }

  CachedProjectionOperator& GetProjectionOperatorCache(const int facetInd) { 
    FacetData &r_data = m_FacetData[facetInd];
    return r_data.ProjectionOperator; 
  }

  const CachedProjectionOperator& GetProjectionOperatorCache(const int facetInd) const { 
    const FacetData &data = m_FacetData[facetInd];
    return data.ProjectionOperator; 
  }

public:
  /**
   * Returns a pointer to the unnormalised facet normal, computes a facet normal if necessary.
   */
  virtual const float* GetFacetNormalCached(const int facetIndex);

  /**
   * Computes a normalised facet normal by calling GetUnnormalisedFacetNormalCached.
   */
  float* ComputeNormalisedFacetNormalCached(float *p_n, const int facetInd);

  /**
   * \brief Returns the normalised normal on a stored, old surface facet, computes it if necessary
   */
  const float* GetNormalisedOldFacetNormalCached(const int facetInd);

  /**
   * \brief Returns the normalised normal on a stored, old surface facet, assumes (without checking) that the normal is up-to-date.
   */
  const float* GetNormalisedOldFacetNormal(const int facetInd) const;

  virtual void SetNumberOfFacets(const int numFacets);

  /** Returns the facet projection operator, recomputes it if necessary */
  const float* GetFacetProjectionOperatorCached(const int facetIndex);

  /** Returns the facet projection operator, assumes it's up-to-date */
  const float* GetFacetProjectionOperator(const int facetIndex) const;

  bool ProjectOntoFacetCached(float *p_xi, const float x[], const int facetIndex) { return Superclass::ProjectOntoFacet(p_xi, x, GetFacetProjectionOperatorCached(facetIndex)); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void LoadFromXMLPostloadHook(void);

  tledDeformableContactSurfaceImplCPU(const tledMesh &volMesh) { this->ConstructFromSolidMesh(volMesh); }
  tledDeformableContactSurfaceImplCPU(void) {}
  virtual ~tledDeformableContactSurfaceImplCPU(void) {}
  /** @} */
};

/**
 * \brief Shorthand for deformable triangle surface meshes.
 * \ingroup contact
 */
typedef tledDeformableContactSurfaceImplCPU<tledDeformableContactSurfaceImpl<3, tledDeformableContactSurfaceCPU> > tledDeformableContactSurfaceT3CPU;

/**
 * \brief Shorthand for deformable quad surface meshes.
 * \ingroup contact
 */
typedef tledDeformableContactSurfaceImplCPU<tledDeformableContactSurfaceImpl<4, tledDeformableContactSurfaceCPU> > tledDeformableContactSurfaceQ4CPU;

#include "tledDeformableContactSurfaceCPU.tpp"
#endif
