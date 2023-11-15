// =========================================================================
// File:       tledUnstructuredContactManager.h
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
#ifndef tledUnstructuredContactManager_H
#define tledUnstructuredContactManager_H

#include "tledBVH.h"
#include "tledSelfCollisionBVH.h"
#include "tledBVHTraverser.h"
#include "tledBVHTraverserCPU.h"
#include "xmlParser.h"
#include "tledDeformableContactSurface.h"
#include "tledRigidContactSurface.h"
#ifdef _GPU_
#include "tledCUDAHelpers.h"
#endif

#include <vector>

#if defined GPU_GP_CONTACT && !defined _GPU_
#error "Inconsistency: GPU general-purpose contact modelling pipeline is enabled, GPU support isn't. Please reconfigure."
#endif

class tledSolver;
class tledModel;
class tledDeformableRigidContactSolver;
class tledDeformableDeformableContactSolver;

/**
 * \brief Manages resources related to unstructured-mesh contact modelling.
 * \ingroup contact
 */
class tledUnstructuredContactManager {
  /**
   * \name Contact Surfaces
   * @{
   */
private:
  tledDeformableContactSurface *mp_DeformableSurface;
  std::vector<tledRigidContactSurface*> mvp_RigidSurfaces;
  int m_CoordinateHistorySize;

public:
  template <class TRigidSurface>
  const TRigidSurface& GetRigidSurface(const int surfaceIndex) const { return static_cast<const TRigidSurface&>(*mvp_RigidSurfaces[surfaceIndex]); }
  template <class TRigidSurface>
  TRigidSurface& GetRigidSurface(const int surfaceIndex) { return static_cast<TRigidSurface&>(*mvp_RigidSurfaces[surfaceIndex]); }
  int GetNumberOfRigidSurfaces(void) const { return mvp_RigidSurfaces.size(); }
  void AddRigidSurface(const XMLNode &meshSpecRoot);
  void ResetRigidSurfaces(void);

  template <class TDeformableSurface>
  TDeformableSurface& GetDeformableSurface(void) { return static_cast<TDeformableSurface&>(*mp_DeformableSurface); }
  template <class TDeformableSurface>
  const TDeformableSurface& GetDeformableSurface(void) const { return static_cast<const TDeformableSurface&>(*mp_DeformableSurface); }
  /** Builds the deformable contact surface by extracting the surface of the simulation solid mesh. */
  void AddDeformableSurface(const tledMesh &mesh);
  /** Loads a previously exported deformable contact surface */
  void AddDeformableSurface(const XMLNode &rootNode);

  /** Adds a deformable surface consisting of the surface extracted from the simulation's solid mesh and a membrane */
  void AddDeformableMembraneSurface(const tledSurface &membrane, const tledMesh &mesh);
  void AddDeformableMembraneSurface(const XMLNode &rootNode);  

  /** @{ */
  /**
   * \brief Controls how many time steps to go back for "old" node positions on dynamic surface types (default 10).
   */  
  int GetCoordinateHistorySize(void) const { return m_CoordinateHistorySize; }
  void SetCoordinateHistorySize(const int numTSteps) { m_CoordinateHistorySize = numTSteps; }
  /** @} */
  /** @} */

  /**
   * \name Bounding Volume Hierarchies
   * @{
   */
private:
  std::string m_BVType;
  tledSelfCollisionBVH *mp_DeformableBVH;
  std::vector<tledBVH*> mvp_RigidSurfaceBVHs;
  float m_CloseDistance, m_SafetyMargin;
  float m_BVMargin;

private:
  void _CreateBVHs(const bool hasMembrane);

protected:
  /** Sets the BVH used for contact search on the deformable contact surface. Memory must be dynamically allocated and will be released by the manager's destructor. */
  void SetDeformableBVH(tledSelfCollisionBVH *p_bvh) { mp_DeformableBVH = p_bvh; }

  /** Factory function called during initialisation (w/o XML import) */
  virtual tledSelfCollisionBVH* CreateDeformableBVH(void);

  /** Factory function called during initialisation from XML */
  virtual tledSelfCollisionBVH* LoadDeformableBVH(const XMLNode bvhNode);

  /** Sets the BVH used for contact search on a specific rigid contact surface. Memory must be dynamically allocated and will be released by the manager's destructor. */
  void SetRigidBVH(tledBVH *p_bvh, const int rigidSurfaceIndex);

public:  
  const std::string& GetBVType(void) const { return m_BVType; }
  void SetBVType(const std::string &type) { m_BVType = type; }

  /** 
   * \brief Setter for BV safety margin
   * 
   * By default determined from mesh densities.
   */
  void SetBVMargin(const float margin) { m_BVMargin = margin; }

  /** Returns the guaranteed usable margin between geometry BV bounds. */ 
  float GetBVMargin(void) const { return m_BVMargin; }

  tledSelfCollisionBVH& GetDeformableBVH(void) { return *mp_DeformableBVH; }
  const tledSelfCollisionBVH& GetDeformableBVH(void) const { return *mp_DeformableBVH; }

  const tledBVH& GetRigidBVH(const int rigidSurfaceIndex) const { return *mvp_RigidSurfaceBVHs[rigidSurfaceIndex]; }
  /** @} */

  /**
   * \name Solver Parameters
   * @{
   */
private:
  bool m_UseGPU;
  bool m_DoMultiPassContactHandling;

#ifdef _GPU_
private:
  void _InitGPUSetup(void);
#endif

public:
  bool UseGPU(void) const { return m_UseGPU; }

  /** @{ */
  /** True if the user has requested multi-pass contact handling. */
  bool DoMultiPassContactHandling(void) const { return m_DoMultiPassContactHandling; }
  void SetDoMultiPassContactHandling(const bool v) { m_DoMultiPassContactHandling = v; }  
  /** @} */
  /** @} */

  /**
   * \name Dependencies
   * @{
   */
private:
  tledSolver *mp_Solver;
  tledModel *mp_Model;

public:
  tledSolver& GetSolver(void) { return *mp_Solver; }
  const tledSolver& GetSolver(void) const { return *mp_Solver; }
  tledModel& GetModel(void) { return *mp_Model; }
  const tledModel& GetModel(void) const { return *mp_Model; }
  /** @} */

  /**
   * \name Initialisation and Update
   * @{
   */
public:
  /** Updater for BVHs, surfaces, subclassed for different setups */
  class Updater;

private:
  Updater *mp_Updater;

protected:
  /** Factory for the inner-most updater, responsible for updating the deformable-geometry surface object as well as its BVH */
  virtual Updater* CreateRootUpdater(void);

public:
  /** Updates contact surfaces and BVHs */
  void Update(void);

  /** 
   * \brief Add an updater to the update chain.
   *
   * The memory for the updater object should be allocated dynamically, it will be freed up by the manager's destructor.
   */
  void InstallUpdater(Updater *p_updater);
  /** @} */

  /**
   * \name Contact Response Calculation
   * @{
   */
private:
  float m_Dt;
  tledDeformableDeformableContactSolver *mp_DeformableDeformableSolver;
  std::vector<tledDeformableRigidContactSolver*> mvp_DeformableRigidContactSolvers;

protected:
  /** Dynamic allocation expected, deallocation performed by manager. */
  void SetDeformableContactSolver(tledDeformableDeformableContactSolver *p_solver) { mp_DeformableDeformableSolver = p_solver; }

  /** Order in which solvers are added must match order of rigid contact surfaces. Dynamic allocation expected, deallocation performed by manager. */
  void AddRigidContactSolver(tledDeformableRigidContactSolver *p_solver) { mvp_DeformableRigidContactSolvers.push_back(p_solver); }
  
  /** Solver factory: hasMembrane indicates whether contact solvers capable of handling stand-alone membranes are required. Initialisation of solver objects must also be performed by any overriding function. */
  virtual void CreateSolvers(const bool hasMembrane);

public:
  /** \return Solver time step */
  float GetDt(void) const { return m_Dt; }

  /** \brief Distance at which a node is considered dangerously close to a surface, and at which rate constraints kick in. */
  float GetCloseDistance(void) const { return m_CloseDistance; }

  /** \brief Distance at which a node is considered dangerously close to a surface, and at which rate constraints kick in. (default: BV margin/2) */
  void SetCloseDistance(const float critDist) { m_CloseDistance = critDist; }  

  /** \brief A minimum safety margin guaranteed in all contacts */
  float GetSafetyMargin(void) const { return m_SafetyMargin; }
  void SetSafetyMargin(const float eps) { m_SafetyMargin = eps; }

  /** \return true if user wants deformable-deformable (and self-collisions) handled. */
  bool DoDeformableDeformableContactHandling(void) const { return mp_DeformableDeformableSolver != NULL; }  

  bool ComputeDeformableDeformableContactResponses(float *p_F, const float UNext[], const float UCurr[]); 
  bool ComputeDeformableRigidContactResponses(float *p_F, const float UNext[], const float UCurr[]); 

#ifdef GPU_GP_CONTACT
  bool ComputeDeformableDeformableContactResponses(float4 *dp_R, const float4 *dpc_UNext, const float4 *dpc_UCurr);
  bool ComputeDeformableRigidContactResponses(float4 *dp_R, const float4 *dpc_UNext, const float4 *dpc_UCurr); 
#endif

  /** Saves surfaces for next iteration, etc. */
  virtual void FinishContactHandling(void);
  /** @} */

  /**
   * \name Export
   * @{
   */
private:
  class _XMLImporter;
  class _XMLExporter;

public:
  static const char* GetManagerXMLTag(void) { return "UnstructuredContactManager"; }

  /** Export to XML for repeated calculations on identical geometry */
  XMLNode ExportToXML(void) const;
  /** @} */  

  /**
   * \name Construction, Destruction, Initialisation
   * @{
   */
public:
  /** Computes all BVHs etc. from scratch for the given contact geometry and setup represented by the model */
  virtual void Init(tledModel &r_model, tledSolver &r_solver, const bool useGPU);

  tledUnstructuredContactManager(void);

  /** Creates the object and immediately initialises it */
  tledUnstructuredContactManager(tledModel &r_model, tledSolver &r_solver, const bool useGPU);

  virtual ~tledUnstructuredContactManager(void);
  /** @} */
};

class tledUnstructuredContactManager::Updater {
  /**
   * \name Update Chain
   * @{
   */
private:
  Updater *mp_NextUpdater;
  
public:
  void SetNextInChain(Updater *p_next) { mp_NextUpdater = p_next; }
  /** @} */

  /**
   * \name Execution
   * @{ 
   */
public:
  virtual void Update(tledUnstructuredContactManager &r_manager);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  Updater(void) : mp_NextUpdater(NULL) {}
  virtual ~Updater(void);
  /** @} */
};
  
#endif
