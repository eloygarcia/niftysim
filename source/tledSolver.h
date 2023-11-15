// =========================================================================
// File:       tledSolver.h
// Purpose:    Main finite element object base class
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    July 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================
#ifndef tledSolver_H
#define tledSolver_H

#include "tledMesh.h"
#include "tledModel.h"
#include "tledContactManager.h"

#include <vector>
#include <string>
#include <cstdio>

/**
 * \defgroup solver Solver
 */

class tledShellSolver;

/**
 * \brief Solver base class
 * \ingroup solver
 */
class tledSolver
{
private:
  tledShellSolver *mp_ShellSolver;
  tledContactManager* Contacts;

protected:
  /** Setter for the shell solver object. The object is assumed to reside in heap memory, which is deallocated on destruction. */
  void SetShellSolver(tledShellSolver *p_solver) { mp_ShellSolver = p_solver; }

public:
   tledSolver() : mp_ShellSolver(NULL) {}
  virtual ~tledSolver();
   
   virtual void Init(tledModel* Model) = 0;
  
  /**
   * \name Boundary Conditions
   * @{
   */
public:
   virtual void SetFixed(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ) = 0;
   virtual void SetDisps(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ) = 0;
   virtual void SetDisps(std::vector<int>* IndX, std::vector<float>* UX, std::vector<int>* IndY, std::vector<float>* UY,
            std::vector<int>* IndZ, std::vector<float>* UZ) = 0;
   virtual void SetDisps(std::vector<float>* UX, std::vector<float>* UY, std::vector<float>* UZ) = 0;

   virtual void SetExtForces(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ) = 0;
   virtual void SetExtForces(std::vector<int>* IndX, std::vector<float>* FX, std::vector<int>* IndY, std::vector<float>* FY, std::vector<int>* IndZ, std::vector<float>* FZ) = 0;
   virtual void SetExtForces(std::vector<float>* FX, std::vector<float>* FY, std::vector<float>* FZ) = 0;

  /** Copies the current effective external forces to a buffer (must be DOF x Number of Nodes in size), mainly intended for unit testing. */
  virtual float* GetAllExtForces(float *p_dst) const = 0;
  /** @} */

   virtual void PerformStep() = 0;
   virtual void GetDivergence(bool* Div) = 0;
   /** Unsets the divergence flag */
   virtual void UnsetDivergence(void) = 0;
   virtual void GetForces(std::vector<int>* NodeInd, std::vector<float>* Forces) = 0;
   virtual float* GetAllForces() = 0;
   virtual void GetDisps(std::vector<int>* NodeInd, std::vector<float>* Disps) = 0;
   virtual float* GetAllDisps() = 0;
   virtual float* GetAllNextDisps() = 0;
   virtual void GetNodeVMStress(float* NodeSVM) = 0; // Get nodal averaged Von Mises equivalent stresses
   virtual void GetNodeVMStrain(float* NodeEVM) = 0; // Get nodal averaged Von Mises equivalent strains
   virtual void GetGaussPtSPKStress(float* S) = 0; // Get 2nd Piola-Kirchhoff stresses for each Gauss Pt
   virtual void GetGaussPtGreenStrain(float* E) = 0; // Get Green's strain for each Gauss Pt
   virtual tledMesh* GetMesh() = 0;
   virtual const tledMesh* GetMesh() const = 0;
   virtual void GetMassVector(float* mass) const = 0;
   virtual void PrepareOutput(void) = 0;
   virtual void PrintNodalForces() = 0;
   virtual void PrintDispNodalForces() = 0;
   virtual void PrintDispNodalForceSums() = 0;
   virtual void PrintKineticEnergy();
   virtual void PrintStrainEnergy();
   virtual void PrintNodalDisps() = 0;
   virtual void PrintNodalForceSums() = 0;
   
  /**
   * \name Contacts
   * @{
   */
public:
   virtual void SetContactManager(tledContactManager* contacts) { Contacts = contacts; }
   tledContactManager& GetContactManager(void) { return *Contacts; }
   const tledContactManager& GetContactManager(void) const { return *Contacts; }
  /** @} */
   
   virtual void InitialiseSolutionVariables(void) = 0;
   virtual void InitialiseConstraints(void) = 0;
   virtual void SetTimeStep(double dt) = 0;
   virtual void ComputeStresses(void) = 0;
   
   virtual float GetStrainEnergy(void) = 0;
   virtual float GetKineticEnergy(void) = 0;
   virtual void SetAllDisps(float* U) = 0;
   
   virtual void SetElementMatParams(int el, std::vector<float> params) = 0;
   virtual void SetMultipleElementMatParams(std::vector<int> el, std::vector<float>* params) = 0;
   
   virtual void SetGeometry(std::vector<float> NodeCoords) = 0; // Set new node coords (NB: topology remains unchanged)

  /**
   * \name Shell and Membrane Elements
   * @{
   */
public:
  /** Query if there is a membrane component to the simulation */
  bool HasMembrane(void) const { return mp_ShellSolver != NULL; }

  /** \brief R/W access to shell/membrane solver backend. */
  tledShellSolver& GetShellSolver(void) { return *mp_ShellSolver; }
  /** \brief R/O access to shell/membrane solver backend. */
  const tledShellSolver& GetShellSolver(void) const { return *mp_ShellSolver; }
  /** @} */
};

#endif	// tledSolver_H
