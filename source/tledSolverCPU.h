// =========================================================================
// File:       tledSolverCPU.h
// Purpose:    Main finite element object - CPU version
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledSolverCPU_H
#define tledSolverCPU_H

#include "tledSolver.h"
#include "tledMaterial.h"
#include "tledElement.h"
#include "tledElementT4.h"
#include "tledElementH8.h"
#include "tledElementT4ANP.h"
#include "tledTimeStepperCPU.h"

/**
 * \brief Main TLED solver with CPU execution
 * \ingroup solver
 */
class tledSolverCPU : public tledSolver
{
public:
   tledSolverCPU();
   virtual ~tledSolverCPU();

   virtual void SetFixed(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ);
   virtual void SetDisps(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ);
   virtual void SetDisps(std::vector<int>* IndX, std::vector<float>* UX, std::vector<int>* IndY, std::vector<float>* UY,
            std::vector<int>* IndZ, std::vector<float>* UZ);
   virtual void SetDisps(std::vector<float>* UX, std::vector<float>* UY, std::vector<float>* UZ);
   virtual void SetExtForces(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ);
   virtual void SetExtForces(std::vector<int>* IndX, std::vector<float>* FX, std::vector<int>* IndY, std::vector<float>* FY,
            std::vector<int>* IndZ, std::vector<float>* FZ);
   virtual void SetExtForces(std::vector<float>* FX, std::vector<float>* FY, std::vector<float>* FZ);
   virtual float* GetAllExtForces(float *p_dst) const;
   virtual void PerformStep();
   virtual void GetDivergence(bool* Div);
   virtual void UnsetDivergence(void);
   virtual void GetForces(std::vector<int>* NodeInd, std::vector<float>* Forces);
   virtual float* GetAllForces() {return F;}
   virtual void GetDisps(std::vector<int>* NodeInd, std::vector<float>* Disps);
   virtual float* GetAllDisps() {return mp_TimeStepper->GetCurrentDisplacements();} 
   virtual float* GetAllNextDisps() {return mp_TimeStepper->GetNextDisplacements();}
   virtual void GetNodeVMStress(float* NodeSVM);
   virtual void GetNodeVMStrain(float* NodeEVM);
   virtual void GetGaussPtSPKStress(float* S);
   virtual void GetGaussPtGreenStrain(float* E);
   virtual tledMesh* GetMesh() {return Mesh;}
   virtual const tledMesh* GetMesh() const {return Mesh;}
   virtual void GetMassVector(float* mass) const;
   virtual void PrepareOutput() {}
   virtual void PrintNodalForces();
   virtual void PrintDispNodalForces();
   virtual void PrintDispNodalForceSums();
   virtual void PrintNodalDisps();
   virtual void PrintNodalForceSums();   

  virtual void SetContactManager(tledContactManager* contacts);

   virtual void InitialiseSolutionVariables(void);
   virtual void InitialiseConstraints(void);
   virtual void SetTimeStep(double dt);
   
   virtual float GetStrainEnergy(void);
   virtual float GetKineticEnergy(void);
   virtual void SetAllDisps(float* U);
   virtual void ComputeStresses(void); // Note: invalid with viscoelastic materials, since the deformation history is not taken into account
   
   virtual void SetElementMatParams(int el, std::vector<float> params);
   virtual void SetMultipleElementMatParams(std::vector<int> el, std::vector<float>* params);
   
   virtual void SetGeometry(std::vector<float> NodeCoords);

   tledTimeStepperCPU& GetTimeStepper(void) { return *mp_TimeStepper; }
   const tledTimeStepperCPU& GetTimeStepper(void) const { return *mp_TimeStepper; }


  /**
   * \name Misc. Simulation Setup Queries
   * @{
   */
protected:
  const double& GetAlpha(void) const { return m_Alpha; }
  /** Direct (R/O) access to the internal diagonal mass matrix. */
  const float* GetMassInternal(void) const { return M; }

public:
  int GetNumberOfDOFs(void) const { return NDOF; }
  int GetNumberOfNodes(void) const { return NumNodes; }
  int GetNumberOfElements(void) const { return NumEls; }
  const double& GetDt(void) const { return Dt; }
  /** @} */

  /** 
   * \name Time-stepping sub-routines
   * @{
   */
protected:
  void ApplyFixed();
   void ApplyDisps();
   void ApplyExtForces();
   virtual void ComputeNewForces();
   virtual void ComputeNewForcesANP();
   void ComputeNewDisps();
   /** @} */

   /**
    * \name Setup
    * @{
    */
protected:
   void CompileMass(void);

   /** Solver object takes ownership of the tledTimeStepper object, assumes it was dynamically allocated and destroys it upon its own destruction. */
   void SetTimeStepper(tledTimeStepperCPU *p_stepper) { mp_TimeStepper = p_stepper; }

   /** Time-ODE solver setup function */
  virtual void InstantiateTimeStepper(void);

public:
  virtual void Init(tledModel* Model);
  /** @} */

  /**
   * \name Access to Computation Resources
   * @{
   */
protected:
  tledElement** GetElements(void) { return Elements; }
  float* GetInternalForceBuffer(void) { return F; }
  float* GetExternalForceBuffer(void) { return R; }
  float* GetPressureBuffer(void) { return Pa; }
  const float* GetNodeVolumes(void) const { return Va; }
  /** @} */     

private:
   tledMesh* Mesh;	// Geometric information
   int NumNodes;
   int NumEls;
   tledElement** Elements;	// Array of element objects
   int NDOF;		// Number of DOFs per node
   float* A;		// Central difference coefficients
   float* B;
   float* C;
   float* mp_EffectiveF;
   tledTimeStepperCPU *mp_TimeStepper;
   float* F;		// Global nodal forces
   float* R;		// External nodal loads
   float* M;		// Nodal mass vector (corresponding to diagonalised mass matrix)
   bool Divergence;	// Flag to indicate divergence state
   float* Va;		// Nodal volumes
   float* Pa;		// Nodal pressures
   bool ANP;		// Flag to indicate use of ANP formulation
   double Dt;           // Simulation time step
   double m_Alpha;        // Structural damping coefficient

   std::vector<int>* IndFX;	// Indices of fixed nodes
   std::vector<int>* IndFY;
   std::vector<int>* IndFZ;
   std::vector<int>* IndDX;	// Indices of displaced nodes
   std::vector<int>* IndDY;
   std::vector<int>* IndDZ;
   std::vector<int>* IndRX;	// Indices of externally loaded nodes
   std::vector<int>* IndRY;
   std::vector<int>* IndRZ;
   std::vector<float>* UDX;	// Displacement vals for displaced nodes
   std::vector<float>* UDY;
   std::vector<float>* UDZ;
   std::vector<float>* RX;	// Force vals for externally loaded nodes
   std::vector<float>* RY;
   std::vector<float>* RZ;
};

#endif //tledSolverCPU_H
