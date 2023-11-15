// =========================================================================
// File:       tledSolverGPU.h
// Purpose:    Main finite element object - GPU version
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    April 2011
// 
// Copyright (c) 2011, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifdef _GPU_

#ifndef tledSolverGPU_H
#define tledSolverGPU_H

#include "tledSolver.h"
#include "tledTimeStepperGPU.h"
#include "tledCUDAHelpers.h"

/**
 * \brief Main TLED solver with GPU execution
 * \ingroup solver
 */
class tledSolverGPU : public tledSolver
{
public:
   tledSolverGPU();
   virtual ~tledSolverGPU();
   
   virtual void Init(tledModel* Model);
   virtual void SetFixed(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ);
   virtual void SetDisps(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ) {;}
   virtual void SetDisps(std::vector<int>* IndX, std::vector<float>* UX, std::vector<int>* IndY, std::vector<float>* UY,
                         std::vector<int>* IndZ, std::vector<float>* UZ);
   virtual void SetDisps(std::vector<float>* UX, std::vector<float>* UY, std::vector<float>* UZ) {;}
   virtual void SetExtForces(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ) {;}
   virtual void SetExtForces(std::vector<int>* IndX, std::vector<float>* FX, std::vector<int>* IndY, std::vector<float>* FY,
                             std::vector<int>* IndZ, std::vector<float>* FZ);
   virtual void SetExtForces(std::vector<float>* FX, std::vector<float>* FY, std::vector<float>* FZ) {;}
   virtual float* GetAllExtForces(float *p_dst) const;
   virtual void PerformStep();
   virtual void GetDivergence(bool* Div);
   virtual void UnsetDivergence(void);
   virtual void GetForces(std::vector<int>* NodeInd, std::vector<float>* Forces);
   virtual float* GetAllForces();
   virtual void GetDisps(std::vector<int>* NodeInd, std::vector<float>* Disps);
   virtual float* GetAllDisps();
   virtual float* GetAllNextDisps();
   virtual void GetNodeVMStress(float* NodeSVM);
   virtual void GetNodeVMStrain(float* NodeEVM);
   virtual void GetGaussPtSPKStress(float* S) {;}
   virtual void GetGaussPtGreenStrain(float* E) {;}
   virtual tledMesh* GetMesh() {return Mesh;}
   virtual const tledMesh* GetMesh() const {return Mesh;}
   virtual void GetMassVector(float* mass) const;
   virtual void PrepareOutput();
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

  /** Note: invalid with viscoelastic materials, since the deformation history is not taken into account */
  virtual void ComputeStresses(void); 
   
   virtual void SetElementMatParams(int el, std::vector<float> params);
   virtual void SetMultipleElementMatParams(std::vector<int> el, std::vector<float>* params);
   
   virtual void SetGeometry(std::vector<float> NodeCoords);

  /**
   * \name Element Parameters
   * @{
   */
private:
   std::vector<int> m_ElementElSetIndex; 
   std::vector<float> m_ElSetDensities; 

protected:
  int GetElementElementSetIndex(const int elementIndex) const { return m_ElementElSetIndex[elementIndex]; }
  float GetElementSetDensity(const int elementSetIndex) const { return m_ElSetDensities[elementSetIndex]; }
  float GetElementDensity(const int elementIndex) const { return this->GetElementSetDensity(this->GetElementElementSetIndex(elementIndex)); }
  /** @} */

  /**
   * \name Access to GPU on-device memory
   * @{
   */
public:
  const float4* GetAllOnDeviceCurrentDisplacements(void) const;
  const float4* GetAllOnDeviceNextDisplacements(void) const;
  /** @} */
   
protected:
   // FUNCTIONS
   void InstantiateTimeStepper(void);
   void ComputeElementT4Variables(tledModel* Model, std::vector<std::vector<int2> >* NInd);
   void ComputeElementT4ANPVariables(tledModel* Model, std::vector<std::vector<int2> >* NInd);
   void ComputeElementH8Variables(tledModel* Model, std::vector<std::vector<int2> >* NInd);
   void ComputeLEparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa); // Linear elastic
   void ComputeNHparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa); // Neo-Hookean
   void ComputeTIparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa); // Transversely isotropic
   void ComputeABparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa); // Arruda-Boyce
   void ComputePYparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa); // Polynomial
   void ComputeViscparams(float* inMatViscParams, int Ni, int Nv, int2* N, float2* Ai, float2* Av); // Viscoelastic
   double ComputeHexVol(double x[8][3]);
   void ComputeH8DhDx(double x[8][3], double V, double fDhDx[8][3]);
   void ComputeBmat(double x[8][3], double B[8][3]);
   void ComputeCIJK(double C[8][8][8]);
   void ComputeNodalVolumes();
   void UpdateMaterialParams(int el, std::vector<float> params);
   void UpdateHGParams(int el);
   void UpdateElementT4Geometry();
   void UpdateElementH8Geometry();
   
   // HOST VARIABLES
   float4* h_DhDx;      // Shape function derivatives (global)
   int4* h_EInd;	// Element node indices
   float4* h_CD;	// Central difference coefficients
   float3* h_F;		// Global nodal forces
   int2* h_NodeMap;	// Lookup array for FCds
   int2* h_FCds;	// Lookup array for FEl
   float4* h_HG;	// Hourglass parameters
   float4* h_Uload;	// Displaced DOFs (defined as float4 from the outset since it will point to pinned memory)
   float4* h_R;		// External nodal forces
   float4* h_Vol_MType_K;	// Element volumes (.x), material types (.y), and bulk moduli (.z)
   float4* h_ElasticParams;	// Element material params
   int2* h_NumProny;		// Number of isochoric and volumetric Prony series terms
   float2* h_bkwdEulerIso;	// Backward Euler integration coefficients - Iso
   float2* h_bkwdEulerVol;	// Backward Euler integration coefficients - Vol
   float3* h_SPKa;         // Second Piola-Kirchhoff stress tensor [S11,S22,S33]
   float3* h_SPKb;         // Second Piola-Kirchhoff stress tensor [S12,S23,S13]
   float3* h_Ca;           // Right Cauchy-Green deformation tensor [C11,C22,C33]
   float3* h_Cb;           // Right Cauchy-Green deformation tensor [C12,C23,C13]
   float* h_Pa;		// Nodal pressures (ANP)
   float* h_Va;		// Nodal volumes (ANP)
   bool h_Divergence;

   // DEVICE VARIABLES
   float4* d_DhDx;
   int4* d_EInd;
   float4* d_CD;
   float4* d_HG;
   float4* d_F;
   float4* d_FEl;
   float4* d_Uload;
   float4* d_R;
   int2* d_FCds;
   int2* d_NodeMap;
   int4* d_FixMask;
   int4* d_DispMask;
   float4* d_Vol_MType_K;
   float4* d_ElasticParams;
   int2* d_NumProny;
   float2* d_bkwdEulerIso;
   float2* d_bkwdEulerVol;
   float4* d_StressStateIso;	// State variables for viscoelastic materials
   float4* d_StressStateVol;
   float4* d_SPKa;
   float4* d_SPKb;
   float4* d_Ca;
   float4* d_Cb;
   float* d_Pa;
   float* d_Va;
   bool* d_Divergence;
   float4* d_NodeCds;
   tledGPUContacts* d_Contacts;
   
   // UTILITY VARIABLES
   tledMesh* Mesh;         // Geometric information
   const char* EType;	// Element type
   int NumNodes;           // Number of nodes
   int NumEls;		// Number of elements
   int NDOF;		// Number of DOFs per node
   double Dt;		// Time step
   double Ttot;            // Total analysis time
   float alpha;            // Damping coefficient
   float* M;		// Mass vector
   float* ElementVols;  // Element volumes (note: may be different from values in h_Vol_MType_K.x for H8 elements)
   int FCdsLength;         // Length of FCds array
   int numElasticParamsF4; // Number of float4 variables required to store each el's elastic mat params
   int2 maxNumViscTerms;   // Max number of isochoric and volumetric visco terms for any element in the mesh
   bool ANP;		// Flag to indicate use of ANP formulation
   float2* HGLame;      // Element Lame params, [Lambda,Mu], required for hourglass control
   float HGKappa; // Hourglass control parameter
   float* SVM;    // List of element Von Mises effective stresses
   float* EVM;    // List of element Von Mises effective strains
  tledTimeStepperGPU *mp_TimeStepper;
};



#endif // tledSolverGPU_H

#endif // _GPU_
