// =========================================================================
// File:       tledUnstructuredContactManager.cpp
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
#include "tledUnstructuredContactManager.h"
#include "tledDeformableContactSurfaceCPU.h"
#include "tledRigidContactSurfaceCPU.h"
#include "tledStaticBVH.h"
#include "tledSelfCollisionBVH.h"
#include "tledRigidMotionBVHUpdater.h"
#include "tledSelfCollisionBVHCreator.h"
#include "tledDeformableDeformableContactSolver.h"
#include "tledDeformableRigidContactSolver.h"
#include "tledMovingRigidContactSurfaceCPU.h"
#include "tledModel.h"
#include "tledSolver.h"
#include "tledAABB.h"
#include "tledOBB.h"
#include "tledShellSolver.h"
#include "tledDeformableMembraneContactSurfaceCPU.h"
#include "tledParallel.h"
#include "tledXMLExporter.h"
#include "tledXMLImporter.h"
#include "tledDeformableDeformableContactSolverCPU.h"
#include "tledDeformableRigidContactSolverCPU.h"

#ifdef GPU_GP_CONTACT
#include "tledDeformableContactSurfaceGPU.h"
#include "tledMovingRigidContactSurfaceGPU.h"
#include "tledStaticBVHGPU.h"
#include "tledSolverGPU.h"
#endif

#include <cstdlib>
#include <iostream>
#include <limits>
#include <cassert>
#include <cmath>
#include <typeinfo>

class _BasicUpdater : public tledUnstructuredContactManager::Updater {
public:
  virtual void Update(tledUnstructuredContactManager &r_manager) {
#ifdef GPU_GP_CONTACT
    if (r_manager.UseGPU()) {
      r_manager.GetDeformableSurface<tledDeformableContactSurfaceGPU>().Update(static_cast<const tledSolverGPU&>(r_manager.GetSolver()).GetAllOnDeviceNextDisplacements());
    } else {
#endif
      r_manager.GetDeformableSurface<tledDeformableContactSurfaceCPU>().Update(r_manager.GetSolver().GetAllNextDisps());
#ifdef GPU_GP_CONTACT
    }
#endif
    r_manager.GetDeformableBVH().Update();    
    tledUnstructuredContactManager::Updater::Update(r_manager);
  }

public:
  virtual ~_BasicUpdater(void) {}
};

class _SetupNotSupportedUpdater : public tledUnstructuredContactManager::Updater {
public:
   virtual void Update(tledUnstructuredContactManager &r_manager) {
     tledFatalError("Setup (currently) not supported.");
   }

public:
  virtual ~_SetupNotSupportedUpdater(void) {}
};

class _MembraneContactSurfaceUpdater : public tledUnstructuredContactManager::Updater {
private:
  std::vector<float> m_MembraneElementThicknesses;

public:
   virtual void Update(tledUnstructuredContactManager &r_manager) {
     tledDeformableMembraneContactSurfaceCPU &r_surface = r_manager.GetDeformableSurface<tledDeformableMembraneContactSurfaceCPU>();

     r_manager.GetSolver().GetShellSolver().ComputeAllElementThicknesses(&m_MembraneElementThicknesses.front(), r_manager.GetSolver().GetAllDisps());
     r_surface.SetMembraneFacetThickenesses(&m_MembraneElementThicknesses.front());

     tledUnstructuredContactManager::Updater::Update(r_manager);
   }

public:
  _MembraneContactSurfaceUpdater(tledDeformableMembraneContactSurfaceCPU &r_defSurface, const tledSolver &mainSolver) {
    m_MembraneElementThicknesses.insert(m_MembraneElementThicknesses.end(), mainSolver.GetShellSolver().GetSurface().GetNumberOfFacets(), std::numeric_limits<float>::quiet_NaN());
    mainSolver.GetShellSolver().GetAllInitialElementThicknesses(&m_MembraneElementThicknesses.front());
    r_defSurface.SetMembraneFacetThickenesses(&m_MembraneElementThicknesses.front());
  }

  virtual ~_MembraneContactSurfaceUpdater(void) {}
};

class _MovingRigidSurfaceUpdater : public tledUnstructuredContactManager::Updater {
private:
  tledMovingRigidContactSurface *mp_Surface;
  tledDynamicBVH *mp_BVH;

public:
  virtual void Update(tledUnstructuredContactManager &r_manager) {
    mp_Surface->Update();
    mp_BVH->Update();
    tledUnstructuredContactManager::Updater::Update(r_manager);
  }

public:
  _MovingRigidSurfaceUpdater(tledMovingRigidContactSurface &r_surface, tledDynamicBVH &r_bvh) : mp_Surface(&r_surface), mp_BVH(&r_bvh) {}  
  virtual ~_MovingRigidSurfaceUpdater(void) {}
};

class tledUnstructuredContactManager::_XMLExporter : public tledXMLExporter<tledUnstructuredContactManager> {
public:
  typedef tledXMLExporter<tledUnstructuredContactManager> Superclass;

protected:
  virtual void WriteBody(void);

public:
  virtual const char* GetRootElementName(void) const { return tledUnstructuredContactManager::GetManagerXMLTag(); }

public:
  virtual ~_XMLExporter(void) {}
};

class tledUnstructuredContactManager::_XMLImporter : public tledXMLImporter<tledUnstructuredContactManager> {
public:
  typedef tledXMLImporter<tledUnstructuredContactManager> Superclass;

private:
  const tledModel &mc_Model;
  const bool mc_UseMembrane;
  const bool mc_UseGPU;  

private:
  template <class TBV>
  void _ImportBVHsWithBVType(void);  

public:
  virtual void Import(void);

public:
  _XMLImporter(const tledModel &model, const bool useMembrane, const bool useGPU) : mc_Model(model), mc_UseMembrane(useMembrane), mc_UseGPU(useGPU) {}
  virtual ~_XMLImporter(void) {}
};

void tledUnstructuredContactManager::Updater::Update(tledUnstructuredContactManager &r_manager) {
  if (mp_NextUpdater != NULL) mp_NextUpdater->Update(r_manager);
}

tledUnstructuredContactManager::Updater::~Updater(void) {
  if (mp_NextUpdater != NULL) delete mp_NextUpdater;
}

void tledUnstructuredContactManager::_XMLExporter::WriteBody() {
  this->CreateTextNode("BVType", this->GetInput().GetBVType());
  this->CreateNumericNode("CloseDistance", this->GetInput().GetCloseDistance());
  this->CreateNumericNode("SafetyMargin", this->GetInput().GetSafetyMargin());
  this->CreateNumericNode("BVMargin", this->GetInput().GetBVMargin());
  this->CreateNumericNode("DoMultiPassContactHandling", this->GetInput().DoMultiPassContactHandling());
  
  this->GetRootNode().addChild(this->GetInput().GetDeformableSurface<tledDeformableContactSurface>().ExportToXML());
  this->GetRootNode().addChild(this->GetInput().GetDeformableBVH().ExportToXML());  

  for (int s = 0; s < this->GetInput().GetNumberOfRigidSurfaces(); s++) {
    this->GetRootNode().addChild(this->GetInput().GetRigidBVH(s).ExportToXML());
  }
}

void tledUnstructuredContactManager::_XMLImporter::Import() {
  if (this->GetRootNode().nChildNode("BVType") > 0) this->GetOutput().SetBVType(this->GetUniqueChild("BVType", true).getText());  
  this->GetOutput().SetCloseDistance(this->GetNumericElementValue<float>(this->GetUniqueChild("CloseDistance", true)));
  this->GetOutput().SetSafetyMargin(this->GetNumericElementValue<float>(this->GetUniqueChild("SafetyMargin", true)));
  this->GetOutput().SetBVMargin(this->GetNumericElementValue<float>(this->GetUniqueChild("BVMargin", true)));
  this->GetOutput().SetDoMultiPassContactHandling(this->GetNumericElementValue<bool>(this->GetUniqueChild("DoMultiPassContactHandling", true)));  

  if (mc_UseMembrane) {
    this->GetOutput().AddDeformableMembraneSurface(this->GetUniqueChild(tledDeformableMembraneContactSurfaceXMLExporter<tledDeformableMembraneContactSurfaceT3CPU>::GetDeformableMembraneContactSurfaceXMLTag(), true)); 
  } else {
    this->GetOutput().AddDeformableSurface(this->GetUniqueChild(tledDeformableContactSurfaceXMLExporter<tledDeformableContactSurfaceT3CPU>::GetDeformableContactSurfaceXMLTag(), true)); 
  }
  this->GetOutput().SetDeformableBVH(this->GetOutput().LoadDeformableBVH(this->GetUniqueChild(tledSelfCollisionBVHXMLExporter<tledSelfCollisionBVHImpl<tledDeformableContactSurfaceT3CPU, tledAABB<2> > >().GetRootElementName(), true)));
  this->GetOutput().GetDeformableBVH().GetUpdater().Init();
  if (this->GetOutput().GetDeformableSurface<tledDeformableContactSurface>().GetNumberOfFacets() != this->GetOutput().GetDeformableBVH().GetNumberOfLeafs()) {
    tledLogErrorStream(tledHelper::FatalError() << "Inconsistency in BVH to surface dimensions on deformable contact surface: " << this->GetOutput().GetDeformableSurface<tledDeformableContactSurface>().GetNumberOfFacets() << "/" << this->GetOutput().GetDeformableBVH().GetNumberOfLeafs() << ". Does precomputation file match simulation file?");
  }  

  for (int si = 0; si < mc_Model.GetNumberOfRigidContactSurfaces(); si++) {
    this->GetOutput().AddRigidSurface(mc_Model.GetRigidContactSurfaceDefinition(si));
  }

  if (this->GetOutput().GetNumberOfRigidSurfaces() != this->GetRootNode().nChildNode("BVH")) {
    tledFatalError("Number of precomputed rigid contact surfaces does not match number in simulation file. Does precomputation file match simulation file?");
  }

  for (int s = 0; s < this->GetOutput().GetNumberOfRigidSurfaces(); s++) {
    const tledRigidContactSurface &grSurf = this->GetOutput().GetRigidSurface<tledRigidContactSurface>(s);

    if (grSurf.IsMoving()) {
      if (!this->GetOutput().UseGPU()) {
	tledMovingRigidContactSurfaceCPU &r_movSurf = this->GetOutput().GetRigidSurface<tledMovingRigidContactSurfaceCPU>(s);
	tledDynamicBVH *p_bvh;

	r_movSurf.SetTotalNumberOfSteps(int(std::ceil(this->GetOutput().GetModel().GetTotalTime()/this->GetOutput().GetModel().GetTimeStep())));
	this->GetOutput().SetRigidBVH(p_bvh = tledDynamicBVH::CreateBVH(grSurf, this->GetRootNode().getChildNode("BVH", s), mc_UseGPU), s); 
	this->GetOutput().InstallUpdater(new _MovingRigidSurfaceUpdater(r_movSurf, *p_bvh));      
      } else {
#ifdef GPU_GP_CONTACT
	tledMovingRigidContactSurfaceGPU &r_movSurf = this->GetOutput().GetRigidSurface<tledMovingRigidContactSurfaceGPU>(s);
	tledDynamicBVH *p_bvh;

	r_movSurf.SetTotalNumberOfSteps(int(std::ceil(this->GetOutput().GetModel().GetTotalTime()/this->GetOutput().GetModel().GetTimeStep())));
	this->GetOutput().SetRigidBVH(p_bvh = tledDynamicBVH::CreateBVH(grSurf, this->GetRootNode().getChildNode("BVH", s), mc_UseGPU), s); 
	this->GetOutput().InstallUpdater(new _MovingRigidSurfaceUpdater(r_movSurf, *p_bvh));      	
#else
	tledFatalFeatureNotEnabledError;
#endif
      }
    } else {
      this->GetOutput().SetRigidBVH(tledStaticBVH::CreateBVH(grSurf, this->GetRootNode().getChildNode("BVH", s), mc_UseGPU), s);
    }

    if (this->GetOutput().GetRigidSurface<tledRigidContactSurface>(s).GetNumberOfFacets() != this->GetOutput().GetRigidBVH(s).GetNumberOfLeafs()) {
      tledLogErrorStream(tledHelper::FatalError() << "Inconsistency in BVH to surface dimensions on rigid contact surface " << s << ". Does precomputation file match simulation file?");
    }  
  } /* for rigid surfaces */
}

tledSelfCollisionBVH* tledUnstructuredContactManager::LoadDeformableBVH(const XMLNode xmlNode) {
  return tledSelfCollisionBVH::CreateBVH(this->GetDeformableSurface<tledDeformableContactSurface>(), xmlNode, this->UseGPU());
}

tledSelfCollisionBVH* tledUnstructuredContactManager::CreateDeformableBVH() {
  return tledSelfCollisionBVH::CreateBVH(this->GetDeformableSurface<tledDeformableContactSurface>(), this->GetBVType(), this->GetBVMargin(), this->GetBVMargin() - 1.001f*this->GetCloseDistance(), this->UseGPU());
}

void tledUnstructuredContactManager::_CreateBVHs(const bool hasMembrane) {
  mp_DeformableBVH = this->CreateDeformableBVH();
  for (int sInd = 0; sInd < this->GetNumberOfRigidSurfaces(); sInd++) {
    const tledRigidContactSurface &surface = this->GetRigidSurface<tledRigidContactSurface>(sInd);

    if (surface.IsMoving()) {
      if (this->UseGPU()) {
#ifdef GPU_GP_CONTACT
	tledMovingRigidContactSurfaceGPU &r_movSurf = this->GetRigidSurface<tledMovingRigidContactSurfaceGPU>(sInd);

	r_movSurf.SetTotalNumberOfSteps(int(std::ceil(GetModel().GetTotalTime()/GetModel().GetTimeStep())));
	mvp_RigidSurfaceBVHs.push_back(tledDynamicBVH::CreateBVH(surface, this->GetBVType(), this->GetBVMargin(), this->UseGPU()));
	this->InstallUpdater(new _MovingRigidSurfaceUpdater(r_movSurf, static_cast<tledDynamicBVH&>(*mvp_RigidSurfaceBVHs.back())));
#else
	tledFatalFeatureNotEnabledError;
#endif
      } else {
	tledMovingRigidContactSurfaceCPU &r_movSurf = this->GetRigidSurface<tledMovingRigidContactSurfaceCPU>(sInd);

	r_movSurf.SetTotalNumberOfSteps(int(std::ceil(GetModel().GetTotalTime()/GetModel().GetTimeStep())));
	mvp_RigidSurfaceBVHs.push_back(tledDynamicBVH::CreateBVH(surface, this->GetBVType(), this->GetBVMargin(), this->UseGPU()));
	this->InstallUpdater(new _MovingRigidSurfaceUpdater(r_movSurf, static_cast<tledDynamicBVH&>(*mvp_RigidSurfaceBVHs.back())));
      }
    } else {
      mvp_RigidSurfaceBVHs.push_back(tledStaticBVH::CreateBVH(surface, this->GetBVType(), this->GetBVMargin(), this->UseGPU()));
    }
  }
}

void tledUnstructuredContactManager::CreateSolvers(const bool hasMembrane) {
  if (this->GetModel().DoDeformableDeformableContactHandling()) {
    this->SetDeformableContactSolver(tledDeformableDeformableContactSolver::CreateContactSolver(*this));
    if (this->GetModel().DoSelfCollisionContactHandling() != mp_DeformableDeformableSolver->DoSelfCollision()) mp_DeformableDeformableSolver->ToggleDoSelfCollision();    
  } else this->SetDeformableContactSolver(NULL);
      
  for (int sInd = 0; sInd < this->GetNumberOfRigidSurfaces(); sInd++) {
    this->AddRigidContactSolver(tledDeformableRigidContactSolver::CreateContactSolver(*this, sInd));
  }

  if (mp_DeformableDeformableSolver != NULL) mp_DeformableDeformableSolver->Init();
  for (std::vector<tledDeformableRigidContactSolver*>::iterator ip_s = mvp_DeformableRigidContactSolvers.begin(); ip_s < mvp_DeformableRigidContactSolvers.end(); ip_s++) {
    (*ip_s)->Init();
  }
}

XMLNode tledUnstructuredContactManager::ExportToXML() const {
  _XMLExporter exporter;

  exporter.SetInput(*this);
  exporter.Export();

  return exporter.GetRootNode();
}

tledUnstructuredContactManager::tledUnstructuredContactManager() : m_BVType("Unintialised"), m_UseGPU(false), m_DoMultiPassContactHandling(false), mp_Solver(NULL), mp_Model(NULL) {
  mp_DeformableDeformableSolver = NULL;
  mp_DeformableBVH = NULL;
  m_CoordinateHistorySize = -1;
  m_CloseDistance = m_SafetyMargin = m_BVMargin = std::numeric_limits<float>::quiet_NaN();
  mp_Updater = NULL;
}

tledUnstructuredContactManager::tledUnstructuredContactManager(tledModel &r_model, tledSolver &r_solver, const bool useGPU) {
  this->Init(r_model, r_solver, useGPU);
}

tledUnstructuredContactManager::Updater* tledUnstructuredContactManager::CreateRootUpdater() {
  return new _BasicUpdater();
}

void tledUnstructuredContactManager::Init(tledModel &r_model, tledSolver &r_solver, const bool useGPU) {
  const bool hasMembrane = r_model.GetNumberOfShellElementSets() > 0 && !r_model.DoShellUseMeshSurface();
  const bool doLoadFromXML = r_model.HasPrecomputedData(this->GetManagerXMLTag());

  m_UseGPU = useGPU; 
  mp_Solver = &r_solver; 
  mp_Model = &r_model;

  mp_Updater = NULL;
  m_BVMargin = m_CloseDistance = std::numeric_limits<float>::quiet_NaN();  
  m_CoordinateHistorySize = 10;

  assert(this->GetModel().DoDeformableDeformableContactHandling() || this->GetModel().GetNumberOfRigidContactSurfaces() > 0);
  m_BVType = this->GetModel().GetBoundingVolumeType();

#ifndef GPU_GP_CONTACT
  if (UseGPU()) this->InstallUpdater(new _SetupNotSupportedUpdater());
  else {
#endif
    this->InstallUpdater(this->CreateRootUpdater());
#ifndef GPU_GP_CONTACT
  }
#endif    
  
  if (doLoadFromXML) {
    _XMLImporter importer(this->GetModel(), hasMembrane, this->UseGPU());

    importer.SetOuputObject(*this);
    importer.SetRootNode(this->GetModel().GetPrecomputedData(this->GetManagerXMLTag()));
    importer.Import();
    
    if ((int)this->GetDeformableSurface<tledDeformableContactSurface>().GetVolume2SurfaceNodeMap().size() != this->GetModel().GetMesh()->GetNumNodes()) {
      tledFatalError("Node number mismatch between contact surface and solid mesh. Does precomputation file match simulation file?");
    }
  } else {
    if (hasMembrane) {
      tledSurface *p_membrane = this->GetModel().GetGenericShellMesh();
    
      this->AddDeformableMembraneSurface(*p_membrane, *this->GetModel().GetMesh());
    
      delete p_membrane;
    } else {
      this->AddDeformableSurface(*this->GetModel().GetMesh());
    }  

    for (int si = 0; si < this->GetModel().GetNumberOfRigidContactSurfaces(); si++) {
      this->AddRigidSurface(this->GetModel().GetRigidContactSurfaceDefinition(si));
    }
  }

  if (hasMembrane) this->InstallUpdater(new _MembraneContactSurfaceUpdater(this->GetDeformableSurface<tledDeformableMembraneContactSurfaceCPU>(), r_solver)); 
  if (!std::isnan(this->GetModel().GetDeformableFrictionCoefficient())) this->GetDeformableSurface<tledContactSurface>().SetFrictionCoefficient(this->GetModel().GetDeformableFrictionCoefficient());
  if (!doLoadFromXML) this->GetDeformableSurface<tledDeformableContactSurface>().Init();

  {
    std::vector<float> mass(this->GetSolver().GetMesh()->GetNumNodes());

    this->GetSolver().GetMassVector(&mass.front());
    this->GetDeformableSurface<tledDeformableContactSurface>().InitNodeMasses(&mass.front());
  }

  if (!doLoadFromXML) {
    float defBVTau, modRateDist, modOff;

    modOff = r_model.GetContactSafetyMargin();
    modRateDist = r_model.GetRateConstraintDistance();

    if (modOff < 0 || modRateDist < 0) {
      tledFatalError("Rate constraint distance, contact offset distance have to be >= 0.");
    }

    this->GetDeformableSurface<tledDeformableContactSurface>().ComputeDiameters();    

    if (std::isnan(modOff) && std::isnan(modRateDist)) {
      defBVTau = 1e-2f*(this->GetDeformableSurface<tledDeformableContactSurface>().GetMinDiameter() + this->GetDeformableSurface<tledDeformableContactSurface>().GetMaxDiameter());
      modRateDist = defBVTau/2;
      modOff = 1e-2f*modRateDist;
    } else if (std::isnan(modOff)) {
      modOff = std::min(0.1f*modRateDist, 1e-4f*(this->GetDeformableSurface<tledDeformableContactSurface>().GetMinDiameter() + this->GetDeformableSurface<tledDeformableContactSurface>().GetMaxDiameter()));      
      defBVTau = std::max(modRateDist*2, 1e-3f*(this->GetDeformableSurface<tledDeformableContactSurface>().GetMinDiameter() + this->GetDeformableSurface<tledDeformableContactSurface>().GetMaxDiameter()));
    } else {
      assert(std::isnan(modRateDist));
      modRateDist = 10*modOff;
      defBVTau = std::max(modRateDist*2, 1e-3f*(this->GetDeformableSurface<tledDeformableContactSurface>().GetMinDiameter() + this->GetDeformableSurface<tledDeformableContactSurface>().GetMaxDiameter()));
    }

    this->SetBVMargin(defBVTau);    
    this->SetCloseDistance(modRateDist);
    this->SetSafetyMargin(modOff);
  }
  
  m_DoMultiPassContactHandling = this->GetModel().DoMultiPassContactHandling();

  if (!doLoadFromXML) _CreateBVHs(hasMembrane);

  assert(!std::isnan(GetBVMargin()));
  m_Dt = float(r_model.GetTimeStep());

  this->CreateSolvers(hasMembrane);

#ifdef GPU_GP_CONTACT
  if (this->UseGPU()) {
    _InitGPUSetup();
  }
#endif
}

bool tledUnstructuredContactManager::ComputeDeformableDeformableContactResponses(float *p_F, const float UNext[], const float UCurr[]) { 
  return static_cast<tledDeformableDeformableContactSolverCPU*>(mp_DeformableDeformableSolver)->ComputeContactResponses(p_F, UNext, UCurr); 
}

bool tledUnstructuredContactManager::ComputeDeformableRigidContactResponses(float *p_F, const float UNext[], const float UCurr[]) { 
  bool hadContacts = false;

  assert(!this->UseGPU());
  for (std::vector<tledDeformableRigidContactSolver*>::iterator i_solver = mvp_DeformableRigidContactSolvers.begin(); i_solver < mvp_DeformableRigidContactSolvers.end(); i_solver++) {
    hadContacts |= static_cast<tledDeformableRigidContactSolverCPU*>(*i_solver)->ComputeContactResponses(p_F, UNext, UCurr);
  }

  return hadContacts;
}

tledUnstructuredContactManager::~tledUnstructuredContactManager() {  
  if (mp_Updater != NULL) delete mp_Updater;

  if (mp_DeformableDeformableSolver != NULL) delete mp_DeformableDeformableSolver;
  for (std::vector<tledDeformableRigidContactSolver*>::iterator ip_s = mvp_DeformableRigidContactSolvers.begin(); ip_s < mvp_DeformableRigidContactSolvers.end(); ip_s++) {
    delete *ip_s;
  }

  if (mp_DeformableBVH != NULL) delete mp_DeformableBVH;
  for (std::vector<tledBVH*>::iterator ip_bvh = mvp_RigidSurfaceBVHs.begin(); ip_bvh < mvp_RigidSurfaceBVHs.end(); ip_bvh++) {
    if (*ip_bvh != NULL) delete *ip_bvh;
  }

  if (mp_DeformableBVH != NULL) delete mp_DeformableSurface;
  for (std::vector<tledRigidContactSurface*>::iterator ip_rigSurf = mvp_RigidSurfaces.begin(); ip_rigSurf < mvp_RigidSurfaces.end(); ip_rigSurf++) {
    delete *ip_rigSurf;
  }
}

void tledUnstructuredContactManager::SetRigidBVH(tledBVH *p_bvh, const int surfaceIndex) {
  assert(surfaceIndex < this->GetNumberOfRigidSurfaces());
  while (surfaceIndex >= (int)mvp_RigidSurfaceBVHs.size()) mvp_RigidSurfaceBVHs.push_back(NULL);
  mvp_RigidSurfaceBVHs[surfaceIndex] = p_bvh;
}

void tledUnstructuredContactManager::Update() {
  mp_Updater->Update(*this); 
}

void tledUnstructuredContactManager::InstallUpdater(Updater *p_updater) {
  p_updater->SetNextInChain(mp_Updater);
  mp_Updater = p_updater;
}

void tledUnstructuredContactManager::FinishContactHandling(void) {
#ifdef GPU_GP_CONTACT
  if (this->UseGPU()) {
    this->GetDeformableSurface<tledDeformableContactSurfaceGPU>().Update(static_cast<const tledSolverGPU&>(this->GetSolver()).GetAllOnDeviceNextDisplacements());
  } else {
#endif
    this->GetDeformableSurface<tledDeformableContactSurfaceCPU>().Update(this->GetSolver().GetAllNextDisps());
#ifdef GPU_GP_CONTACT
  }
#endif
  this->GetDeformableSurface<tledDeformableContactSurface>().Save();
}

void tledUnstructuredContactManager::AddDeformableSurface(const tledMesh &mesh) {
  mp_DeformableSurface = tledDeformableContactSurface::CreateSurface(mesh, this->UseGPU());
}

void tledUnstructuredContactManager::AddDeformableMembraneSurface(const tledSurface &membrane, const tledMesh &mesh) {
  mp_DeformableSurface = tledDeformableMembraneContactSurface::CreateSurface(mesh, membrane, this->UseGPU());
}

void tledUnstructuredContactManager::AddDeformableMembraneSurface(const XMLNode &rootNode) {
  mp_DeformableSurface = tledDeformableMembraneContactSurface::CreateSurface(rootNode, this->UseGPU());
  mp_DeformableSurface->SetCoordinateHistorySize(this->GetCoordinateHistorySize());
}

void tledUnstructuredContactManager::AddDeformableSurface(const XMLNode &rootNode) {
  mp_DeformableSurface = tledDeformableContactSurface::CreateSurface(rootNode, this->UseGPU());
  mp_DeformableSurface->SetCoordinateHistorySize(this->GetCoordinateHistorySize());
}

void tledUnstructuredContactManager::ResetRigidSurfaces() {
  for (int s = 0; s < this->GetNumberOfRigidSurfaces(); s++) this->GetRigidSurface<tledRigidContactSurface>(s).ResetNodes();
}

void tledUnstructuredContactManager::AddRigidSurface(const XMLNode &meshSpecRoot) {
  std::vector<int> slaveNodeInds;

  mvp_RigidSurfaces.push_back(tledRigidContactSurface::CreateSurface(meshSpecRoot, this->UseGPU()));

  if (meshSpecRoot.nChildNode("Motion") > 0) {
    if (!this->UseGPU()) {
      static_cast<tledMovingRigidContactSurfaceCPU*>(mvp_RigidSurfaces.back())->SetHistoryLength(this->GetCoordinateHistorySize());
    } else {
#ifdef GPU_GP_CONTACT
      static_cast<tledMovingRigidContactSurfaceGPU*>(mvp_RigidSurfaces.back())->SetHistoryLength(this->GetCoordinateHistorySize());
#else
      tledFatalFeatureNotEnabledError;
#endif
    }
  }

  if (meshSpecRoot.nChildNode("SlaveNodes") > 0) {
    XMLNode slv = meshSpecRoot.getChildNode("SlaveNodes");
    std::istringstream iss(slv.getAttribute("NumNodes"));
    int numSlaveNodes;

    iss >> numSlaveNodes;  
    if (!iss.fail() && numSlaveNodes > 0) {
      slaveNodeInds = GetXMLTextAsVector<int>(slv);
      if (slaveNodeInds.size() == 1 && numSlaveNodes > 1) slaveNodeInds = tledSequenceGenerator::MakeSequence(slaveNodeInds[0], slaveNodeInds[0] + numSlaveNodes);
    } else {
      tledFatalError("Need a valid NumNodes attribute in <SlaveNodes> tag.");
    }
  } else {
    if (mp_DeformableSurface == NULL) {
      tledFatalError("Need a deformable surface before first rigid contact surface can be instantiated.");
    }
    slaveNodeInds = tledSequenceGenerator::MakeSequence(0, mp_DeformableSurface->GetNumberOfNodes());
  }

  mvp_RigidSurfaces.back()->SetSlaveNodeIndices(slaveNodeInds, mp_DeformableSurface->GetNumberOfNodes());
}
