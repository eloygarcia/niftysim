// =========================================================================
// File:       tledContactManager.h
// Purpose:    Manage all contact objects in the model
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    July 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifndef tledContactManager_H
#define tledContactManager_H

#include "tledContactCylinder.h"
#include "tledContactUSProbe.h"
#include "tledContactPlate.h"
#include "tledModel.h"
#include "tledUnstructuredContactManager.h"

#ifdef _GPU_
struct tledGPUContacts
{
   tledGPUContactCylinder** Cyls;
   int NumContactCyls;
   tledGPUContactUSProbe** Prbs;
   int NumContactPrbs;
   tledGPUContactPlate** Plts;
   int NumContactPlts;
};
#endif // _GPU_

/** 	 
 * \defgroup contact Contact Modelling Module
 */ 	 
	  	 
/** 	 
 * \brief Manager for contact objects
 * \ingroup contact 	 
 */
class tledContactManager
{
  /**
   * \name Construction, Destruction
   * @{
   */
public:
   tledContactManager(void) : Model(NULL) {;}
   tledContactManager(tledModel* model);
   virtual ~tledContactManager(void);
   /** @} */

  /**
   * \name Analytic Contact Surfaces
   * @{
   */
public:
   /** Updates configuration of analytical contact surfaces to current time step (time = TR) */
   void Update(double TR);

   tledContactCylinder* GetContactCyl(int cylNum) {return &(Cyls[cylNum]);}
   tledContactUSProbe* GetContactPrb(int prbNum) {return &(Prbs[prbNum]);}
   tledContactPlate* GetContactPlt(int pltNum) {return &(Plts[pltNum]);}

   int GetNumContactCyls() {return numContactCyls;}
   int GetNumContactPrbs() {return numContactPrbs;}
   int GetNumContactPlts() {return numContactPlts;}
   int GetNumContactObjects() {return numContactCyls+numContactPrbs+numContactPlts;}

   // Set contact cylinder params
   void SetContactCylDisp(int num, std::vector<float> disp);
   void SetContactCylRadChange(int num, float dr);
   void SetContactCylOrigin(int num, std::vector<float> orig);
   void SetContactCylAxis(int num, std::vector<float> axis);
   void SetContactCylRadius(int num, float r);
   void SetContactCylLength(int num, float l);

   // Set contact US probe params
   void SetContactPrbDisp(int num, std::vector<float> disp);
   void SetContactPrbRadChange(int num, float dr);
   void SetContactPrbOrigin(int num, std::vector<float> orig);
   void SetContactPrbAxis(int num, std::vector<float> axis);
   void SetContactPrbRadius(int num, float r);
   void SetContactPrbLength(int num, float l);

   // Set contact plate params
   void SetContactPltDisp(int num, std::vector<float> disp);
   void SetContactPltStartCrnrA(int num, std::vector<float> A);
   void SetContactPltStartCrnrB(int num, std::vector<float> B);
   void SetContactPltStartCrnrC(int num, std::vector<float> C);
   
#ifdef _GPU_
   tledGPUContacts* GetContactsDevicePointer(void);
#endif // _GPU_
  /** @} */

  /**
   * \name General Queries
   * @{
   */
public:
  tledModel& GetModel(void) { return *Model; }
  /** @} */

  /**
   * \name Unstructured-Mesh Contact (Only available in CPU builds)
   * @{
   */
private:
  tledUnstructuredContactManager *mp_UnstructuredContacts;

protected:
  /** Overridable factory function for unstructured contact managers, called by InitialiseUnstructuredContacts, no init should be performed by overriding functions. */
  virtual tledUnstructuredContactManager* CreateUnstructuredContactManager(void);

public:
  tledUnstructuredContactManager& GetUnstructuredContactManager(void) { return *mp_UnstructuredContacts; }
  const tledUnstructuredContactManager& GetUnstructuredContactManager(void) const { return *mp_UnstructuredContacts; }

  /** \brief Initialises unstructured contacts if such are requested by the model */
  void InitialiseUnstructuredContacts(tledSolver &r_solver, const bool useGPU);

  /** \return true if unstructured mesh contacts are to be taken into account. */
  bool DoUnstructuredContacts(void) const { return mp_UnstructuredContacts != NULL; }
  /** @} */

  /**
   * \name XML Export
   * @{
   */
public:
  XMLNode ExportToXML(void) const { return mp_UnstructuredContacts != NULL? this->GetUnstructuredContactManager().ExportToXML() : XMLNode(); }
  /** @} */
  
private:
   int numContactCyls;
   int numContactPrbs;
   int numContactPlts;
   std::vector<tledContactCylinder> Cyls;
   std::vector<tledContactUSProbe> Prbs;
   std::vector<tledContactPlate> Plts;
   tledModel* Model;
   
#ifdef _GPU_
   tledGPUContacts h_GPUContacts;
   tledGPUContacts* d_GPUContacts;
#endif // _GPU_
};

#endif // tledContactManager_H

