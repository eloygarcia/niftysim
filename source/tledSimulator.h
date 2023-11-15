// =========================================================================
// File:       tledSimulator.h
// Purpose:    Simulation object that includes a solver and supporting
//             structures, and executes a complete simulation
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    December 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#ifndef tledSimulator_H
#define tledSimulator_H

#include "tledHelper.h"
#include "tledSolver.h"
#include "tledSolverCPU.h"
#include "tledConstraintManager.h"
#include "tledConstraint.h"
#include "tledModel.h"
#include "tledTimer.h"
#include "tledSolutionWriter.h"
#include "tledContactManager.h"

#ifdef _GPU_
#include "tledSolverGPU.h"
#include "tledSolverGPU_ROM.h"
#endif   // _GPU_

#ifndef _NO_USING_NAMESPACE_STD_
/* Legacy library behaviour of using namespace std is now deprecated and can be turned off in CMake. */
using namespace std; 
#endif

/**
 * \defgroup simulator Simulator
 * \brief Main interface to TLED implementation
 */

/** 
 * \brief Main simulation class performing the actual time-stepping loop.
 * \ingroup simulator
 */
class tledSimulator
{
  /**
   * \name Construction, Destruction
   * @{
   */
private:
  void _BasicInit(tledModel *p_model);

public:
   tledSimulator(void);
   virtual ~tledSimulator(void);

  /** 
   * The simulation is initialised with data in the model file.
   * However, many of the simulation parameters may be modified prior
   * to running a simulation. This allows multiple simulations to
   * be run with a given model.
   */
   tledSimulator(tledModel* model, const bool useGPU = false, const bool verbose = false, const bool doTiming = false);

   /**
    * Only performs the basic setup, requires a separate call to Init prior to execution of the actual simulation.
    */
   tledSimulator(tledModel &r_model);

   /**
    * Performs the full object initialisation
    */
   virtual void Init(const bool requestGPU, const bool verbose, const bool timing);
   /** @} */
   
  /**
   * \name Processing
   * @{
   */
protected:
  /** Initialisation of solution variables, etc. */
  virtual void PerformMainLoopInit(void);

  /** Simulation main loop. Returns 0 on success, non-zero otherwise. */
  virtual int MainLoop(void);
  
public:
   /** Main simulation loop. Returns 0 on success, non-zero otherwise. */
   virtual int Simulate(void);

  /** Deletes all constraints. */
   void ResetConstraints(void);
  /** @} */

  /**
   * \name General Setup Queries
   * @{
   */
public:
  bool UseGPU(void) const { return m_UseGPU; }
  int GetError(void) const {return error;}
  double GetSimulationTime(void) const { return Timer->GetDuration(); } // Elapsed time in ms  

  /** Outputs additional diagnostic messages. */
  bool IsVerbose(void) const { return Verbose; }

  /** Performs time benchmarking */
  bool DoesTiming(void) const { return Timing; }
  /** @} */

   /**
    * \name Functions for setting (changing) simulation params...
    * @{
    */
public:
   void SetSimulationTimeStep(double dt);
   const double& GetSimulationTimeStep(void) const { return Dt; }

   void SetSimulationTotalTime(double t);
   const double& GetSimulationTotalTime(void) const { return T; }

   int GetNumberOfSimulationTimeSteps(void) const { return NumSteps; }
   /** @} */

  /**
   * \name Simulator Components
   * @{
   */
protected:
  /** Overridable solver factory function. Initialisation performed by Init. At time of call: only model is available. */
  virtual tledSolver* CreateSolver(void);

  /** Overridable contact manager factory function. */
  virtual tledContactManager* CreateContactManager(void);

  /** Overridable constraint manager factory function. At time of call: only model, solver are available, no contacts. Initialisation performed outside factory function. */
  virtual tledConstraintManager* CreateConstraintManager(void);

  tledTimer& GetTimer(void) { return *Timer; }
  const tledTimer& GetTimer(void) const { return *Timer; }
  
public:
   tledSolver* GetSolver(void) {return Solver;}
   const tledSolver* GetSolver(void) const {return Solver;}

   void SetModel(tledModel &r_model) { Model = &r_model; }
   tledModel* GetModel(void) {return Model;}
   const tledModel* GetModel(void) const {return Model;}

   tledConstraintManager* GetConstraintManager(void) {return Constraints;}
   const tledConstraintManager* GetConstraintManager(void) const {return Constraints;}

   tledContactManager* GetContactManager(void) {return Contacts;}
   const tledContactManager* GetContactManager(void) const {return Contacts;}
   /** @} */
   
  /**
   * \name I/O
   * @{
   */
protected:
  tledSolutionWriter& GetSolutionWriter(void) { return *Writer; }
  const tledSolutionWriter& GetSolutionWriter(void) const { return *Writer; }

public:
  /** Sets the file prefix where the output files are stored. Defaults to tledModel::GetDirectory. Requires object initialisation prior to call. */
  void SetIOFilePrefix(const std::string &pfx) { this->GetSolutionWriter().SetFilePrefix(pfx); }
  
  /** File prefix where the output files are written. */
  const std::string& GetIOFilePrefix(void) const { return this->GetSolutionWriter().GetFilePrefix(); }
  /** @} */

private:
   tledModel* Model;
   tledSolver* Solver;
   tledConstraintManager* Constraints;
   tledContactManager* Contacts;
   tledSolutionWriter* Writer;
   tledTimer* Timer;
   
   double Dt;
   double T;
   int NumSteps;
   
   bool Verbose;
   bool Timing;
   bool m_IsInitialised;
   bool m_UseGPU;
   int error;
};

#endif // tledSimulator
