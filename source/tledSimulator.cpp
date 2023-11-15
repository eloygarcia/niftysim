// =========================================================================
// File:       tledSimulator.cpp
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

#include "tledSimulator.h"
#include "tledMatrixFunctions.h"
#include "tledSubModelManager.h"

#ifdef _USE_BOOST_
#include "tledParallelSolverCPU.h"
#endif

using namespace std;

tledSimulator::~tledSimulator()
{
   if (Solver)
      delete Solver;
   if (Constraints)
      delete Constraints;
   if (Writer)
      delete Writer;
   if (Timer)
      delete Timer;
   if (Contacts)
      delete Contacts;
}

tledSolver* tledSimulator::CreateSolver() {
  tledSolver *p_solver = NULL;

#ifdef _GPU_
  if (this->UseGPU())
   {
      if (this->GetModel()->GetROM() == true)
	p_solver = new tledSolverGPU_ROM;
      else
	p_solver = new tledSolverGPU;
   }
   else
#endif // _GPU_
   {
#ifdef _USE_BOOST_
     if (this->GetModel()->GetMesh()->GetNumEls() > 1024) p_solver = new tledParallelSolverCPU;
     else p_solver = new tledSolverCPU;
#else
     p_solver = new tledSolverCPU;
#endif
   }

   return p_solver;
}

tledContactManager* tledSimulator::CreateContactManager() {
  return new tledContactManager(this->GetModel());
}

tledConstraintManager* tledSimulator::CreateConstraintManager() {
  return new tledConstraintManager();
}

void tledSimulator::_BasicInit(tledModel *p_model) {
  error = 0;
  Verbose = false;
  Timing = false;
  Model = p_model;

  // Compute NumSteps =========================================
  Dt = this->GetModel()->GetTimeStep();
  T = this->GetModel()->GetTotalTime();
  NumSteps = (int)(T/Dt + 0.5);
}

void tledSimulator::Init(const bool requestGPU, const bool verbose, const bool doTiming) {
  m_UseGPU = requestGPU;
  Verbose = verbose;
  Timing = doTiming;

#ifndef _GPU_
  if (this->UseGPU()) {
    tledLogErrorStream(tledHelper::Warning() << "GPU solver was requested but it was not enabled at compile time. Proceeding in CPU mode.");
    m_UseGPU = false;
  }
#endif
   
  // Initialise Solver ===========================================
  Solver = this->CreateSolver();
      
  if (Verbose)
    cout << "Beginning precomp...";
  Solver->Init(Model);
  if (Verbose)
    cout << "Done!" << endl;
      
  // Define model constraints ======================
  if (Verbose)
    cout << "Assembling constraints...";

  Constraints = this->CreateConstraintManager();
  Constraints->Init(*Solver, *Model);

  if (Verbose)
    cout << Constraints->GetNumConstraints() << " constraints constructed" << endl;
         
  // Define contacts ==============================
  if (Verbose)
    cout << "Assembling contacts...";
  if ((!this->UseGPU())&(this->GetModel()->GetNumContactObjects()>0))
    {
      tledLogErrorStream(tledHelper::Warning() << "Sorry, contact modelling only available in GPU mode.");
      error = 1;
    }
  Contacts = this->CreateContactManager();
  Solver->SetContactManager(Contacts);
  Contacts->InitialiseUnstructuredContacts(*this->GetSolver(), this->UseGPU());
  if (Verbose)
    cout << Contacts->GetNumContactObjects() << " contacts constructed" << endl;
      
  // Setup output requests =========================================
  Writer = new tledSolutionWriter(Solver,Model,NumSteps);
   
  // Instantiate Timer =========================================
  Timer = new tledTimer();  

  m_IsInitialised = true;
}

tledSimulator::tledSimulator(tledModel* model, const bool requestGPU, const bool verbose, const bool doTiming) {
  _BasicInit(model);
  this->Init(requestGPU, verbose, doTiming);
}

tledSimulator::tledSimulator(tledModel &r_model) {
  Writer = NULL;
  Timer = NULL;
  Solver = NULL;
  Contacts = NULL;
  Constraints = NULL;
  m_IsInitialised = false;

  _BasicInit(&r_model);  
}

tledSimulator::tledSimulator() {
  Model = NULL;
  Writer = NULL;
  Timer = NULL;
  Solver = NULL;
  Contacts = NULL;
  Constraints = NULL;
  m_IsInitialised = false;
}

void tledSimulator::PerformMainLoopInit() {
  Solver->InitialiseSolutionVariables();
  Writer->InitialiseHistories();
   
  Solver->SetFixed(Constraints->GetFixInd(0),Constraints->GetFixInd(1),Constraints->GetFixInd(2));
  Solver->SetDisps(Constraints->GetDispInd(0),Constraints->GetDispInd(1),Constraints->GetDispInd(2));
  Solver->SetExtForces(Constraints->GetForceInd(0),Constraints->GetForceInd(1),Constraints->GetForceInd(2));
}

int tledSimulator::MainLoop() {
   bool Divergence = false;
   int failStep = 0;
   int progress = 0;

   if (Timing)
      Timer->StartTimer();
   
   for (int step = 0; step < NumSteps+1; step++)
   {
      if (Verbose)
      {
         if ((float)(step+1)/(float)NumSteps*100 >= progress + 10)
         {
            progress += 10;
            cout << "\tstep: " << step+1 << "\t(" << progress << "%)" << endl;
         }
      }
      // Check for divergence
      if ((step+1)%(int)(NumSteps*0.01 + 1) == 0)	// Check every 1% of simulation time
      {
         Solver->GetDivergence(&Divergence);
         if (Divergence == true)
         {
            failStep = step+1;
            break;
         }
      }
      Solver->SetDisps(Constraints->GetDispInd(0),Constraints->GetDispVal(0,step,Dt,T),
                       Constraints->GetDispInd(1),Constraints->GetDispVal(1,step,Dt,T),
                       Constraints->GetDispInd(2),Constraints->GetDispVal(2,step,Dt,T));
      Solver->SetExtForces(Constraints->GetForceInd(0),Constraints->GetForceVal(0,step,Dt,T),
                           Constraints->GetForceInd(1),Constraints->GetForceVal(1,step,Dt,T),
                           Constraints->GetForceInd(2),Constraints->GetForceVal(2,step,Dt,T));
#ifdef _GPU_
      Contacts->Update(Dt*(step+1)/T);
#endif
      Solver->PerformStep();

      // Save results
      Writer->SaveSolution(step);
   }
   
   if (Timing)
      Timer->StopTimer();
   
   if (Divergence) {
      cout << "Failed!" << endl;
      cout << "\n!!!Divergence after " << failStep << " steps -> reduce time step.\n" << endl;
      cout << "Simulation terminated" << endl;
      return 1;
   } else {
      if (Verbose) cout << "Done!" << endl;

      return 0;
   }
}

int tledSimulator::Simulate() {
  int mainLoopExitCode;

  if (Model == NULL) {
    tledFatalError("Have no model, cannot proceed.");
  }

  if (!m_IsInitialised) {
    tledLogErrorStream(tledHelper::Warning() << "Simulate called w/o prior call to Init. Initialising with default values: GPU = off, verbose = off, timing = off.");
    this->Init(false, false, false);
  }

  this->PerformMainLoopInit();

   if (Verbose)
   {
      cout << "\nSIMULATION DATA" << endl;
      cout << "Time Step:\t" << Dt << endl;
      cout << "Time Total:\t" << T << endl;
      cout << "# Steps:\t" << NumSteps << endl;
      cout << "\nSimulating..." << endl;
   }
  
   mainLoopExitCode = this->MainLoop();

   // Write solution results to file =============================
   Writer->WriteSolution();
   
   return mainLoopExitCode;
}

void tledSimulator::SetSimulationTimeStep(double dt)
{
   Dt = dt;
   NumSteps = (int)(T/Dt + 0.5);
   Solver->SetTimeStep(Dt);
   Writer->SetNumSteps(NumSteps);
}

void tledSimulator::SetSimulationTotalTime(double t)
{
   T = t;
   NumSteps = (int)(T/Dt + 0.5);
}

void tledSimulator::ResetConstraints(void)
{
   Solver->InitialiseConstraints();
   if (Constraints != NULL) delete Constraints;
   Constraints = new tledConstraintManager();
}



