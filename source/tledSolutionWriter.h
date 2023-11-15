// =========================================================================
// File:       tledSolutionWriter.h
// Purpose:    Compile model solutions for saving
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    Dec 2009
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifndef tledSolutionWriter_H
#define tledSolutionWriter_H

#include "tledSolver.h"
#include "tledModel.h"

#include <vector>
#include <string>

/**
 * \brief Periodically stores simulation variables and has the ability to write them to a file.
 * \ingroup export
 */
class tledSolutionWriter
{
  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledSolutionWriter(void);
  tledSolutionWriter(tledSolver* Solver, tledModel* Model, int numSteps);
  ~tledSolutionWriter(void);
  /** @} */

  /**
   * \name Saving and F/S Output
   * @{
   */
public:
  void SaveSolution(int step);
  
  /** Sets the file prefix where the output files are stored. Defaults to tledModel::GetDirectory */
  void SetFilePrefix(const std::string &pfx);
  
  /** File prefix where the output files are stored. */
  const std::string& GetFilePrefix(void) const { return m_FilePrefix; }

  /** Writes displacements and internal forces to a file U.txt and F.txt, respectively, in the directory holding the simulation XML. */
  void WriteSolution(void);
  void InitialiseHistories(void);
  /** @} */

  /**
   * \name Save Setting Overrides and Queries
   * @{
   */
public:
  void SetDoSaveU(const bool doSave) { saveDisplacements = doSave; }
  void SetDoSaveF(const bool doSave) { saveForces = doSave; }
  void SetDoSaveEKinetic(const bool doSave) { saveEKin = doSave; }
  void SetDoSaveEStrain(const bool doSave) { saveEStrain = doSave; }

  bool DoSaveU(void) const { return saveDisplacements; }
  bool DoSaveF(void) const { return saveForces; }
  bool DoSaveEKinetic(void) const { return saveEKin; }
  bool DoSaveEStrain(void) const { return saveEStrain; }

  /** Total step count setter. NB: this also resets the solution histories */
  void SetNumSteps(const int numSteps); 

  /** Save frequency setter (in N. of steps). Does not reallocate memory, requires a call to SetNumSteps after change. */
  void SetFrequency(const int frq) { m_Freq = frq; }

  int GetFrequency(void) const { return m_Freq; }
  /** @} */

private:
  tledSolver* solver;
  tledModel* model;
  bool saveDisplacements;
  bool saveForces;
  bool saveStresses;
  bool saveStrains;
  bool saveEKin;
  bool saveEStrain;
  std::vector<float>* U; // Transient displacement solutions
  std::vector<float>* F; // Transient internal force solutions
  std::vector<float>* S; // Transient 2nd Piola-Kichhoff stress solutions
  std::vector<float>* E; // Transient Green's strain solutions
  std::vector<float> m_EKin; // Total kinetic energy
  std::vector<float> m_EStrain; // Total strain energy
  std::string m_FNameU;   // File name for U
  std::string m_FNameF;   // File name for F
  std::string m_FNameS;   // File name for S
  std::string m_FNameE;   // File name for E
  int m_Freq;  // Save the solution every freq steps
  int nSvs;   // Number of saved time steps
  int Scnt; // Counter for solution saves
  std::string m_FilePrefix;
};


#endif // tledSolutionWriter_H
