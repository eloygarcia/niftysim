// =========================================================================
// File:       niftysim.cpp
// Purpose:    Main function definition for NiftySim
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

#include "tledSolver.h"
#include "tledModel.h"
#include "tledTimer.h"
#include "tledShellSolver.h"
#include "tledSimulator.h"
#include "tledSubModelManager.h"
#include "tledCUDAHelpers.h"

#ifdef _Visualisation_
#include "tledModelViewer.h"
#include "tledVTKMeshExporter.h"
#include "tledVTKSurfaceExporter.h"
#endif // _Visualisation_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <cstdio>
#include <cstring>

struct InputParameters
{
   InputParameters()
   {
      xmlFName = NULL;	// Required param -> no default
      printNForces = false;
      printNDispForces = false;
      printNDispForcesSums = false;
      printNDisps = false;
      printStrainEnergy = false;
      printKineticEnergy = false;
      plotModel = false;
      sportMode = false;
      doTiming = false;
      verbose = false;
#ifdef _Visualisation_
      meshExportPath = NULL;
      shellExportPath = NULL;
      subMeshExportPrefix = NULL;
      doExportNodeCompression = false;
#endif
      precompImportFName = NULL;
      precompExportFName = NULL;
      outputPrefix = NULL;
#ifdef _GPU_
      cudaDeviceId = -1;
#endif
   }
   char* xmlFName;
   bool printNForces;	// Print all nodal forces
   bool printNDispForces;	// Print loaded nodal forces
   bool printNDispForcesSums;	// Print sums of loaded nodal forces
   bool printNDisps;	// Print all nodal displacements
   bool printStrainEnergy;
   bool printKineticEnergy;
   bool plotModel; // Plot the results
   bool sportMode; // Use GPU solver
   bool doTiming; // Report execution time   
   bool verbose;
#ifdef _GPU_
  int cudaDeviceId; // CUDA device to be used for computation. Only evaluated in sport mode.
#endif
#ifdef _Visualisation_  
   bool doExportNodeCompression; // Remove unreferenced nodes from exported VTK meshes.
   char const *meshExportPath; // VTK mesh/displacement export
   char const *shellExportPath; // VTK shell-mesh/displacement export
   char const *subMeshExportPrefix;  // Prefix for VTK export of submeshes
#endif
   char const *precompImportFName; // Path to XML file holding precomputed geometry data used for initialisation
   char const *precompExportFName; // Destination path for precomputation XML
   char const *outputPrefix;  // Prefix for U.txt, F.txt, etc.
};

void Usage(char *name);
int parseInputs(int argc, char **argv, InputParameters* inParams);
void PrintAnalysisSettings(tledModel* Model, InputParameters* inParams);
void PrintIntro();
void PlotMeshes(tledSimulator* Simulator, InputParameters* inParams);
void PrintResults(tledSolver* Solver, InputParameters* inParams);


int main(int argc, char **argv)
{
   // Parse inputs ================================================
   InputParameters* inParams = new InputParameters;
   int inputResult = parseInputs(argc,argv,inParams);
   if (inputResult == 1)
   {
      std::cout << "Exiting" << std::endl;
      return EXIT_FAILURE;
   }
   else if (inputResult == 2)
   {
      std::cout << "Exiting" << std::endl;
      return EXIT_SUCCESS;
   }
   
   if (inParams->verbose)
      PrintIntro();

   // Load Model ===========================================
   tledModel* Model = new tledModel(inParams->xmlFName);
   if (Model->GetError() > 0)
   {
      std::cout << "Exiting" << std::endl;
      return EXIT_FAILURE;
   }

   if (inParams->verbose)
      PrintAnalysisSettings(Model,inParams);
   
#ifdef _GPU_
   if (inParams->cudaDeviceId < 0) {
#ifndef _CUDA_5PLUS_SDK
     inParams->cudaDeviceId = cutGetMaxGflopsDeviceId();
#else
     inParams->cudaDeviceId = gpuGetMaxGflopsDeviceId();
#endif
   }
   cudaSetDevice(inParams->cudaDeviceId);   
#endif

   // Import of precomputed data ===========================================
   if (inParams->precompImportFName != NULL) {
     Model->SetPrecomputedGeometryFilename(inParams->precompImportFName);
   }

   // Construct Simulator ===========================================
   tledSimulator* Simulator = new tledSimulator(Model,inParams->sportMode,inParams->verbose,inParams->doTiming);
   int error = Simulator->GetError();
   if (error > 0)
   {
      std::cerr << "!!! Problems when initialising the simulation\n--> exiting." << std::endl;
      return EXIT_FAILURE;
   }

   // Check if only precomputation was requested, if so just export XML and exit ===========================================
   if (inParams->precompExportFName != NULL) {
     std::cout << "Exporting precomputable data for repeated simulations to file " << inParams->precompExportFName << std::endl;
     Model->StartXMLExport();
     Model->GetExportXMLRootNode().addChild(Simulator->GetContactManager()->ExportToXML());
     if (Model->WritePrecomputableDataXML(inParams->precompExportFName)) return EXIT_SUCCESS;
     else {
       std::cerr << "Error exporting data:\n";
       std::perror(inParams->precompExportFName);

       return EXIT_FAILURE;
     }
   }

   if (inParams->outputPrefix != NULL) Simulator->SetIOFilePrefix(inParams->outputPrefix);
   
   // Simulate! ===========================================
   error = Simulator->Simulate();
   if (error > 0)
   {
      std::cerr << "!!! Problems when running the simulation\n--> exiting." << std::endl;
      return EXIT_FAILURE;
   }
   
   Simulator->GetSolver()->PrepareOutput();

   if (inParams->doTiming)
      std::cout << "Execution time: " << Simulator->GetSimulationTime() << " ms" << std::endl;

#ifdef _Visualisation_
   if (inParams->meshExportPath && Model->GetMesh()->GetNumEls() > 0) {
     tledVTKMeshExporter exporter;

     exporter.SetFileName(inParams->meshExportPath);
     if (inParams->doExportNodeCompression) exporter.SetDoNodeCompression(true);
     exporter.SetMesh(*Model->GetMesh());
     exporter.AddNodeVectorAttribute("displacements", Simulator->GetSolver()->GetAllDisps());
     
     if (!exporter.Write()) {
       std::cerr << "Export of mesh/final displacements to " << inParams->meshExportPath << " failed:\n";
       std::perror(inParams->meshExportPath);
     }
   } 

   if (inParams->shellExportPath && &Simulator->GetSolver()->GetShellSolver() != NULL && &Simulator->GetSolver()->GetShellSolver().GetSurface() != NULL && Simulator->GetSolver()->GetShellSolver().GetSurface().GetNumberOfFacets() > 0) {
     tledVTKSurfaceExporter exporter(inParams->shellExportPath);

     if (inParams->doExportNodeCompression) exporter.SetDoNodeCompression(true);
     exporter.SetMesh(Simulator->GetSolver()->GetShellSolver().GetSurface());
     exporter.AddNodeVectorAttribute("displacements", Simulator->GetSolver()->GetAllDisps());
     if (!exporter.Write()) {
       std::cerr << "Export of shell mesh/final displacements to " << inParams->shellExportPath << " failed:\n";
       std::perror(inParams->shellExportPath);
     }
   }

   if (inParams->subMeshExportPrefix && Model->HasSubModels()) {
     for (int sm = 0; sm < Model->GetSubModelManager().GetNumberOfSubMeshes(); sm++) if (Model->GetSubModelManager().GetNumberOfSubMeshElements(sm) > 0) {
	 tledMesh *p_subMesh = Model->GetSubModelManager().ExportSubMesh(sm, *Model->GetMesh());
	 std::vector<float> subMeshDisps(3*Model->GetSubModelManager().GetNumberOfSubMeshNodes(sm));
	 tledVTKMeshExporter exporter;
	 std::ostringstream oss;

	 Model->GetSubModelManager().ExportSubMeshNodeAttributes(&subMeshDisps.front(), Simulator->GetSolver()->GetAllDisps(), sm);

	 oss << inParams->subMeshExportPrefix << sm << ".vtk";
	 exporter.SetFileName(oss.str());
	 exporter.SetMesh(*p_subMesh);
	 exporter.AddNodeVectorAttribute("displacements", &subMeshDisps.front());

	 if (!exporter.Write()) {
	   std::cerr << "Export of sub-mesh/final displacements to " << oss.str() << " failed:\n";
	   std::perror(oss.str().c_str());
	 }

	 delete p_subMesh;
       }
   }
#endif

   // Print out solution results =============================
   PrintResults(Simulator->GetSolver(),inParams);
   
   // Plot the mesh ===============================================
   PlotMeshes(Simulator,inParams);

   if (inParams->verbose)
      std::cout << "\nProgram ended successfully" << std::endl;

   // Clean up ====================================================
   delete inParams;
   delete Model;
   delete Simulator;
#if defined _GPU_ && !defined _CUDA_3MINUS_SDK
   if (inParams->sportMode) {
     tledCheckCUDAErrors(cudaDeviceReset());
   }
#endif
      
   return EXIT_SUCCESS;

} // END MAIN() ****************************************************

void Usage(char *name)
{
   std::cout << "=====================================================================" << std::endl;
   std::cout << "Usage :" << std::endl;
   std::cout << name << " -x <file> [options]" << std::endl;
   std::cout << "*** [INPUTS] ***" << std::endl;
   std::cout << "\t-x\tXML file name <file>" << std::endl;
   std::cout << "*** [OPTIONS] ***" << std::endl;
   std::cout << "\t-output-prefix <PATH PFX>: Write outputs to given prefix (U.txt, F.txt, etc.).\n";
   std::cout << "\t-precompute <XML DESTINATION PATH>: Precomputes data for fast repeated simulations involving contacts with different material settings.\n";
   std::cout << "\t-initialise-with <XML PATH>: Loads precomputed data for quick initialisation from an XML file created with the -precompute switch.\n";
#ifdef _GPU_
   std::cout << "\t-sport\t\t\tUse GPU execution" << std::endl;
   std::cout << "\t-device <DEVICE ID>\tUse CUDA device with given ID for computation (use \"nvidia-smi -L\" to determine). By default uses fastest installed GPU." << std::endl;
#endif
   std::cout << "\t-print\t\t\tPrint solution results <NF,LNF,LNFS,ND,ES,EK>" << std::endl;
   std::cout << "\t\tNF:\tNodal forces" << std::endl;
   std::cout << "\t\tLNF:\tLoaded nodal forces" << std::endl;
   std::cout << "\t\tLNFS:\tLoaded nodal force sums" << std::endl;
   std::cout << "\t\tND:\tNodal displacements" << std::endl;
   std::cout << "\t\tES:\tStrain energy" << std::endl;
   std::cout << "\t\tEK:\tKinetic energy" << std::endl;
#ifdef _Visualisation_
   /* Disabled due to VTK 6/CUDA incompat., restore when fixed upstream. */
   std::cout << "\t-plot\t\t\t\tDisplay the model results" << std::endl;
   std::cout << "\t-export-mesh <PATH>:\t\tWrites the mesh and final displacements to a VTK unstructured grid file.\n";
   std::cout << "\t-export-membrane <PATH>:\tWrites the membrane/shell mesh and final displacements to a VTK unstructured grid file.\n";
   std::cout << "\t-export-submeshes <PREFIX>:\tExports the sub-meshes (need to be defined in simulation XML) with their individual final displacements.\n"
	<< "\t\t\t\t\tExport paths are given by <PREFIX>0.vtk, .., <PREFIX>N.vtk where N is the number of sub-meshes defined in the input XML.\n";
   std::cout << "\t-compress-exported-nodes:\tRemoves any unreferenced vertices (and corresponding attributes) from the exported meshes.\n";
#endif
   std::cout << "\t-t\tReport execution time (ms)" << std::endl;
   std::cout << "\t-v\tVerbose output" << std::endl;
   std::cout << "\t-help\tDisplay this message" << std::endl;
   std::cout << "=====================================================================" << std::endl;
   return;	
}

int parseInputs(int argc, char **argv, InputParameters* inParams)
{
   // Returns:
   // 0 to indicate normal program continuation
   // 1 to indicate erroneous program termination
   // 2 to indicate normal program termination
   for (int i = 1; i < argc; i++)
   {
      if (std::strcmp(argv[i], "-help")==0 || 
         std::strcmp(argv[i], "-Help")==0 || 
         std::strcmp(argv[i], "-HELP")==0 || 
         std::strcmp(argv[i], "-h")==0 || 
         std::strcmp(argv[i], "--h")==0)
      {
         Usage(argv[0]);
         return 2;
      }
      else if (string(argv[i]) == "-output-prefix") inParams->outputPrefix = argv[++i];
      else if (std::strcmp(argv[i], "-x") == 0) {
         inParams->xmlFName = argv[++i];
      }
#ifdef _Visualisation_
      else if (string(argv[i]) == "-export-mesh") inParams->meshExportPath = argv[++i];
      else if (string(argv[i]) == "-export-membrane") inParams->shellExportPath = argv[++i];
      else if (string(argv[i]) == "-export-submeshes") inParams->subMeshExportPrefix = argv[++i];
      else if (string(argv[i]) == "-compress-exported-nodes") inParams->doExportNodeCompression = true;
#endif
      else if (string(argv[i]) == "-device") {
#ifdef _GPU_
	if (i + 1 < argc) {
	  std::istringstream iss(argv[++i]);

	  inParams->cudaDeviceId = *std::istream_iterator<int>(iss);
	  if (iss.fail()) {
	    std::cerr << "Could not parse device ID " << argv[i] << std::endl;
	    std::abort();
	  }
	} else {
	  std::cerr << "Device ID required.\n";
	  std::abort();	  
	}
#else 
	std::cerr << "\"-device\" switch is only valid in CUDA enabled builds.\n";
	std::abort();
#endif
      } else if (std::strcmp(argv[i], "-sport") == 0) {
#ifdef _GPU_
         inParams->sportMode = true;
#else // _GPU_
	 std::cerr << "Sorry, GPU execution requires CUDA\n--> using CPU execution." << std::endl;
#endif // _GPU_
      } else if (std::string(argv[i]) == "-precompute") {
	if (i + 1 < argc) {
	  inParams->precompExportFName = argv[++i];
	} else {
	  std::cerr << "-precompute requires a path argument.\n";
	  std::abort();
	}
      } else if (std::string(argv[i]) == "-initialise-with") {
	if (i + 1 < argc) {
	  inParams->precompImportFName = argv[++i];
	} else {
	  std::cerr << "-initialise-with requires a path argument.\n";
	  std::abort();
	}
      } else if (std::strcmp(argv[i], "-print") == 0) {
	if (i + 1 < argc) {
	  char* temp = argv[++i];
	  if (std::strcmp(temp, "NF") == 0)
            inParams->printNForces = 1;
	  else if (std::strcmp(temp, "LNF") == 0)
            inParams->printNDispForces = 1;
	  else if (std::strcmp(temp, "LNFS") == 0)
            inParams->printNDispForcesSums = 1;
	  else if (std::strcmp(temp, "ND") == 0)
            inParams->printNDisps = 1;
	  else if (std::strcmp(temp, "ES") == 0) 
	    inParams->printStrainEnergy = true;
	  else if (std::strcmp(temp, "EK") == 0) 
	    inParams->printKineticEnergy = true;
	  else
	    {
	      std::cerr << "!!!Value " << temp << " invalid for parameter " << argv[i-1] << "." << std::endl;
	      Usage(argv[0]);
	      return 1;
	    }
	} else {
	  std::cerr << "Require a type argument with the \"-print\" switch.\n";
	  std::abort();
	}
      }
      else if (std::strcmp(argv[i], "-plot") == 0)
      {
#ifdef _Visualisation_
         inParams->plotModel = true;
#else // _Visualisation_
	 std::cout << "Sorry, visualisation requires VTK." << std::endl;
#endif // _Visualisation_
      }
      else if (std::strcmp(argv[i], "-t") == 0)
      {
         inParams->doTiming = true;
      }
      else if (std::strcmp(argv[i], "-v") == 0)
      {
         inParams->verbose = true;
      }
      else
      {
         std::cerr << "!!!Parameter " << argv[i] << " unknown." << std::endl;
         Usage(argv[0]);
         return 1;
      }
   }

   // Check inputs
   if (inParams->xmlFName==NULL)
   {
      std::cerr << "XML file not specified!" << std::endl;
      Usage(argv[0]);
      return 1;
   }

   return 0;
}

void PrintIntro()
{
   std::cout << "\n" << std::endl;
   std::cout << "\t*********************************************************" << std::endl;
   std::cout << "\t*\t\t\t\t\t\t\t*" << std::endl;
   std::cout << "\t*\tNifty Sim\t\t\t\t\t*" << std::endl;
   std::cout << "\t*\tNonlinear explicit finite element solver\t*" << std::endl;
   std::cout << "\t*\t\t\t\t\t\t\t*" << std::endl;
   std::cout << "\t*\tZeike Taylor\t\t\t\t\t*" << std::endl;
   std::cout << "\t*\tUniversity of Sheffield\t\t\t\t*\n";
   std::cout << "\t*\tDept. of Mechanical Engineering\t\t\t*" << std::endl;
   std::cout << "\t*\tz.a.taylor@sheffield.ac.uk\t\t\t*" << std::endl;
   std::cout << "\t*\tStian Johnsen\t\t\t\t\t*\n";
   std::cout << "\t*\tUniversity College London\t\t\t*\n";
   std::cout << "\t*\tCentre for Medical Image Computing\t\t*" << std::endl;
   std::cout << "\t*\tstian.johnsen.09@ucl.ac.uk\t\t\t*\n";
   std::cout << "\t*\t\t\t\t\t\t\t*" << std::endl;
   std::cout << "\t*********************************************************" << std::endl;
   std::cout << "\n" << std::endl;
}

void PrintAnalysisSettings(tledModel* Model, InputParameters* inParams)
{
   if (inParams->sportMode)
      std::cout << "\n---- GPU EXECUTION ----\n" << std::endl;
   else
      std::cout << "\n---- CPU EXECUTION ----" << std::endl;
   std::cout << "\nMODEL DATA" << std::endl;
   std::cout << "Model:\t\t" << Model->GetFileName() << std::endl;
   std::cout << "Nodes:\t\t" << Model->GetNumNodes() << "\t(" << Model->GetNumNodes()*3 << " DOF)" << std::endl;
   std::cout << "Elements:\t" << Model->GetNumEls() << std::endl;
   if (Model->GetROM() == true)
   {
      if (!inParams->sportMode)
      {
#ifdef _GPU_
         std::cout << "\nReduced order modelling only available in GPU mode\n--> using GPU execution." << std::endl;
         inParams->sportMode = true;
         std::cout << "\nUSING REDUCED ORDER MODEL" << std::endl;
         std::cout << "(" << Model->GetNumBasisVectors() << " basis vectors)" << std::endl;
#else // _GPU_
         std::cout << "\nSorry, reduced order modelling only available in GPU mode\n--> using full model." << std::endl;
#endif // _GPU_
      }
      else
      {
         std::cout << "\nUSING REDUCED ORDER MODEL" << std::endl;
         std::cout << "(" << Model->GetNumBasisVectors() << " basis vectors)" << std::endl;
      }
   }
   std::cout << std::endl;
}

void PlotMeshes(tledSimulator* Simulator, InputParameters* inParams)
{
#ifdef _Visualisation_
   if (inParams->plotModel == true)
   {
      if (inParams->verbose)
         std::cout << "\nPlotting results" << std::endl;
      
      tledModelViewer* mv = new tledModelViewer(Simulator);
      mv->CreateVTKModels();
      mv->DisplayModel();
      
      delete mv;
   }
#endif // _Visualisation_
}

void PrintResults(tledSolver* Solver, InputParameters* inParams)
{
   if (inParams->printNForces == true)
      Solver->PrintNodalForces();
   if (inParams->printNDispForces == true)
      Solver->PrintDispNodalForces();
   if (inParams->printNDispForcesSums == true)
      Solver->PrintDispNodalForceSums();
   if (inParams->printNDisps == true)
      Solver->PrintNodalDisps();
   if (inParams->printKineticEnergy) 
     Solver->PrintKineticEnergy();
   if (inParams->printStrainEnergy)
     Solver->PrintStrainEnergy();
}


