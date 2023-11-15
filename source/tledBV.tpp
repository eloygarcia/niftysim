// =========================================================================
// File:       tledBV.tpp
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

#ifdef _GPU_

template <>
void tledBV<2>::InitGPU(GPUBV &r_dst);

template <>
void tledBV<4>::InitGPU(GPUBV &r_dst);

#endif
