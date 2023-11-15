#ifndef _tledDeviceConstantMemoryVariables_CU_
#define _tledDeviceConstantMemoryVariables_CU_

__device__ __constant__ int c_NumNodes;
__device__ __constant__ int c_NumEls;
__device__ __constant__ int c_NPE;	// Nodes per element
__device__ __constant__ int2 c_maxNumViscTerms;
__device__ __constant__ float3 c_gamma; // Central difference coeffs (ROM)
__device__ __constant__ int c_numContactCyls; // Number of contact cylinders
__device__ __constant__ int c_numContactPrbs; // Number of contact ultrasound probes
__device__ __constant__ int c_numContactPlts; // Number of contact plates

#endif // _tledDeviceConstantMemoryVariables_CU_
