#ifndef _tledDeviceTextures_CU_
#define _tledDeviceTextures_CU_

#include "tledCUDAHelpers.h"

#ifndef TLED_USE_TEXTURE_OBJECTS

    texture<float4,1,cudaReadModeElementType> txDhDx;
    texture<int4,1,cudaReadModeElementType> txEInd;
    texture<float4,1,cudaReadModeElementType> txHG;
    texture<int2,1,cudaReadModeElementType> txNodeMap;
    texture<int2,1,cudaReadModeElementType> txFCds;
    texture<float4,1,cudaReadModeElementType> txUprev;
    texture<float4,1,cudaReadModeElementType> txUcurr;
    texture<float,1,cudaReadModeElementType> txUnextV;
    texture<float4,1,cudaReadModeElementType> txFEl;
    texture<float4,1,cudaReadModeElementType> txCD;
    texture<int4,1,cudaReadModeElementType> txFixMask;
    texture<int4,1,cudaReadModeElementType> txDispMask;
    texture<float4,1,cudaReadModeElementType> txUload;
    texture<float4,1,cudaReadModeElementType> txF4Fext;
    texture<float4,1,cudaReadModeElementType> txVol_MType_K;
    texture<float4,1,cudaReadModeElementType> txElasticParams;
    texture<int2,1,cudaReadModeElementType> txNumProny;
    texture<float2,1,cudaReadModeElementType> txEulerIso;
    texture<float2,1,cudaReadModeElementType> txEulerVol;
    texture<float4,1,cudaReadModeElementType> txStateIso;
    texture<float4,1,cudaReadModeElementType> txStateVol;
    texture<float,1,cudaReadModeElementType> txVa;
    texture<float,1,cudaReadModeElementType> txPa;
    texture<float,1,cudaReadModeElementType> txMfc;
    texture<float4,1,cudaReadModeElementType> txNodeCds;
    texture<int,1,cudaReadModeElementType> txCylSlvMask;
    texture<float4,1,cudaReadModeElementType> txCylOrig;
    texture<float4,1,cudaReadModeElementType> txCylAxis;
    texture<int,1,cudaReadModeElementType> txPrbSlvMask;
    texture<float4,1,cudaReadModeElementType> txPrbOrig;
    texture<float4,1,cudaReadModeElementType> txPrbAxis;
    texture<int,1,cudaReadModeElementType> txPltSlvMask;
    texture<float4,1,cudaReadModeElementType> txPltA;
    texture<float4,1,cudaReadModeElementType> txPltB;
    texture<float4,1,cudaReadModeElementType> txPltC;

#else

    __device__ cudaTextureObject_t txDhDx;
    __device__ cudaTextureObject_t txEInd;
    __device__ cudaTextureObject_t txHG;
    __device__ cudaTextureObject_t txNodeMap;
    __device__ cudaTextureObject_t txFCds;
    __device__ cudaTextureObject_t txUprev;
    __device__ cudaTextureObject_t txUcurr;
    __device__ cudaTextureObject_t txUnextV;
    __device__ cudaTextureObject_t txFEl;
    __device__ cudaTextureObject_t txCD;
    __device__ cudaTextureObject_t txFixMask;
    __device__ cudaTextureObject_t txDispMask;
    __device__ cudaTextureObject_t txUload;
    __device__ cudaTextureObject_t txF4Fext;
    __device__ cudaTextureObject_t txVol_MType_K;
    __device__ cudaTextureObject_t txElasticParams;
    __device__ cudaTextureObject_t txNumProny;
    __device__ cudaTextureObject_t txEulerIso;
    __device__ cudaTextureObject_t txEulerVol;
    __device__ cudaTextureObject_t txStateIso;
    __device__ cudaTextureObject_t txStateVol;
    __device__ cudaTextureObject_t txVa;
    __device__ cudaTextureObject_t txPa;
    __device__ cudaTextureObject_t txMfc;
    __device__ cudaTextureObject_t txNodeCds;
    __device__ cudaTextureObject_t txCylSlvMask;
    __device__ cudaTextureObject_t txCylOrig;
    __device__ cudaTextureObject_t txCylAxis;
    __device__ cudaTextureObject_t txPrbSlvMask;
    __device__ cudaTextureObject_t txPrbOrig;
    __device__ cudaTextureObject_t txPrbAxis;
    __device__ cudaTextureObject_t txPltSlvMask;
    __device__ cudaTextureObject_t txPltA;
    __device__ cudaTextureObject_t txPltB;
    __device__ cudaTextureObject_t txPltC;

#endif // TLED_USE_TEXTURE_OBJECTS

#endif // _tledDeviceTextures_CU_
