#pragma once

#include <cuda_runtime.h>

/** Texture compatibility functions for CUDA texture objects */

#ifdef TLED_USE_TEXTURE_OBJECTS
template <class T>
__device__ inline T tledTexFetch(
    cudaTextureObject_t tex,
    int idx);

template <>
__device__ inline float tledTexFetch<float>(
    cudaTextureObject_t tex,
    int idx)
{
    return tex1Dfetch<float>(tex, idx);
}
template <>
__device__ inline float2 tledTexFetch<float2>(
    cudaTextureObject_t tex,
    int idx)
{
    return tex1Dfetch<float2>(tex, idx);
}
template <>
__device__ inline float4 tledTexFetch<float4>(
    cudaTextureObject_t tex,
    int idx)
{
    return tex1Dfetch<float4>(tex, idx);
}


template <>
__device__ inline int tledTexFetch<int>(
    cudaTextureObject_t tex,
    int idx)
{
    return tex1Dfetch<int>(tex, idx);
}
template <>
__device__ inline int2 tledTexFetch<int2>(
    cudaTextureObject_t tex,
    int idx)
{
    return tex1Dfetch<int2>(tex, idx);
}
template <>
__device__ inline int4 tledTexFetch<int4>(
    cudaTextureObject_t tex,
    int idx)
{
    return tex1Dfetch<int4>(tex, idx);
}
#else

template <typename T>
__device__ inline T tledTexFetch(
    texture<T,1,cudaReadModeElementType> tex,
    int idx)
{
    return tex1Dfetch(tex, idx);
}

#endif

namespace tled {
#ifdef TLED_USE_TEXTURE_OBJECTS

template <class T>
inline cudaTextureObject_t CreateTextureObject(
    const void *devicePtr,
    size_t bytes)
{
    cudaResourceDesc resDesc{};
    // memset(&resDesc, 0, sizeof(resDesc));

    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = const_cast<void*>(devicePtr);
    resDesc.res.linear.sizeInBytes = bytes;
    resDesc.res.linear.desc = cudaCreateChannelDesc<T>();

    // cudaTextureDesc texDesc;
    // memset(&texDesc, 0, sizeof(texDesc));
    cudaTextureDesc texDesc{};
    
    texDesc.readMode = cudaReadModeElementType;

    cudaTextureObject_t texObj = 0;

    // cudaCreateTextureObject(
    //     &texObj,
    //     &resDesc,
    //     &texDesc,
    //     nullptr);

    cudaError_t err = cudaCreateTextureObject(
        &texObj,
        &resDesc,
        &texDesc,
        nullptr);

    if (err != cudaSuccess)
    {
        fprintf(stderr,
                "cudaCreateTextureObject failed: %s\n",
                cudaGetErrorString(err));
        return 0;
    }

    return texObj;
}

inline void DestroyTextureObject(
    cudaTextureObject_t texObj)
{
    if (texObj != 0)
        cudaDestroyTextureObject(texObj);
}
#endif
}

#ifdef TLED_USE_TEXTURE_OBJECTS
    template <typename T>
    inline cudaError_t tledBindTexture(
        cudaTextureObject_t &tex,
        const T *devPtr,
        size_t bytes)
    {
        tex = tled::CreateTextureObject<T>(
            devPtr,
            bytes);

    return (tex != 0)
        ? cudaSuccess
        : cudaErrorUnknown;
        // return cudaSuccess;
    }


    template <typename T>
    inline cudaError_t tledBindTextureSymbol(
        cudaTextureObject_t &deviceSymbol,
        const T *devPtr,
        size_t bytes)
    {
        cudaTextureObject_t hostHandle =
            tled::CreateTextureObject<T>(
                devPtr,
                bytes);

        return cudaMemcpyToSymbol(
            deviceSymbol,
            &hostHandle,
            sizeof(cudaTextureObject_t));
    }
#else
    template <typename T>
    inline cudaError_t tledBindTexture(
        texture<T, 1, cudaReadModeElementType> &tex,
        const T *devPtr,
        size_t bytes)
    {
        return cudaBindTexture(
            nullptr,
            tex,
            devPtr,
            bytes);
    }
#endif
