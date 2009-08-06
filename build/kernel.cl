/*
 * Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and 
 * international Copyright laws.  Users and possessors of this source code 
 * are hereby granted a nonexclusive, royalty-free license to use this code 
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE 
 * OR PERFORMANCE OF THIS SOURCE CODE.  
 *
 * U.S. Government End Users.   This source code is a "commercial item" as 
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of 
 * "commercial computer  software"  and "commercial computer software 
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995) 
 * and is provided to the U.S. Government only as a commercial end item.  
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
 * source code with only those rights set forth herein. 
 *
 * Any use of this source code in individual and commercial software must 
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

#define BLOCKDIM 256

// Macros to simplify shared memory addressing
#define SX(i) sharedPos[i+ get_local_size(0)*get_local_id(1)]
// This macro is only used the multithreadBodies (MT) versions of kernel code below
#define SX_SUM(i,j) sharedPos[i+ get_local_size(0)*j]    // i + blockDimx * j

/**
 * @param Fi current force on atom i (from previous serial calls)
 * @param ai atom position (xyz) and charge (w)
 * @param aj atom position (xyz) and charge (w)
 *
 * @return updated force on atom i (Fi)
 */
float4 interaction(float4 Fi, float4 ai, float4 aj)
{
    float4 r;

    // r_ij  [3 FLOPS]
    r.x = ai.x - aj.x;
    r.y = ai.y - aj.y;
    r.z = ai.z - aj.z;
    r.w = 0;

    // distSqr = dot(r_ij, r_ij) + EPS^2  [6 FLOPS]
    float distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
    distSqr += 0.0025;

    // invDistCube =1/distSqr^(3/2)  [4 FLOPS (2 mul, 1 sqrt, 1 inv)]
    float invDist = rsqrt(distSqr);

    //float distSixth = distSqr * distSqr * distSqr;
    //float invDistCube = 1.0f / sqrt(distSixth);
    
    // s = m_j * invDistCube [1 FLOP]
    float e = 0.5 * 332.0716 * ai.w * aj.w * invDist;

    // a_i =  a_i + s * r_ij [6 FLOPS]
    Fi.x += r.x * e;
    Fi.y += r.y * e;
    Fi.z += r.z * e;
    Fi.w += e;

    return Fi;
}

// This is the "tile_calculation" function from the GPUG3 article.
float4 gravitation(float4 myPos, float4 accel, __local float4 sharedPos[])
{
    
    // The CUDA 1.1 compiler cannot determine that i is not going to 
    // overflow in the loop below.  Therefore if int is used on 64-bit linux 
    // or windows (or long instead of long long on win64), the compiler
    // generates suboptimal code.  Therefore we use long long on win64 and
    // long on everything else. (Workaround for Bug ID 347697)
#ifdef _Win64
    unsigned long long i = 0;
#else
    unsigned long i = 0;
#endif

    // Here we unroll the loop

    // Note that having an unsigned int loop counter and an unsigned
    // long index helps the compiler generate efficient code on 64-bit
    // OSes.  The compiler can't assume the 64-bit index won't overflow
    // so it incurs extra integer operations.  This is a standard issue
    // in porting 32-bit code to 64-bit OSes.

    for (unsigned int counter = 0; counter < get_local_size(0); ) 
    {
        accel = interaction(accel, SX(i++), myPos); 
	counter++;
    }

    return accel;
}


// WRAP is used to force each block to start working on a different 
// chunk (and wrap around back to the beginning of the array) so that
// not all multiprocessors try to read the same memory locations at 
// once.
#define WRAP(x,m) (((x)<m)?(x):(x-m))  // Mod without divide, works on values from 0 up to 2m

float4 computeAtomForce(float4 bodyPos, __global float4* positions, int numBodies, __local float4 sharedPos[])
{
  float4 force = {0.0f, 0.0f, 0.0f, 0.0f};

  int threadIdxx = get_local_id(0);
  int threadIdxy = get_local_id(1);
  int blockIdxx = get_group_id(0);
  int blockIdxy = get_group_id(1);
  int gridDimx = get_num_groups(0);
  int blockDimx = get_local_size(0);
  int blockDimy = get_local_size(1);
  int p = blockDimx;
  int q = blockDimy;
  int n = numBodies;
  int numTiles = n / (p * q);
  
  for (int tile = blockIdxy; tile < numTiles + blockIdxy; tile++) 
  {
    sharedPos[threadIdxx+ blockDimx * threadIdxy] = 
      positions[WRAP(blockIdxx + q * tile + threadIdxy, gridDimx) * p
      + threadIdxx];

    // __syncthreads();
    barrier(CLK_LOCAL_MEM_FENCE);
    // This is the "tile_calculation" function from the GPUG3 article.
    force = gravitation(bodyPos, force, sharedPos);

    // __syncthreads();
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  // When the numBodies / thread block size is < # multiprocessors (16 on G80), the GPU is 
  // underutilized.  For example, with a 256 threads per block and 1024 bodies, there will only 
  // be 4 thread blocks, so the GPU will only be 25% utilized. To improve this, we use multiple 
  // threads per body.  We still can use blocks of 256 threads, but they are arranged in q rows 
  // of p threads each.  Each thread processes 1/q of the forces that affect each body, and then 
  // 1/q of the threads (those with threadIdx.y==0) add up the partial sums from the other 
  // threads for that body.  To enable this, use the "--p=" and "--q=" command line options to 
  // this example. e.g.: "nbody.exe --n=1024 --p=64 --q=4" will use 4 threads per body and 256 
  // threads per block. There will be n/p = 16 blocks, so a G80 GPU will be 100% utilized.

  // We use a bool template parameter to specify when the number of threads per body is greater 
  // than one, so that when it is not we don't have to execute the more complex code required!
  SX_SUM(threadIdxx, threadIdxy).x = force.x;
  SX_SUM(threadIdxx, threadIdxy).y = force.y;
  SX_SUM(threadIdxx, threadIdxy).z = force.z;
  SX_SUM(threadIdxx, threadIdxy).w = force.w;


  barrier(CLK_LOCAL_MEM_FENCE);//__syncthreads();

  // Save the result in global memory for the integration step
  if ( get_local_id(0) == 0) 
  {
    for (int i = 1; i < blockDimy; i++) 
    {
      force.x += SX_SUM(get_local_id(0),i).x;
      force.y += SX_SUM(get_local_id(0),i).y;
      force.z += SX_SUM(get_local_id(0),i).z;
      force.w += SX_SUM(get_local_id(0),i).w;
    }
  }

  return force;
}

__kernel void electrostaticKernel(
    __global float4 *positions, 
    __global float4 *gradients, 
    __global int numAtoms, 
    __local float4 sharedPos[])
{
    int threadIdxx = get_local_id(0);
    int threadIdxy = get_local_id(1);
    int blockIdxx = get_group_id(0);
    int blockIdxy = get_group_id(1);
    int gridDimx = get_num_groups(0);
    int blockDimx = get_local_size(0);
    int blockDimy = get_local_size(1);

    int index = (blockIdxx * blockDimx) + threadIdxx;
    float4 pos = positions[index];   

    float4 force = computeAtomForce(pos, positions, numAtoms, sharedPos);

    gradients[index] = force;
}

