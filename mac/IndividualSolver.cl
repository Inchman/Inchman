/*
* Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
*
* Please refer to the NVIDIA end user license agreement (EULA) associated
* with this source code for terms and conditions that govern your use of
* this software. Any use, reproduction, disclosure, or distribution of
* this software and related documentation outside the terms of the EULA
* is strictly prohibited.
*
*/

#include "helpers.cl"

// types for keys
/*
typedef uint type_key;
typedef uint2 type_key2;
typedef uint4 type_key4;
*/
typedef ulong type_key;
typedef ulong2 type_key2;
typedef ulong4 type_key4;

// macros to extract index/key values from combined pair
// keys are the ones we sort  .. lowest 32 bits
#define get_key(kvp) (kvp & 0xffffffff)
#define get_index(kvp) (kvp >> 32)
#define KEYS_MULT  4294967296
#define KEY_INDEX_PAIR(key, index) (key + index *KEYS_MULT)
#define KILL_PARTICLE(kvp) (kvp | 0xffffffff)
//----------------------------------------------------------------------------
// Scans each warp in parallel ("warp-scan"), one element per thread.
// uses 2 numElements of shared memory per thread (64 = elements per warp)
//----------------------------------------------------------------------------
#define WARP_SIZE 32
type_key scanwarp(type_key val, volatile __local type_key* sData, int maxlevel)
{
    // The following is the same as 2 * RadixSort::WARP_SIZE * warpId + threadInWarp = 
    // 64*(threadIdx.x >> 5) + (threadIdx.x & (RadixSort::WARP_SIZE - 1))
    int localId = get_local_id(0);
    int idx = 2 * localId - (localId & (WARP_SIZE - 1));
    sData[idx] = 0;
    idx += WARP_SIZE;
    sData[idx] = val;     

    if (0 <= maxlevel) { sData[idx] += sData[idx - 1]; }
    if (1 <= maxlevel) { sData[idx] += sData[idx - 2]; }
    if (2 <= maxlevel) { sData[idx] += sData[idx - 4]; }
    if (3 <= maxlevel) { sData[idx] += sData[idx - 8]; }
    if (4 <= maxlevel) { sData[idx] += sData[idx -16]; }

    return sData[idx] - val;  // convert inclusive -> exclusive
}

//----------------------------------------------------------------------------
// scan4 scans 4*RadixSort::CTA_SIZE numElements in a block (4 per thread), using 
// a warp-scan algorithm
//----------------------------------------------------------------------------
type_key4 scan4(type_key4 idata, __local type_key* ptr)
{    
    
    uint idx = get_local_id(0);

    type_key4 val4 = idata;
    type_key sum[3];
    sum[0] = val4.x;
    sum[1] = val4.y + sum[0];
    sum[2] = val4.z + sum[1];
    
    type_key val = val4.w + sum[2];
    
    val = scanwarp(val, ptr, 4);
    barrier(CLK_LOCAL_MEM_FENCE);

    if ((idx & (WARP_SIZE - 1)) == WARP_SIZE - 1)
    {
        ptr[idx >> 5] = val + val4.w + sum[2];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (idx < WARP_SIZE)
      ptr[idx] = scanwarp(ptr[idx], ptr, 2);
    
    barrier(CLK_LOCAL_MEM_FENCE);

    val += ptr[idx >> 5];

    val4.x = val;
    val4.y = val + sum[0];
    val4.z = val + sum[1];
    val4.w = val + sum[2];

    return val4;
}

type_key4 rank4(type_key4 preds, __local type_key* sMem, __local type_key* numtrue)
{
	int localId = get_local_id(0);
	int localSize = get_local_size(0);

	type_key4 address = scan4(preds, sMem);
	
	if (localId == localSize - 1) 
	{
		numtrue[0] = address.w + preds.w;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	
	type_key4 rank;
	int idx = localId*4;
	rank.x = (preds.x) ? address.x : numtrue[0] + idx - address.x;
	rank.y = (preds.y) ? address.y : numtrue[0] + idx + 1 - address.y;
	rank.z = (preds.z) ? address.z : numtrue[0] + idx + 2 - address.z;
	rank.w = (preds.w) ? address.w : numtrue[0] + idx + 3 - address.w;
	
	return rank;
}

void radixSortBlockKeysOnly(type_key4 *key, uint nbits, uint startbit, __local type_key* sMem, __local type_key* numtrue)
{
  int localId = get_local_id(0);
  int localSize = get_local_size(0);
  
  for(uint shift = startbit; shift < (startbit + nbits); ++shift)
    {
      type_key4 lsb;
      lsb.x = !(((*key).x >> shift) & 0x1);
      lsb.y = !(((*key).y >> shift) & 0x1);
      lsb.z = !(((*key).z >> shift) & 0x1);
      lsb.w = !(((*key).w >> shift) & 0x1);
        
      type_key4 r;
      
      r = rank4(lsb, sMem, numtrue);
      
      // This arithmetic strides the ranks across 4 CTA_SIZE regions
      sMem[(r.x & 3) * localSize + (r.x >> 2)] = (*key).x;
      sMem[(r.y & 3) * localSize + (r.y >> 2)] = (*key).y;
      sMem[(r.z & 3) * localSize + (r.z >> 2)] = (*key).z;
      sMem[(r.w & 3) * localSize + (r.w >> 2)] = (*key).w;
      barrier(CLK_LOCAL_MEM_FENCE);
      
      // The above allows us to read without 4-way bank conflicts:
      (*key).x = sMem[localId];
      (*key).y = sMem[localId +     localSize];
      (*key).z = sMem[localId + 2 * localSize];
      (*key).w = sMem[localId + 3 * localSize];
      
      barrier(CLK_LOCAL_MEM_FENCE);
    }
}

// this is uint4 since each thread works on four keys..
__kernel void radixSortBlocksKeysOnly(__global type_key4* keysIn, 
				      __global type_key4* keysOut,
				      uint nbits,
				      uint startbit,
				      uint numElements, 
				      uint totalBlocks,
				      __local type_key* sMem)
{
	int globalId = get_global_id(0);
	__local type_key numtrue[1];

	type_key4 key;
	key = keysIn[globalId];
	
	barrier(CLK_LOCAL_MEM_FENCE);
	
	radixSortBlockKeysOnly(&key, nbits, startbit, sMem, numtrue);
	
	keysOut[globalId] = key;
}

//----------------------------------------------------------------------------
// Given an array with blocks sorted according to a 4-bit radix group, each 
// block counts the number of keys that fall into each radix in the group, and 
// finds the starting offset of each radix in the block.  It then writes the radix 
// counts to the counters array, and the starting offsets to the blockOffsets array.
//
// Template parameters are used to generate efficient code for various special cases
// For example, we have to handle arrays that are a multiple of the block size 
// (fullBlocks) differently than arrays that are not. "loop" is used when persistent 
// CTAs are used. 
//
// By persistent CTAs we mean that we launch only as many thread blocks as can 
// be resident in the GPU and no more, rather than launching as many threads as
// we have elements. Persistent CTAs loop over blocks of elements until all work
// is complete.  This can be faster in some cases.  In our tests it is faster
// for large sorts (and the threshold is higher on compute version 1.1 and earlier
// GPUs than it is on compute version 1.2 GPUs.
//                                
//----------------------------------------------------------------------------
__kernel void findRadixOffsets(__global type_key2* keys,
			       __global uint* counters,
			       __global uint* blockOffsets,
			       uint startbit,
			       uint numElements,
			       uint totalBlocks,
			       __local type_key* sRadix1)
{
	__local uint  sStartPointers[16];

    uint groupId = get_group_id(0);
    uint localId = get_local_id(0);
    uint groupSize = get_local_size(0);

    type_key2 radix2;

    radix2 = keys[get_global_id(0)];
        

    sRadix1[2 * localId]     = (radix2.x >> startbit) & 0xF;
    sRadix1[2 * localId + 1] = (radix2.y >> startbit) & 0xF;

    // Finds the position where the sRadix1 entries differ and stores start 
    // index for each radix.
    if(localId < 16) 
    {
        sStartPointers[localId] = 0; 
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if((localId > 0) && (sRadix1[localId] != sRadix1[localId - 1]) ) 
    {
        sStartPointers[sRadix1[localId]] = localId;
    }
    if(sRadix1[localId + groupSize] != sRadix1[localId + groupSize - 1]) 
    {
        sStartPointers[sRadix1[localId + groupSize]] = localId + groupSize;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if(localId < 16) 
    {
        blockOffsets[groupId*16 + localId] = sStartPointers[localId];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Compute the sizes of each block.
    if((localId > 0) && (sRadix1[localId] != sRadix1[localId - 1]) ) 
    {
        sStartPointers[sRadix1[localId - 1]] = 
            localId - sStartPointers[sRadix1[localId - 1]];
    }
    if(sRadix1[localId + groupSize] != sRadix1[localId + groupSize - 1] ) 
    {
        sStartPointers[sRadix1[localId + groupSize - 1]] = 
            localId + groupSize - sStartPointers[sRadix1[localId + groupSize - 1]];
    }
        

    if(localId == groupSize - 1) 
    {
        sStartPointers[sRadix1[2 * groupSize - 1]] = 
            2 * groupSize - sStartPointers[sRadix1[2 * groupSize - 1]];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if(localId < 16) 
    {
        counters[localId * totalBlocks + groupId] = sStartPointers[localId];
    }
}

// a naive scan routine that works only for array that
// can fit into a single block, just for debugging purpose,
// not used in the sort now
__kernel void scanNaive(__global uint *g_odata, 
                        __global uint *g_idata, 
                        uint n,
                        __local uint* temp)
{

    int localId = get_local_id(0);

    int pout = 0;
    int pin = 1;

    // Cache the computational window in shared memory
    temp[pout*n + localId] = (localId > 0) ? g_idata[localId-1] : 0;

    for (int offset = 1; offset < n; offset *= 2)
    {
        pout = 1 - pout;
        pin  = 1 - pout;
        barrier(CLK_LOCAL_MEM_FENCE);

        temp[pout*n+localId] = temp[pin*n+localId];

        if (localId >= offset)
            temp[pout*n+localId] += temp[pin*n+localId - offset];
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    g_odata[localId] = temp[pout*n+localId];
}

//----------------------------------------------------------------------------
// reorderData shuffles data in the array globally after the radix offsets 
// have been found. On compute version 1.1 and earlier GPUs, this code depends 
// on RadixSort::CTA_SIZE being 16 * number of radices (i.e. 16 * 2^nbits).
// 
// On compute version 1.1 GPUs ("manualCoalesce=true") this function ensures
// that all writes are coalesced using extra work in the kernel.  On later
// GPUs coalescing rules have been relaxed, so this extra overhead hurts 
// performance.  On these GPUs we set manualCoalesce=false and directly store
// the results.
//
// Template parameters are used to generate efficient code for various special cases
// For example, we have to handle arrays that are a multiple of the block size 
// (fullBlocks) differently than arrays that are not.  "loop" is used when persistent 
// CTAs are used. 
//
// By persistent CTAs we mean that we launch only as many thread blocks as can 
// be resident in the GPU and no more, rather than launching as many threads as
// we have elements. Persistent CTAs loop over blocks of elements until all work
// is complete.  This can be faster in some cases.  In our tests it is faster
// for large sorts (and the threshold is higher on compute version 1.1 and earlier
// GPUs than it is on compute version 1.2 GPUs.
//----------------------------------------------------------------------------
__kernel void reorderDataKeysOnly(__global type_key  *outKeys, 
                                  __global type_key2  *keys, 
                                  __global uint  *blockOffsets, 
                                  __global uint  *offsets, 
                                  __global uint  *sizes, 
                                  uint startbit,
                                  uint numElements,
                                  uint totalBlocks,
                                  __local type_key2* sKeys2)
{
    __local uint sOffsets[16];
    __local uint sBlockOffsets[16];

    __local type_key *sKeys1 = (__local type_key*)sKeys2; 

    uint groupId = get_group_id(0);

    uint globalId = get_global_id(0);
    uint localId = get_local_id(0);
    uint groupSize = get_local_size(0);

    sKeys2[localId]   = keys[globalId];

    if(localId < 16)  
    {
        sOffsets[localId]      = offsets[localId * totalBlocks + groupId];
        sBlockOffsets[localId] = blockOffsets[groupId * 16 + localId];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // the radices probably can stay uint..
    type_key radix = (sKeys1[localId] >> startbit) & 0xF;
    uint globalOffset = sOffsets[radix] + localId - sBlockOffsets[radix];

    if (globalOffset < numElements)
    {
        outKeys[globalOffset]   = sKeys1[localId];
    }

    radix = (sKeys1[localId + groupSize] >> startbit) & 0xF;
    globalOffset = sOffsets[radix] + localId + groupSize - sBlockOffsets[radix];

    if (globalOffset < numElements)
    {
        outKeys[globalOffset]   = sKeys1[localId + groupSize];
    }
}

/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */



//All three kernels run 512 threads per workgroup
//Must be a power of two

////////////////////////////////////////////////////////////////////////////////
// Scan codelets
////////////////////////////////////////////////////////////////////////////////
#if(1)
    //Naive inclusive scan: O(N * log2(N)) operations
    //Allocate 2 * 'size' local memory, initialize the first half
    //with 'size' zeros avoiding if(pos >= offset) condition evaluation
    //and saving instructions
    inline uint scan1Inclusive(uint idata, __local uint *l_Data, uint size){
        uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
        l_Data[pos] = 0;
        pos += size;
        l_Data[pos] = idata;

        for(uint offset = 1; offset < size; offset <<= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            uint t = l_Data[pos] + l_Data[pos - offset];
            barrier(CLK_LOCAL_MEM_FENCE);
            l_Data[pos] = t;
        }

        return l_Data[pos];
    }

    inline uint scan1Exclusive(uint idata, __local uint *l_Data, uint size){
        return scan1Inclusive(idata, l_Data, size) - idata;
    }

#else
    #define LOG2_WARP_SIZE 5U
    #define      WARP_SIZE (1U << LOG2_WARP_SIZE)

    //Almost the same as naiveScan1 but doesn't need barriers
    //assuming size <= WARP_SIZE
    inline uint warpScanInclusive(uint idata, __local uint *l_Data, uint size){
        uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
        l_Data[pos] = 0;
        pos += size;
        l_Data[pos] = idata;

        for(uint offset = 1; offset < size; offset <<= 1)
            l_Data[pos] += l_Data[pos - offset];

        return l_Data[pos];
    }

    inline uint warpScanExclusive(uint idata, __local uint *l_Data, uint size){
        return warpScanInclusive(idata, l_Data, size) - idata;
    }

    inline uint scan1Inclusive(uint idata, __local uint *l_Data, uint size){
        if(size > WARP_SIZE){
            //Bottom-level inclusive warp scan
            uint warpResult = warpScanInclusive(idata, l_Data, WARP_SIZE);

            //Save top elements of each warp for exclusive warp scan
            //sync to wait for warp scans to complete (because l_Data is being overwritten)
            barrier(CLK_LOCAL_MEM_FENCE);
            if( (get_local_id(0) & (WARP_SIZE - 1)) == (WARP_SIZE - 1) )
                l_Data[get_local_id(0) >> LOG2_WARP_SIZE] = warpResult;

            //wait for warp scans to complete
            barrier(CLK_LOCAL_MEM_FENCE);
            if( get_local_id(0) < (WORKGROUP_SIZE_SCAN / WARP_SIZE) ){
                //grab top warp elements
                uint val = l_Data[get_local_id(0)];
                //calculate exclsive scan and write back to shared memory
                l_Data[get_local_id(0)] = warpScanExclusive(val, l_Data, size >> LOG2_WARP_SIZE);
            }

            //return updated warp scans with exclusive scan results
            barrier(CLK_LOCAL_MEM_FENCE);
            return warpResult + l_Data[get_local_id(0) >> LOG2_WARP_SIZE];
        }else{
            return warpScanInclusive(idata, l_Data, size);
        }
    }

    inline uint scan1Exclusive(uint idata, __local uint *l_Data, uint size){
        return scan1Inclusive(idata, l_Data, size) - idata;
    }
#endif


//Vector scan: the array to be scanned is stored
//in work-item private memory as uint4
inline uint4 scan4Inclusive(uint4 data4, __local uint *l_Data, uint size){
    //Level-0 inclusive scan
    data4.y += data4.x;
    data4.z += data4.y;
    data4.w += data4.z;

    //Level-1 exclusive scan
    uint val = scan1Inclusive(data4.w, l_Data, size / 4) - data4.w;

    return (data4 + (uint4)val);
}

inline uint4 scan4Exclusive(uint4 data4, __local uint *l_Data, uint size){
    return scan4Inclusive(data4, l_Data, size) - data4;
}


////////////////////////////////////////////////////////////////////////////////
// Scan kernels
////////////////////////////////////////////////////////////////////////////////
__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE_SCAN, 1, 1)))
void scanExclusiveLocal1(
    __global uint4 *d_Dst,
    __global uint4 *d_Src,
    uint size,
    __local uint* l_Data
){
    //Load data
    uint4 idata4 = d_Src[get_global_id(0)];

    //Calculate exclusive scan
    uint4 odata4  = scan4Exclusive(idata4, l_Data, size);

    //Write back
    d_Dst[get_global_id(0)] = odata4;
}

//Exclusive scan of top elements of bottom-level scans (4 * THREADBLOCK_SIZE)
__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE_SCAN, 1, 1)))
void scanExclusiveLocal2(
    __global uint *d_Buf,
    __global uint *d_Dst,
    __global uint *d_Src,
    uint N,
    uint arrayLength,
    __local uint* l_Data
){
    //Load top elements
    //Convert results of bottom-level scan back to inclusive
    //Skip loads and stores for inactive work-items of the work-group with highest index(pos >= N)
    uint data = 0;
    if(get_global_id(0) < N)
    data =
        d_Dst[(4 * WORKGROUP_SIZE_SCAN - 1) + (4 * WORKGROUP_SIZE_SCAN) * get_global_id(0)] + 
        d_Src[(4 * WORKGROUP_SIZE_SCAN - 1) + (4 * WORKGROUP_SIZE_SCAN) * get_global_id(0)];

    //Compute
    uint odata = scan1Exclusive(data, l_Data, arrayLength);

    //Avoid out-of-bound access
    if(get_global_id(0) < N)
        d_Buf[get_global_id(0)] = odata;
}

//Final step of large-array scan: combine basic inclusive scan with exclusive scan of top elements of input arrays
__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE_SCAN, 1, 1)))
void uniformUpdate(
    __global uint4 *d_Data,
    __global uint *d_Buf
){
    __local uint buf[1];

    uint4 data4 = d_Data[get_global_id(0)];

    if(get_local_id(0) == 0)
        buf[0] = d_Buf[get_group_id(0)];

    barrier(CLK_LOCAL_MEM_FENCE);
    data4 += (uint4)buf[0];
    d_Data[get_global_id(0)] = data4;
}

// ----------------------------------------------------------------------------------------------------------------
__kernel void padKeys(__global type_key *d_keys, uint nkeys) {
  uint tid = get_global_id(0);

  // pad by setting the cell number to max .. will keep the id though!
  d_keys[nkeys-1-tid] = KILL_PARTICLE(d_keys[nkeys-1-tid]);
} // pad keys

// This one computes the boundaries for each cell [the cell number is the key value..]
__kernel void computeCellIndices(__global type_key *d_keys,
				 __global uint *d_cellIndicesLeft,
				 __global uint *d_cellIndicesRight,
				 uint nkeys, uint ncells)
{  
  uint tid = get_global_id(0);


  // read in keys  
  uint myKey = get_key(d_keys[tid]);

  
  // write first values
  // we need to pad with zeros until we actually have an occupied cell
  if (tid == 0) {
    uint i;
    for (i=0; i<myKey; i++)
      {
	d_cellIndicesLeft[i]=0;
	d_cellIndicesRight[i]=0;
      }
    d_cellIndicesLeft[myKey]=0;
  }
  
  
  // check right neighbour
  if (tid < (nkeys-1)) {
    uint rightKey = get_key(d_keys[tid+1]);
    
    if (rightKey > myKey) {
      // set cell boundaries
      d_cellIndicesRight[myKey] = tid;
      uint i;
      
      // and fill all indices up (up to rightKey)
      for (i=1; i< (rightKey-myKey); i++) {       
	d_cellIndicesLeft[myKey+i] = tid+1;
	d_cellIndicesRight[myKey+i] = tid+1;
      }
      d_cellIndicesLeft[rightKey] = tid+1;      
    }
  }
  
  // do final
  if (tid==(nkeys-1)) {
    d_cellIndicesRight[myKey] = tid;

    uint i;
    // we fill up the last ones with pointing to the next particle
    // (this is needed for the reaction routines..)
    for (i=1; i<(ncells-myKey); i++) {
      d_cellIndicesLeft[myKey+i] = tid+1;
      d_cellIndicesRight[myKey+i] = tid+1;
    }      
  }    
} // compute cell indices


// computes the diffusivity/drift
// <<<! computeDiffusivityDrift !>>>

/*inline void computeDiffusivityDrift(_species_t *actualSpecies, int speciesIndex)
{
    // define the species structs
    _species_t A = {0,0,0,0,0}; // will be filled in later for each individual

    _species_t *_allSpecies[1]; // contains list of all species
    _allSpecies[0] = &A;

    // fill in values for this particle
    _allSpecies[speciesIndex] = actualSpecies;

    // compute diffusivity/drift from user-defined method
    #define MeanX _allSpecies[speciesIndex]->MeanX
    #define DiffusivityX _allSpecies[speciesIndex]->DiffusivityX
    #define DiffusivityY _allSpecies[speciesIndex]->DiffusivityY
    #define DriftX _allSpecies[speciesIndex]->DriftX
    #define DriftY _allSpecies[speciesIndex]->DriftY
    DiffusivityX = 1.;
    DiffusivityY = 2.;
    DriftX = 1.5;
    DriftY = 3.;
    #undef MeanX
    #undef DiffusivityX
    #undef DiffusivityY
    #undef DriftX
    #undef DriftY
}*/


__kernel void computeDiffusionConstants(__global type_key *d_keys,
					__global uint *d_cellIndicesLeft,
					__global uint *d_cellIndicesRight,
					__global Real *d_properties,
					const uint propertiesOffset,
					__global Real *d_diffusionConstantsX,
					__global Real *d_diffusionConstantsY,
					__global Real *d_rx,
					__global Real *d_ry,
					__global uint *d_state,
					const Real PhysicalSimTime /* exposed to user */,
					const int speciesIndex,
					__global Real *d_sumMoments,
					__global Real *d_sumStates)
{
  const size_t _cellIndex = get_global_id_2d;
  const size_t _sIndex = get_global_id_for_species_2d(speciesIndex);

  Real _diffx = 0.;
  Real _diffy = 0.;
  Real _rx = 0.;
  Real _ry = 0.;

  uint _particleCount = 0.; // to store in state

  // go through list of species for this cell
  for (uint i=d_cellIndicesLeft[_cellIndex]; i<=d_cellIndicesRight[_cellIndex]; i++) {
    // fetch key-value pair
    type_key kvp = d_keys[i];

    // check if particle actually belongs to my cell
    if (get_key(kvp) == _cellIndex) {
      _particleCount++;

      uint index = get_index(kvp);

      // get the diffusivity/drift from the user-defined method
      /*_species_t _actualSpecies = {0,
                                  d_properties[index + 0 * propertiesOffset],
                                  d_properties[index + 1 * propertiesOffset],
                                  d_properties[index + 2 * propertiesOffset],
                                  d_properties[index + 3 * propertiesOffset]};*/
      _species_t _actualSpecies = {
                                   d_state[get_global_id_for_species_2d(speciesIndex)],
                                   d_state+get_global_area_2d*speciesIndex,
                                   0,
                                   d_properties[index + 0 * propertiesOffset],
                                   d_properties[index + 1 * propertiesOffset],
                                   d_properties[index + 2 * propertiesOffset],
                                   d_properties[index + 3 * propertiesOffset]
                                   };

      computeDiffusivityDrift(&_actualSpecies, speciesIndex);

      // check if drift/diffusion is higher than any
      _diffx = fmax(_actualSpecies.DiffusivityX, _diffx);
      _diffy = fmax(_actualSpecies.DiffusivityY, _diffy);
      _rx    = maxmag(_actualSpecies.DriftX, _rx);
      _ry    = maxmag(_actualSpecies.DriftY, _ry);

      /*
      if (d_properties[index + 0 * propertiesOffset] > _diffx)
        _diffx = d_properties[index + 0 * propertiesOffset];
      if (d_properties[index + 1 * propertiesOffset] > _diffy)
        _diffy = d_properties[index + 1 * propertiesOffset];
      if (fabs(d_properties[index + 2 * propertiesOffset]) > fabs(_rx))
        _rx = d_properties[index + 2 * propertiesOffset];
      if (fabs(d_properties[index + 3 * propertiesOffset]) > fabs(_ry))
        _ry = d_properties[index + 3 * propertiesOffset];
        */
    } // is in cell
  } // go through species

  // write values
  /*
  d_diffusionConstantsX[_cellIndex] = diffx;
  d_rx[_cellIndex] = rx;
  d_diffusionConstantsY[_cellIndex] = diffy;
  d_ry[_cellIndex] = ry;
  d_state[get_global_id_for_species_2d(speciesIndex)] = particleCount;
  */
  // todo : remove me
  //diffx = 1.;
  //diffy = 1.;
  //rx = 0.;
  //ry = 0.;
  // remove me end

  d_diffusionConstantsX[_sIndex] = _diffx;
  d_rx[_sIndex] = _rx;
  d_diffusionConstantsY[_sIndex] = _diffy;
  d_ry[_sIndex] = _ry;
  d_state[_sIndex] = _particleCount;
}

// Kernel for individual diffusion
__kernel void individual_diffuseX(__global type_key *d_keys,
				  __global uint *d_cellIndicesLeft,
				  __global uint *d_cellIndicesRight,
				  __global Real *d_properties,
                                  __global uint *d_state,
                                  const uint propertiesOffset,
				  int seed,
				  int speciesIndex, Real dt, __global Real *d_errors
)
{
    const Real deltax = PhysicalCellWidth;
    const Real deltay = PhysicalCellHeight;

    // get absolute position on the device
    const size_t device_x = get_global_id(0);
    const size_t device_y = get_global_id(1);   
    const size_t global_id = get_global_id_for_species_2d(speciesIndex);

    // Initialize counters and keys for RNG
    threefry2x32_key_t k = {{global_id, 0xdecafbad+seed}};
    threefry2x32_ctr_t c = {{0, 0xf00dcafe}};
    union {
        threefry2x32_ctr_t c;
        int2 r;
    } u;

    // initialize RNG counter - note: we can't seed here!
    c.v[0]=0;

    // this is a mask - decides which directions we can go from here
    // we use an int object, where 0 means -x, 1 means +x, 2 means
    // -y and 3 means +y
    size_t directionMask[4];

    // decide if we're at an outer boundary where we can't go everywhere
    // boundaryMask for that direction needs to be:
    // 0 (reflective boundary or source) TODO check out the source boundary rule
    // 1 (absorbing boundary, i.e. particle sink)
    // TODO maybe add extra kernel for every BC to avoid conditionals
    directionMask[0] = (device_x == 0) ? BoundaryMaskX1 : 1;
    directionMask[1] = (device_x == GridModelWidth-1) ? BoundaryMaskX2 : 1;
    directionMask[2] = (device_y == 0) ? BoundaryMaskY1 : 1;
    directionMask[3] = (device_y == GridModelHeight-1) ? BoundaryMaskY2 : 1;    

    // main loop. For every particle decide if it stays in the cell or moves.
    // go through list of species for this cell
    const size_t _cellIndex = get_global_id_2d;

    for (uint i=d_cellIndicesLeft[_cellIndex]; i<=d_cellIndicesRight[_cellIndex]; i++) {
      // fetch key-value pair
      type_key kvp = d_keys[i];
      // check if particle actually belongs to my cell
      // we still need that check .. since if particles were killed in this cell they're still floating around here until the list is re-sorted!
      if (get_key(kvp) == _cellIndex) {
	uint key = get_key(kvp);
	uint index = get_index(kvp);

	// values for the direction picking
	// we need only two since we do two independent sweeps!
	// p1 is for right/up, p2 for left/down
	Real p1, p2, ps;

	// fetch drift field and diffusion constant
	Real diffHere, rx, dx;

        // get the diffusivity/drift from the user-defined method
        /*_species_t _actualSpecies = {0,
                                    d_properties[index + 0 * propertiesOffset],
                                    d_properties[index + 1 * propertiesOffset],
                                    d_properties[index + 2 * propertiesOffset],
                                    d_properties[index + 3 * propertiesOffset]};*/
        _species_t _actualSpecies = {
            d_state[get_global_id_for_species_2d(speciesIndex)],
            d_state+get_global_area_2d*speciesIndex,
            0,
            d_properties[index + 0 * propertiesOffset],
            d_properties[index + 1 * propertiesOffset],
            d_properties[index + 2 * propertiesOffset],
            d_properties[index + 3 * propertiesOffset]
        };
        computeDiffusivityDrift(&_actualSpecies, speciesIndex);
        diffHere = _actualSpecies.DiffusivityX;
        rx = _actualSpecies.DriftX;
        //diffHere = d_properties[index + 0 * propertiesOffset];
        //rx = d_properties[index + 2 * propertiesOffset];
	dx = deltax;

        // todo: remove me
        //diffHere = 1.;
        //rx = 0.;
        // end remove me

	// compute probabilities
	ps = 1. - (dt * (2.*diffHere + rx*rx*dt))/(dx*dx);
	p1 = (1.-ps)*(1. + (dx*rx)/(2.*diffHere + rx*rx*dt))/2.;
	p2 = p1 + (1.-ps)*(1. - (dx*rx)/(2.*diffHere + rx*rx*dt))/2.;
	
        // draw two random numbers
        // We only need to draw every second round since we only need one number
        // shouldn't give any divergence I think.. but check!
        float rand;
	/*
        if (! (i%2)) {
            c.v[0]++;
            u.c = threefry2x32(c, k);
            rand = u.r.x*M_RAN_INVM32+0.5;
        } else {
            rand = (u.r.y*M_RAN_INVM32+0.5);
	    }*/
	c.v[0]++;
	u.c = threefry2x32(c, k);
	rand = u.r.x*M_RAN_INVM32+0.5;

        // which direction?
        size_t dir=0;
        dir += (rand<=p1);

	// now dir is 0 (-x) or 1 (+x)
	// compute new key (= cell position)
	key -= directionMask[dir] * (dir==0) * (rand<=p2);
	key += directionMask[dir] * (dir==1) * (rand<=p2);

#ifdef GPGMP_BALLISTIC
	// ballistic scattering
	if ((rand<=p2) && (directionMask[dir]==0)) d_properties[index + 2 * propertiesOffset] = -rx;
#endif // ballistic scattering
	  
	// and store it
	d_keys[i] = KEY_INDEX_PAIR(key, index);
      } // if it's part of my cell
    } // loop over particles
} // diffuseX

// Kernel for individual diffusion
__kernel void individual_diffuseY(__global type_key *d_keys,
				  __global uint *d_cellIndicesLeft,
				  __global uint *d_cellIndicesRight,
				  __global Real *d_properties,
                                  __global uint *d_state,
                                  const uint propertiesOffset,
				  int seed,
				  int speciesIndex, Real dt, __global Real *d_errors)
{
    const Real deltax = PhysicalCellWidth;
    const Real deltay = PhysicalCellHeight;

    // get absolute position on the device
    const size_t device_x = get_global_id(0);
    const size_t device_y = get_global_id(1);   
    const size_t global_id = get_global_id_for_species_2d(speciesIndex);

    // Initialize counters and keys for RNG
    threefry2x32_key_t k = {{global_id, 0xdecafbad+seed}};
    threefry2x32_ctr_t c = {{0, 0xf00dcafe}};
    union {
        threefry2x32_ctr_t c;
        int2 r;
    } u;

    // initialize RNG counter - note: we can't seed here!
    c.v[0]=0;

    // this is a mask - decides which directions we can go from here
    // we use an int object, where 0 means -x, 1 means +x, 2 means
    // -y and 3 means +y
    // todo: should pre-compute the direction mask and store as constant array!!
    size_t directionMask[4];

    // decide if we're at an outer boundary where we can't go everywhere
    // boundaryMask for that direction needs to be:
    // 0 (reflective boundary or source) TODO check out the source boundary rule
    // 1 (absorbing boundary, i.e. particle sink)
    // TODO maybe add extra kernel for every BC to avoid conditionals
    directionMask[0] = (device_x == 0) ? BoundaryMaskX1 : 1;
    directionMask[1] = (device_x == GridModelWidth-1) ? BoundaryMaskX2 : 1;
    directionMask[2] = (device_y == 0) ? BoundaryMaskY1 : 1;
    directionMask[3] = (device_y == GridModelHeight-1) ? BoundaryMaskY2 : 1;    

    // main loop. For every particle decide if it stays in the cell or moves.
    // go through list of species for this cell
    const size_t _cellIndex = get_global_id_2d;
    for (uint i=d_cellIndicesLeft[_cellIndex]; i<=d_cellIndicesRight[_cellIndex]; i++) {
      // fetch key-value pair
      type_key kvp = d_keys[i];
      // check if particle actually belongs to my cell
      // we still need that check .. since if particles were killed in this cell they're still floating around here until the list is re-sorted!
      if (get_key(kvp) == _cellIndex) {
	uint key = get_key(kvp);
	uint index = get_index(kvp);

	// values for the direction picking
	// we need only two since we do two independent sweeps!
	// p1 is for right/up, p2 for left/down
	Real p1, p2, ps;

	// fetch drift field and diffusion constant for each particle
	Real diffHere, rx, dx;
        // get the diffusivity/drift from the user-defined method
        /*_species_t _actualSpecies = {0,
                                    d_properties[index + 0 * propertiesOffset],
                                    d_properties[index + 1 * propertiesOffset],
                                    d_properties[index + 2 * propertiesOffset],
                                    d_properties[index + 3 * propertiesOffset]};*/
        _species_t _actualSpecies = {
            d_state[get_global_id_for_species_2d(speciesIndex)],
            d_state+get_global_area_2d*speciesIndex,
            0,
            d_properties[index + 0 * propertiesOffset],
            d_properties[index + 1 * propertiesOffset],
            d_properties[index + 2 * propertiesOffset],
            d_properties[index + 3 * propertiesOffset]
        };
        computeDiffusivityDrift(&_actualSpecies, speciesIndex);
        diffHere = _actualSpecies.DiffusivityY;
        rx = _actualSpecies.DriftY;

        //diffHere = d_properties[index + 1 * propertiesOffset];
        //rx = d_properties[index + 3 * propertiesOffset];
	dx = deltax;

        // todo: remove me
        //diffHere = 1.;
        //rx = 0.;
        // end remove me

	// compute probabilities
	ps = 1. - (dt * (2.*diffHere + rx*rx*dt))/(dx*dx);
	p1 = (1.-ps)*(1. + (dx*rx)/(2.*diffHere + rx*rx*dt))/2.;
	p2 = p1 + (1.-ps)*(1. - (dx*rx)/(2.*diffHere + rx*rx*dt))/2.;

        // draw two random numbers
        // We only need to draw every second round since we only need one number
        // shouldn't give any divergence I think.. but check!
        float rand;
	/*
        if (! (i%2)) {
            c.v[0]++;
            u.c = threefry2x32(c, k);
            rand = u.r.x*M_RAN_INVM32+0.5;
        } else {
            rand = (u.r.y*M_RAN_INVM32+0.5);
	    }*/
	c.v[0]++;
	u.c = threefry2x32(c, k);
	rand = u.r.x*M_RAN_INVM32+0.5;

        // which direction?
        size_t dir=2;
        dir += (rand<=p1);

	// now dir is 2 (-y) or 3 (+y)
	// compute new key (= cell position)
	key -= GridModelWidth * directionMask[dir] * (dir==2) * (rand<=p2);
	key += GridModelWidth * directionMask[dir] * (dir==3) * (rand<=p2);

#ifdef GPGMP_BALLISTIC
	// ballistic scattering
	if ((rand<=p2) && (directionMask[dir]==0)) d_properties[index + 3 * propertiesOffset] = -rx;
#endif // ballistic scattering

	// and store it
	d_keys[i] = KEY_INDEX_PAIR(key, index);
      } // check if it's part of my cell
    } // loop over particles
} // diffuseY


__kernel void maintainPopulation(__global type_key *d_keys,
				 __global uint *d_cellIndicesLeft,
				 __global uint *d_cellIndicesRight,
				 __global int *d_reductionOffsets,
				 const int speciesIndex,
				 __global uint *d_state,
				 int seed,
				 __global uint *d_newParticles,
				 __global uint *d_domainReductionOffset,
				 __global uchar *d_flags)
{
  
  const size_t _cellIndex = get_global_id_2d;

  // first compute the is-state and compare it to goal (as determined by state array)
  // todo: it would be better to sacrifice some device memory and store the particles
  // of each individual species that are "lost" or "gained" in the Gillespie step..
  // this way we just had to look them up here. Would be quite memory intensive though
  // todo: we need to do this anyway if we produce new particles..
  // todo: we can do all of this in the gillespie kernel ..
  const uint goalState = d_state[get_global_id_for_species_2d(speciesIndex)];
  uint isState = 0;

  for (uint i= d_cellIndicesLeft[_cellIndex]; i<=d_cellIndicesRight[_cellIndex]; i++) {
    // fetch key-value pair
    type_key kvp = d_keys[i];

    // check if particle actually belongs to my cell
    if (get_key(kvp) == _cellIndex) {
      d_reductionOffsets[i] = 0; // 0 in reduction means keep particle
      isState++;
    }
  }
  
  // Initialize counters and keys for RNG
  threefry2x32_key_t k = {{_cellIndex, 0xdecafbad+seed}};
  threefry2x32_ctr_t c = {{0, 0xf00dcafe}};
  union {
    threefry2x32_ctr_t c;
    int2 r;
  } u;
  
  // initialize RNG counter - note: we can't seed here!
  c.v[0]=0;

  // Now do we need to kill particles or generate new ones?
  if (isState > goalState) {
    // need to kill
    uint ri = 0;
    for (uint i=0; i<(isState-goalState); i++) {

      // randomly pick particles until we find one which isn't dead yet
      uint rand;
      type_key kvp;
      uint kp;
      do {
	if (! (ri%2)) {
	  c.v[0]++;
	  u.c = threefry2x32(c, k);
	  rand = (uint) ((u.r.x*M_RAN_INVM32+0.5)*isState); // todo: don't convert int->float->int
	} else {
	  rand = (uint) ((u.r.y*M_RAN_INVM32+0.5)*isState);
	}
	ri++;
	kp = d_cellIndicesLeft[_cellIndex]+rand;
	kvp = d_keys[kp];
      } while (get_key(kvp) == UINT_MAX);

      // we can now kill the particle - which is done by setting its cell index to max
      // we keep the particle id though .. so we can re-use it later!
      d_keys[kp] = KILL_PARTICLE(kvp);

      // and update the reduction offset
      d_reductionOffsets[kp] = 1; // 1 in reduction means throw particle
    }
  } // kill particles

  // Do we need to add new particles?
  d_newParticles[get_global_id_2d] = 0;
  d_domainReductionOffset[get_global_id_2d] = 0;
  // d_flags[_cellIndex] = (isState>0) ? 1 : 0;

  // We set the flag in the next cell
  if (_cellIndex == 0) d_flags[_cellIndex]=0;
  if (_cellIndex < get_global_area_2d-1)
    d_flags[_cellIndex+1] = (isState>0) ? 1 : 0;
    
  if (isState < goalState) {
  //if (true) {
    // if that's the case we first need to save how many we need for this cell
    d_newParticles[get_global_id_2d] = goalState-isState;

    // now we need to make space for them in the particles list
    // if all cells had at least one particle we could just add this to the
    // reduction offset for later. However, there might be cells with no 
    // particles.. so we need to determine how many empty cells before the 
    // occupied one need new particles. Later the already occupied cell makes
    // space in the particle list before its own particles. We need
    // a segmented scan for this.
    d_domainReductionOffset[get_global_id_2d] = goalState-isState;

    // if there is already particles in this cell we set the flag for the
    // segmented scan to reset the counter
  }

  /*
  // Do we need to add new particles?
  if ( (isState < goalState) && isState>0) {
    // add it to the first member of cell
    d_reductionOffsets[kp] += (goalState-isState);
  } else {
    // TODO: what do we do now???
  }
  */
} // maintain population of individuals

// This kernel copies pool particles from the end of the old particle
// list to the new one at the various cell locations.
__kernel void updateNewParticles(__global uint *d_newParticles,
				 __global uint *d_domainReductionOffset,
				 __global uint *d_cellIndicesLeft,
				 __global type_key *d_keysOld,
				 __global type_key *d_keysNew,
				 __global int *d_reductionOffsetsAdding,
				 __global int *d_reductionOffsets,
                                 uint nparticles, uint npNew,
				 __global Real *d_properties,
                                 __global uint *d_state,
				 const uint propertiesOffset,
				 const int speciesIndex,
				 int seed
) {

    const size_t _cellIndex = get_global_id_2d;

    // Initialize counters and keys for RNG
    threefry2x32_key_t k = {{_cellIndex, 0xdecafbad+seed}};
    threefry2x32_ctr_t c = {{0, 0xf00dcafe}};
    uint ri = 0.;

    // initialize RNG counter - note: we can't seed here!
    c.v[0]=0;

    union {
      threefry2x32_ctr_t c;
      int2 r;
    } u;
    
    // todo: we don't need to do all these calculations each time!
    for (uint i=0; i<d_newParticles[_cellIndex]; i++) {
      // We now need to find the index in the new list where we can write the particles
      size_t indexOld = d_cellIndicesLeft[_cellIndex]; // this is the index in the old list

      // this is the index in the new particle list where the new particle will go
      size_t indexNew;
      
      // this is the id of the newly created particle
      int id;

      // we need to check if the cell index points to a valid particle or
      // to the first invalid particle of the list
      if (indexOld < nparticles) {

	indexNew = // the previous index pointing to either
	  // (i) the first particle in this cell (if there was any)
	  // or (ii) the first particle in the next occupied cell
	  indexOld;

        // the combined offset for this particle
        indexNew += d_reductionOffsetsAdding[indexOld];

        // we have to consider the compaction offsets only if there were any active particles previously
        if (nparticles>0) {
          indexNew -= d_reductionOffsets[indexOld];

          // if the particle has previously been invalidated
          // its reduction offset is too high by one..
          // so we need to compensate for that
          int kold= get_key(d_keysOld[indexOld]);
          if (kold == 0xffffffff) indexNew++;
        } // nparticles > 0


        // The index points now to the position in the new list of this particle.
        // The space for the new particles was freed immediately *before* this particle,
        // so we need to substract the number of all these new particles.. We get this number by
        // finding the difference between the d_reductionOffsetsAdding for this particle and the
        // d_reductionOffsetsAdding for the previous one.. todo: this is a hack!!
        uint nnp =  d_reductionOffsetsAdding[indexOld];
        if (indexOld>0) nnp -= d_reductionOffsetsAdding[indexOld-1];

        // we substract this number from the new index and are now at the start of the newly added particles
        indexNew -= nnp;

        // we then need to add the domain index ..
        indexNew += d_domainReductionOffset[_cellIndex];

        // and finally add the running index for this particle
        indexNew += i;

        // phew. Now we still have to work out the "pool" index from where the new particles will
        // be obtained. The pool index is the end of the old list (nparticles) plus the running index for the
        // reduction offsets
        size_t indexPool = nparticles + d_reductionOffsetsAdding[indexOld] - nnp +  d_domainReductionOffset[_cellIndex] + i;

        // and write it
        //d_keysNew[indexNew] = KEY_INDEX_PAIR(666, i);
        id = get_index(d_keysOld[indexPool]);
        d_keysNew[indexNew] = KEY_INDEX_PAIR(_cellIndex, id);
      } // indexOld does point to a valid particle
      else {
        indexNew = indexOld;

        if (nparticles>0) {
          indexNew += d_reductionOffsetsAdding[indexOld-1];
          indexNew -= d_reductionOffsets[indexOld-1];
        } else {
          indexNew = 0;
        }

        indexNew += d_domainReductionOffset[_cellIndex];
        indexNew += i;

        size_t indexPool = nparticles;
        if (nparticles>0) indexPool += d_reductionOffsetsAdding[indexOld-1];
        indexPool +=  d_domainReductionOffset[_cellIndex] + i;

        id = get_index(d_keysOld[indexPool]);
        d_keysNew[indexNew] = KEY_INDEX_PAIR(_cellIndex, id);
      } // index old is last particle

      // finally we need to update the properties of the new particle with the
      // user-defined routine. We do that in the same manner as in the 
      // drift/diffusivity method
      /*_species_t _actualSpecies = {0,
                                  d_properties[id + 0 * propertiesOffset],
                                  d_properties[id + 1 * propertiesOffset],
                                  d_properties[id + 2 * propertiesOffset],
                                  d_properties[id + 3 * propertiesOffset]};*/

        _species_t _actualSpecies = {
            d_state[get_global_id_for_species_2d(speciesIndex)],
            d_state+get_global_area_2d*speciesIndex,
            0,
            d_properties[id + 0 * propertiesOffset],
            d_properties[id + 1 * propertiesOffset],
            d_properties[id + 2 * propertiesOffset],
            d_properties[id + 3 * propertiesOffset]
        };

      // huh .. use ternary operator paired with the comma operator
      // you never heard of the comma operator? neither have i ..
      /*
	This is what we want to achieve in the macro:
	if (! (ri%2)) {
	c.v[0]++;
	u.c = threefry2x32(c, k);
	rand = (u.r.x*M_RAN_INVM32+0.5); 
	} else {
	rand = (u.r.y*M_RAN_INVM32+0.5);
	}
	ri++;
      */
#define RANDOM (!(ri%2)) ? (ri++, c.v[0]++, u.c = threefry2x32(c, k), (u.r.x*M_RAN_INVM32+0.5)) : (ri++, u.r.y*M_RAN_INVM32+0.5)

      // <<<! newIndividualProperties !>>>

#undef RANDOM

      //and write it back
      d_properties[id + 0 * propertiesOffset] = _actualSpecies.DiffusivityX;
      d_properties[id + 1 * propertiesOffset] = _actualSpecies.DiffusivityY;
      d_properties[id + 2 * propertiesOffset] = _actualSpecies.DriftX;
      d_properties[id + 3 * propertiesOffset] = _actualSpecies.DriftY;

    } // loop over particles
}

__kernel void clearReductionOffsetsAdding(__global int *d_reductionOffsetsAdding) {
  d_reductionOffsetsAdding[get_global_id_2d] = 0;
}

__kernel void computeReductionOffsets(__global uint *d_newParticles,
				      __global uint *d_domainReductionOffset,
				      __global uchar *d_flags,
				      __global uint *d_cellIndicesLeft,
				      __global int *d_reductionOffsetsAdding
)
{
  // for each *flagged* cell, it sums the number of new particles in this cell
  // and all cells between this one and the previous flagged one (from segmented scan)
  // and stores the sum in d_reductionOffsetsAdding. d_reductionOffsetsAdding is a particle array
  // that holds for each individual the number of particles that should be added before its own.
  // this routine therefore stores the number of particles to be added at the start of the list of
  // particles in its own cell.
  // only do something if my cell is flagged
  const size_t _cellIndex = get_global_id_2d;

  // Now we need to insert particles before the first particle in this cell if
  // (1) the flag in the *next* cell is set (remember, the counter is reset in the next cell)
  // (2) we're at the end of the domain
  if (_cellIndex < get_global_area_2d-1) {
    if (d_flags[_cellIndex+1] == 1) {
      // we need to insert all particles created in previous cells and the ones
      // created in this cell
      d_reductionOffsetsAdding[d_cellIndicesLeft[_cellIndex]] =
         d_newParticles[_cellIndex] + d_domainReductionOffset[_cellIndex];
    }
  }

  if (_cellIndex == get_global_area_2d-1) {
    // if we're at the end of the domain we need to insert particles here as well.
    // if this cell is empty to start with, cellIndicesLeft will point to the first unused particle..
    // if it's not empty it'll point to the particle in this cell. We just write the index to either..
    d_reductionOffsetsAdding[d_cellIndicesLeft[_cellIndex]] =
       d_newParticles[_cellIndex] + d_domainReductionOffset[_cellIndex];
 }
} // compute combined reduction offset

// 1D helper functions.. sort of pointless :)
inline size_t get_global_thread_id_1d() {return get_global_id(0);}
inline size_t get_local_thread_id_1d() {return get_local_id(0);}
inline size_t get_local_thread_area_1d() {return get_local_size(0);}
inline size_t get_global_thread_area_1d() {return get_global_size(0);}
inline size_t get_block_thread_id_1d() {return get_group_id(0);}


__kernel void scanIndividuals(__global int *d_reductionOffsets,
			      __global int *d_reductionOffsetsBlocks,
			      __global int *d_reductionOffsetsAdding,
			      __global int *d_reductionOffsetsBlocksAdding,
			      __local int *s_local,
			      __local int *s_localAdding)
{
  // get local and global ids
  const uint lid = get_local_thread_id_1d();
  const uint tid = get_global_thread_id_1d();

  // store distance information in local
  s_local[2*lid]   = d_reductionOffsets[2*tid];
  s_local[2*lid+1] = d_reductionOffsets[2*tid+1];
  s_localAdding[2*lid]   = d_reductionOffsetsAdding[2*tid];
  s_localAdding[2*lid+1] = d_reductionOffsetsAdding[2*tid+1];

 // n is the number of local threads
  // we reduce over the padding threads as well.. this is ok
  // since the scatter function will take care of it
  // also the padding threads will always be in the last block
  // and we hence don't corrupt the partial sum for the lower blocks
  const uint n = 2*get_local_thread_area_1d();

  // add up all voids for this block - need another kernel to reduce
  // over all blocks
  // Up-sweep phase (see Harris et al. in GPU Gems)
  int offset = 1;
    
  for (size_t d = n>>1; d > 0; d >>= 1)
    {
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d) {
	int ai = offset*(2*lid+1)-1;
	int bi = offset*(2*lid+2)-1;
	s_local[bi] += s_local[ai];
	s_localAdding[bi] += s_localAdding[ai];
      }
      offset *= 2;
    }
    
  if (lid == 0) { s_local[n - 1] = 0; s_localAdding[n - 1] = 0;  } // clear the last element     

  for (int d = 1; d < n; d *= 2) // traverse down tree & build scan  
    {  
      offset >>= 1;  
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d)                       
	{  
	  int ai = offset*(2*lid+1)-1;  
	  int bi = offset*(2*lid+2)-1;  
	  float t = s_local[ai];  
	  float tAdding = s_localAdding[ai];  
	  s_local[ai] = s_local[bi];  
	  s_local[bi] += t;   
	  s_localAdding[ai] = s_localAdding[bi];  
	  s_localAdding[bi] += tAdding;   
	}  
    }  
  barrier(CLK_LOCAL_MEM_FENCE);

  // write results to device memory    
  // we want an *inclusive* so need to shift to the left
  d_reductionOffsets[2*tid] = s_local[2*lid+1];
  d_reductionOffsetsAdding[2*tid] = s_localAdding[2*lid+1];
  if (lid == get_local_thread_area_1d()-1) {
    // need to add last input array
    // TODO: this is ugly.. try to think of sth better!
    // better to do exclusive scan and finally add the input vector..
    // see: Sengupta, Harris, Zhang et Owens (2007): Scan Primitives for GPU Computing
    d_reductionOffsets[2*tid+1] = s_local[2*lid+1] + d_reductionOffsets[2*tid+1];
    d_reductionOffsetsBlocks[get_block_thread_id_1d()] = d_reductionOffsets[2*tid+1];
    d_reductionOffsetsAdding[2*tid+1] = s_localAdding[2*lid+1] + d_reductionOffsetsAdding[2*tid+1];
    d_reductionOffsetsBlocksAdding[get_block_thread_id_1d()] = d_reductionOffsetsAdding[2*tid+1];
  } else {
    d_reductionOffsets[2*tid+1] = s_local[2*lid+2];  
    d_reductionOffsetsAdding[2*tid+1] = s_localAdding[2*lid+2];  
  }

  
  /* // THIS WOULD BE EXCLUSIVE SCAN
  d_reductionOffsets[2*tid] = s_local[2*lid];
  d_reductionOffsets[2*tid+1] = s_local[2*lid+1];

  if (lid == get_local_thread_area_1d()-1)
  d_reductionOffsetsBlocks[get_block_thread_id_1d()] = s_local[2*lid+1];*/

} // scan over individuals

// reduce over all blocks
// todo : i don't think we need the original offsets here ..
__kernel void scanIndividualsBlocks(__global int *d_reductionOffsets, // remove me
				    __global int *d_reductionOffsetsBlocks,
				    __global int *d_reductionOffsetsAdding, // remove me
				    __global int *d_reductionOffsetsBlocksAdding,
				    __local int *s_local,
				    __local int *s_localAdding)
{
  // get local and global ids
  const uint lid = get_local_thread_id_1d();
  const uint tid = get_global_thread_id_1d();

  // store distance information in local
  s_local[2*lid]   = d_reductionOffsetsBlocks[2*tid];
  s_local[2*lid+1] = d_reductionOffsetsBlocks[2*tid+1];
  s_localAdding[2*lid]   = d_reductionOffsetsBlocksAdding[2*tid];
  s_localAdding[2*lid+1] = d_reductionOffsetsBlocksAdding[2*tid+1];
  
  // n is the number of local threads
  uint n = 2*get_local_thread_area_1d();

  // add up all voids for this block - need another kernel to reduce
  // over all blocks
  // Up-sweep phase (see Harris et al. in GPU Gems)
  int offset = 1;
    
  for (size_t d = n>>1; d > 0; d >>= 1)
    {
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d) {
	int ai = offset*(2*lid+1)-1;
	int bi = offset*(2*lid+2)-1;
	s_local[bi] += s_local[ai];
	s_localAdding[bi] += s_localAdding[ai];
      }
      offset *= 2;
    }
    
  if (lid == 0) { s_local[n - 1] = 0; s_localAdding[n - 1] = 0;  } // clear the last element     

  
  for (int d = 1; d < n; d *= 2) // traverse down tree & build scan  
    {  
      offset >>= 1;  
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d)                       
	{  
	  int ai = offset*(2*lid+1)-1;  
	  int bi = offset*(2*lid+2)-1;  
	  float t = s_local[ai];  
	  s_local[ai] = s_local[bi];  
	  s_local[bi] += t;   
	  float tAdding = s_localAdding[ai];  
	  s_localAdding[ai] = s_localAdding[bi];  
	  s_localAdding[bi] += tAdding;   
	}  
    }  
  barrier(CLK_LOCAL_MEM_FENCE);  

  // write results to device memory  
  // we need *exclusive* scan
  d_reductionOffsetsBlocks[2*tid] = s_local[2*lid];
  d_reductionOffsetsBlocks[2*tid+1] = s_local[2*lid+1];  
  d_reductionOffsetsBlocksAdding[2*tid] = s_localAdding[2*lid];
  d_reductionOffsetsBlocksAdding[2*tid+1] = s_localAdding[2*lid+1];  
} // scan individuals blocks

__kernel void scatterIndividuals(__global type_key *d_keys,
				 __global type_key *d_keysOld,
				 __global int *d_reduceOffsets,
				 __global int *d_reduceOffsetsBlocks,
				 __global int *d_reductionOffsetsAdding,
				 __global int *d_reductionOffsetsBlocksAdding,
				 __global int *d_totalAdded,
				 __global int *d_total,
				 __global type_key *d_lostParticles,
                                 uint nparticles,
                                 __global uint *d_domainReductionOffset,
                                 __global uint *d_newParticles)
{
  // get local and global ids
  const uint lid = get_local_thread_id_1d();
  const uint tid = get_global_thread_id_1d();

  // we need to copy two so that we can use the same block partition as scan
  type_key key1 = d_keysOld[2*tid];
  type_key key2 = d_keysOld[2*tid+1];

  // check if we're still in the particles range
  if (2*tid < nparticles) {
    uint relIndexl = d_reduceOffsets[2*tid] + d_reduceOffsetsBlocks[get_block_thread_id_1d()];
    uint relIndexa = d_reductionOffsetsAdding[2*tid] + d_reductionOffsetsBlocksAdding[get_block_thread_id_1d()];
    d_reduceOffsets[2*tid] = relIndexl;
    d_reductionOffsetsAdding[2*tid] = relIndexa;

    // write only if the state is still alive
    if (get_key(key1) < UINT_MAX) {
      // and store it
      d_keys[2*tid - relIndexl + relIndexa] = key1;
    }  else {
      // write to lost particles area
      d_lostParticles[relIndexl-1] = key1;
    }
  } 

  if (2*tid+1 < nparticles) {
    uint relIndexl = d_reduceOffsets[2*tid+1] + d_reduceOffsetsBlocks[get_block_thread_id_1d()];
    uint relIndexa = d_reductionOffsetsAdding[2*tid+1] + d_reductionOffsetsBlocksAdding[get_block_thread_id_1d()];
    // and write the combined indices for later
    d_reduceOffsets[2*tid+1] = relIndexl;
    d_reductionOffsetsAdding[2*tid+1] = relIndexa;

    // write only if the state is still alive
    if (get_key(key2) < UINT_MAX) {
      // and store it
      d_keys[2*tid+1 - relIndexl + relIndexa] = key2;
    } else {
      // write to lost particles area
      d_lostParticles[relIndexl-1] = key2;
    }
  } 

  // last particles thread writes total kill number to device
  if (2*tid == nparticles-1 || 2*tid +1 == nparticles-1 ) {    
    // index for the last particle
    /*
    uint relIndexl = d_reduceOffsets[nparticles-1] + d_reduceOffsetsBlocks[get_block_thread_id_1d()];
    uint relIndexa = d_reductionOffsetsAdding[nparticles-1] + d_reductionOffsetsBlocksAdding[get_block_thread_id_1d()];*/
    // block index has been added before
    uint relIndexl = d_reduceOffsets[nparticles-1];
    uint relIndexa = d_reductionOffsetsAdding[nparticles-1];
    *d_total = nparticles-relIndexl;
    *d_totalAdded = relIndexa;

    // now we need to check if we're in the last cell .. cause if not we'll need to
    // add all particles that were created after me!
    type_key mkey;
    if (2*tid   == nparticles-1) mkey = d_keysOld[2*tid];
    if (2*tid+1 == nparticles-1) mkey = d_keysOld[2*tid+1];

    if (get_key(mkey)<get_global_area_2d-1) {
        // add all intermediate particles
        (*d_totalAdded) = (*d_totalAdded)
                            + d_domainReductionOffset[get_global_area_2d-1]
                            + d_newParticles[get_global_area_2d-1];
    }
  }
}

// todo: we can probably merge this kernel with the padding kernel..
__kernel void writeLostParticles(__global type_key *d_keys,
				 __global type_key *d_lostParticles,
				 uint nparticles)
{
  // get global id
  const uint tid = get_global_thread_id_1d();

  // and write it
  d_keys[nparticles+tid] = d_lostParticles[tid];
}

// todo: surely this kernel fits somewhere else ..
__kernel void copyTempToKeys(__global type_key *d_tempKeys,
			     __global type_key *d_keys) {
  d_keys[get_global_id(0)] = d_tempKeys[get_global_id(0)];
}

// we need to go 2D for the segmented scan .. since the second-level scan needs to work on all blocks in one block
// which means we can't have more than 1024x1024 blocks
// We need to define our own.. the GPGMP functions only work with the actual grid!!
// We need to make sure that consecutive threads access consecutive memory locations
inline size_t get_global_thread_area_2d() {return get_global_size(0) * get_global_size(1);}
inline size_t get_local_thread_id_2d() {return get_local_id(0) + (get_local_id(1) * get_local_size(0));}
inline size_t get_local_thread_area_2d() {return get_local_size(0) * get_local_size(1);}
inline size_t get_block_thread_id_2d() {return (get_group_id(0)  + (get_group_id(1)  * get_num_groups(0)));}
inline size_t get_global_thread_id_2d() {return get_block_thread_id_2d()*get_local_thread_area_2d() + get_local_thread_id_2d();}

// this performs a complete segmented scan in one block (second-level scan)
// implemented after
// see: Sengupta, Harris, Zhang et Owens (2007): Scan Primitives for GPU Computing
// (SHZO07)
__kernel void segmentedScanBlock(__global uint *d_data,
				 __global uchar *d_flags,
				 __global uchar *d_inputFlags,
				 __local uint *s_data,
				 __local uchar *s_flags,
				 __local uchar *s_inputFlags)
{
  // get local and global ids
  const uint lid = get_local_thread_id_2d();
  const uint tid = get_global_thread_id_2d();

  // store data information in local
  s_data[2*lid]   = d_data[2*tid];
  s_data[2*lid+1] = d_data[2*tid+1];

  // store flag information in local
  s_flags[2*lid]   = d_flags[2*tid];
  s_flags[2*lid+1] = d_flags[2*tid+1];
  s_inputFlags[2*lid]   = d_inputFlags[2*tid];
  s_inputFlags[2*lid+1] = d_inputFlags[2*tid+1];
  
  // scan over whole block (padding threads will be disregarded later)
  const uint n = 2*get_local_thread_area_2d();

  // Up-sweep phase
  int offset = 1; // offset = 2^d (in SHZO07)
  for (size_t d = n>>1; d > 0; d >>= 1)
    // d is the number of threads we need for each sweep
    // starting at n/2
    {
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d) {
	int ai = offset*(2*lid+1)-1; // corresponds to k+2^d-1 in SHZO07
	int bi = offset*(2*lid+2)-1; // corresponds to k+2^(d+1)-1 in SHZO07
	if (s_flags[bi]==0) // not set 
	  s_data[bi] += s_data[ai];
	s_flags[bi] = s_flags[ai] | s_flags[bi];
      }
      offset *= 2;
    }
  
  if (lid == 0) { s_data[n - 1] = 0; } // clear last element     

  for (int d = 1; d < n; d *= 2) // traverse down tree & build scan  
    {  
      offset >>= 1;  
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d)                       
	{  
	  int ai = offset*(2*lid+1)-1;  
	  int bi = offset*(2*lid+2)-1;  
	  float t = s_data[ai];  
	  s_data[ai] = s_data[bi];  
	  
	  if (s_inputFlags[ai+1]==1) // if fi[k+2^d] set (this is the *input* flag vector)
	    s_data[bi] = 0; // x[k+2^(d+1)-1] <- 0
	  else if (s_flags[ai]==1) // if f[k+2^d-1] set
	    s_data[bi] = t; // x[k+2^(d+1)-1] <- t
	  else
	    s_data[bi] += t; // x[k+2^(d+1)-1] <- t + x[k+2^(d+1)-1]
	  
	  s_flags[ai]=0; // unset f[k+2^d-1]
	}  
    }  
  barrier(CLK_LOCAL_MEM_FENCE);  

  // write back to device
  d_data[2*tid] = s_data[2*lid];
  d_data[2*tid+1] = s_data[2*lid+1];
} // segmented scan block

// this performs an upsweep and stores the workgroup results
__kernel void segmentedScanUpSweep(__global uint *d_data,
				   __global uchar *d_flags,
				   __global uchar *d_partialFlags,
				   __global uint *d_dataBlocks,
				   __global uchar *d_flagsBlocks,
				   __global uchar *d_firstBlockFlag,
				   __local uint *s_data,
				   __local uchar *s_flags)
{
  // get local and global ids
  const uint lid = get_local_thread_id_2d();
  const uint tid = get_global_thread_id_2d();

  // store data information in local
  s_data[2*lid]   = d_data[2*tid];
  s_data[2*lid+1] = d_data[2*tid+1];

  // store flag information in local
  s_flags[2*lid]   = d_flags[2*tid];
  s_flags[2*lid+1] = d_flags[2*tid+1];

  // scan over whole block (padding threads will be disregarded later)
  const uint n = 2*get_local_thread_area_2d();

  // Up-sweep phase
  int offset = 1; // offset = 2^d (in SHZO07)
  for (size_t d = n>>1; d > 0; d >>= 1)
    // d is the number of threads we need for each sweep
    // starting at n/2
    {
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d) {
	int ai = offset*(2*lid+1)-1; // corresponds to k+2^d-1 in SHZO07
	int bi = offset*(2*lid+2)-1; // corresponds to k+2^(d+1)-1 in SHZO07
	if (s_flags[bi]==0) // not set 
	  s_data[bi] += s_data[ai];
	s_flags[bi] = s_flags[ai] | s_flags[bi];
      }
      offset *= 2;
    }

  // offset = n = 2*block size now ..

  // and store results
  barrier(CLK_LOCAL_MEM_FENCE);  

  // store partial sum and OR tree
  d_data[2*tid] = s_data[2*lid];
  d_data[2*tid+1] = s_data[2*lid+1];  

  d_partialFlags[2*tid] = s_flags[2*lid];
  d_partialFlags[2*tid+1] = s_flags[2*lid+1];  

  // store block values
  if (lid==0) {
    // last value of data and OR
    d_dataBlocks[get_block_thread_id_2d()] = s_data[2*(get_local_thread_area_2d()-1)+1];
    d_flagsBlocks[get_block_thread_id_2d()] = s_flags[2*(get_local_thread_area_2d()-1)+1];

    // first value of partial OR flag vector
    d_firstBlockFlag[get_block_thread_id_2d()] = s_flags[0];

  }
}

// this performs an upsweep and stores the workgroup results
__kernel void segmentedScanDownSweep(__global uint *d_data,
				     __global uchar *d_flags,
				     __global uchar *d_partialFlags,
				     __global uint *d_dataBlocks,
				     __global uchar *d_flagsBlocks,
				     __local uint *s_data,
				     __local uchar *s_flags,
				     __local uchar *s_inputFlags
)
{
  // get local and global ids
  const uint lid = get_local_thread_id_2d();
  const uint tid = get_global_thread_id_2d();

  // store data information in local - this is the partial sum
  s_data[2*lid]   = d_data[2*tid];
  s_data[2*lid+1] = d_data[2*tid+1];
  
  // store partial flag information in local
  s_flags[2*lid]   = d_partialFlags[2*tid];
  s_flags[2*lid+1] = d_partialFlags[2*tid+1];
  s_inputFlags[2*lid]   = d_flags[2*tid];
  s_inputFlags[2*lid+1] = d_flags[2*tid+1];

  // Down-sweep phase
  const uint n = 2*get_local_thread_area_2d();
  int offset = n; // offset = 2^d (in SHZO07)

  barrier(CLK_LOCAL_MEM_FENCE); // this block fence is important .. otherwise s_data[n-1] might be overwritten again
  if (lid == 0) { s_data[n - 1] = d_dataBlocks[get_block_thread_id_2d()]; } // set it to block scan valua

  for (int d = 1; d < n; d *= 2) // traverse down tree & build scan  
    {  
      offset >>= 1;  
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d)                       
	{  
	  int ai = offset*(2*lid+1)-1;  
	  int bi = offset*(2*lid+2)-1;  
	  float t = s_data[ai];  
	  s_data[ai] = s_data[bi];  
	  
	  if (s_inputFlags[ai+1]==1) // if fi[k+2^d] set (this is the *input* flag vector)
	    s_data[bi] = 0; // x[k+2^(d+1)-1] <- 0
	  else if (s_flags[ai]==1) // if f[k+2^d-1] set
	    s_data[bi] = t; // x[k+2^(d+1)-1] <- t
	  else
	    s_data[bi] += t; // x[k+2^(d+1)-1] <- t + x[k+2^(d+1)-1]
	  
	  s_flags[ai]=0; // unset f[k+2^d-1]
	}  
    }  
  barrier(CLK_LOCAL_MEM_FENCE);    

  // and finally save it
  d_data[2*tid] = s_data[2*lid];
  d_data[2*tid+1] = s_data[2*lid+1];

}

__kernel void segmentedScanBlock_single(__global uint *d_data,
					__global uchar *d_flags,
					__global uchar *d_inputFlags,
					__local uint *s_data,
					__local uchar *s_flags,
					__local uchar *s_inputFlags)
{
  // get local and global ids
  const uint lid = get_local_thread_id_2d();
  const uint tid = get_global_thread_id_2d();

  // store data information in local
  s_data[2*lid]   = d_data[2*tid];
  s_data[2*lid+1] = d_data[2*tid+1];

  // store flag information in local
  s_flags[2*lid]   = d_flags[2*tid];
  s_flags[2*lid+1] = d_flags[2*tid+1];
  s_inputFlags[2*lid]   = d_inputFlags[2*tid];
  s_inputFlags[2*lid+1] = d_inputFlags[2*tid+1];
  
  // scan over whole block (padding threads will be disregarded later)
  const uint n = 2*get_local_thread_area_2d();

  // Up-sweep phase
  int offset = 1; // offset = 2^d (in SHZO07)
  for (size_t d = n>>1; d > 0; d >>= 1)
    // d is the number of threads we need for each sweep
    // starting at n/2
    {
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d) {
	int ai = offset*(2*lid+1)-1; // corresponds to k+2^d-1 in SHZO07
	int bi = offset*(2*lid+2)-1; // corresponds to k+2^(d+1)-1 in SHZO07
	if (s_flags[bi]==0) // not set 
	  s_data[bi] += s_data[ai];
	s_flags[bi] = s_flags[ai] | s_flags[bi];
      }
      offset *= 2;
    }

  if (lid == 0) { s_data[n - 1] = 0; } // clear last element     

  for (int d = 1; d < n; d *= 2) // traverse down tree & build scan  
    {  
      offset >>= 1;  
      barrier(CLK_LOCAL_MEM_FENCE);
      if (lid < d)                       
	{  
	  int ai = offset*(2*lid+1)-1;  
	  int bi = offset*(2*lid+2)-1;  
	  float t = s_data[ai];  
	  s_data[ai] = s_data[bi];  
	  
	  if (s_inputFlags[ai+1]==1) // if fi[k+2^d] set (this is the *input* flag vector)
	    s_data[bi] = 0; // x[k+2^(d+1)-1] <- 0
	  else if (s_flags[ai]==1) // if f[k+2^d-1] set
	    s_data[bi] = t; // x[k+2^(d+1)-1] <- t
	  else
	    s_data[bi] += t; // x[k+2^(d+1)-1] <- t + x[k+2^(d+1)-1]
	  
	  s_flags[ai]=0; // unset f[k+2^d-1]
	}  
    }  
  barrier(CLK_LOCAL_MEM_FENCE);  

  // write back to device
  d_data[2*tid] = s_data[2*lid];
  d_data[2*tid+1] = s_data[2*lid+1];
} // segmented scan block
