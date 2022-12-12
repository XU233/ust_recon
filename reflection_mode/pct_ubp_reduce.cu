/**
 * @file pct_ubp.cu
 * 
 * CUDA code to perform ubp on GPU
 */



/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 */
__global__ void ubp( 
                      float * out,
                      const float * x, 
                      const float * y, 
                      const float * z,
                      const float * xk, 
                      const float * yk, 
                      const float * zk,
					  const float * xe,
					  const float * ye,
                      const unsigned int Nx,
                      const unsigned int log2Nx,
                      const unsigned int Nphi, 
                      const unsigned int Nt,
                      const float fixedDelay, 
                      const float tv, 
                      const float * data ) {
    // Work out which thread we are

    __shared__ float sample[512];
            
    int xIdx = blockIdx.x &(Nx-1);
    int angleIdx = blockIdx.x >> log2Nx;
    int transducerIdx = angleIdx * blockDim.x + threadIdx.x; 
    int globalVolumeIndex = xIdx + blockIdx.y * Nx + blockIdx.z * Nx * gridDim.y;
    
    
    // Get our X and Y coords
    float const xi = x[xIdx];
    float const yi = y[blockIdx.y];
    float const zi = z[blockIdx.z];


    //unsigned long start = 0;

    float const xki = xk[transducerIdx];
    float const yki = yk[transducerIdx];
    float const zki = zk[transducerIdx];
	
	float const xei = xe[angleIdx];
    float const yei = ye[angleIdx];
    
	float const d = sqrt( (xi-xki) * (xi-xki) + (yi-yki) * (yi-yki) + (zi-zki) * (zi-zki));
	float const d_emit = sqrt( (xi-xei) * (xi-xei) + (yi-yei) * (yi-yei) );
    float idxf = rintf( (d + d_emit)* tv) - fixedDelay;
    int idx = __float2int_rd(idxf);
    idx = max( idx, 0);
    idx = min( idx, Nt - 1); 
    idx = threadIdx.x + blockDim.x*idx + blockDim.x*Nt*angleIdx;
    //idx = 0;
    sample[threadIdx.x] = data[idx];
    __syncthreads();
	int nTotalThreads = blockDim.x;	// Total number of active threads

	while(nTotalThreads > 1)
	{
		int halfPoint = (nTotalThreads >> 1);	// divide by two
		// only the first half of the threads will be active.
		if (threadIdx.x < halfPoint)
		{
			
			// when calculating the average, sum and divide
			sample[threadIdx.x] += sample[threadIdx.x + halfPoint];
		}
		__syncthreads();

		nTotalThreads = (nTotalThreads >> 1);	// divide by two.
	}

	// At this point in time, thread zero has the min, max, and average
	// It's time for thread zero to write it's final results.
	// Note that the address structure of pResults is different, because
	// there is only one value for every thread block.

	if (threadIdx.x == 0)
	{
        atomicAdd(&out[globalVolumeIndex], sample[0]);
	}
    //atomicAdd(&out[globalBlockIndex],data[idx]);

            
}
