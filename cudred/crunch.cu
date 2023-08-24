#include <cuda_runtime.h>

void calculateSumCombinations(float* input1, float* input2, float* output, int N)
{
    // Declare device memory pointers
    float *d_input1, *d_input2, *d_output;
    
    // Allocate memory on the GPU
    cudaMalloc((void**)&d_input1, N * sizeof(float));
    cudaMalloc((void**)&d_input2, N * sizeof(float));
    cudaMalloc((void**)&d_output, N * N * sizeof(float));
    
    // Copy input data from host to device
    cudaMemcpy(d_input1, input1, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_input2, input2, N * sizeof(float), cudaMemcpyHostToDevice);
    
    // Define the block and grid dimensions
    dim3 blockSize(16, 16);
    dim3 gridSize((N + blockSize.x - 1) / blockSize.x, (N + blockSize.y - 1) / blockSize.y);
    
    // Launch the kernel function
    sumCombinations<<<gridSize, blockSize>>>(d_input1, d_input2, d_output, N);
    
    // Wait for the GPU computation to complete
    cudaDeviceSynchronize();
    
    // Copy the result from device to host
    cudaMemcpy(output, d_output, N * N * sizeof(float), cudaMemcpyDeviceToHost);
    
    // Free the GPU memory
    cudaFree(d_input1);
    cudaFree(d_input2);
    cudaFree(d_output);
}