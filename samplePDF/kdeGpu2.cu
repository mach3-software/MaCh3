// -*- c++ -*-
#define BLOCK_SIZE 256
#define GRID_SIZE 1000
#define MULTI_DATA 10

//MaCh3 included
#include "manager/gpuUtils.cu"


#include "stdio.h"
#include <iostream>
#include <assert.h>
#include "math_constants.h"
#include <iomanip>

__global__ void calcKde2d(int n_data, int n_mc, 
                          double *d_data_x, double *d_data_y, 
                          double *d_mc_x,
                          double *d_mc_y,
                          double *d_mc_sx,
                          double *d_mc_sy,
                          double *d_mcweights,
                          //double *d_pdfatdata)//, 
                          double *d_datamcweights)
{
  __shared__ float s_mc_x[BLOCK_SIZE];
  __shared__ float s_mc_y[BLOCK_SIZE];
  __shared__ float s_mcweights[BLOCK_SIZE];
  __shared__ float s_mc_sx[BLOCK_SIZE];
  __shared__ float s_mc_sy[BLOCK_SIZE];
  
  float l_data_x[MULTI_DATA];
  float l_data_y[MULTI_DATA];
  float l_pdfatdata[MULTI_DATA];
  
  const int theIndex(blockIdx.x * blockDim.x + threadIdx.x);
  const int n_mcChunks(n_mc/blockDim.x+1);
  const int n_dataChunks(n_data/MULTI_DATA+1);
  const int dataChunk(theIndex%n_dataChunks);
  const int mcChunk(theIndex/n_dataChunks);

  int lo_data_idx(0);
  for (int dataIdx=dataChunk;dataIdx<n_data;dataIdx+=n_dataChunks)
    {
      l_data_x[lo_data_idx]=d_data_x[dataIdx];
      l_data_y[lo_data_idx]=d_data_y[dataIdx];
      l_pdfatdata[lo_data_idx]=0;
      lo_data_idx++;
    }  

  if(mcChunk<n_mcChunks)
    {
      const int firstMCIdx(mcChunk*blockDim.x);
      const int loadMCIdx(firstMCIdx+threadIdx.x);
      if(loadMCIdx<n_mc)
        {
          s_mc_x[threadIdx.x]=d_mc_x[loadMCIdx];
          s_mc_y[threadIdx.x]=d_mc_y[loadMCIdx];
          s_mc_sx[threadIdx.x]=d_mc_sx[loadMCIdx];
          s_mc_sy[threadIdx.x]=d_mc_sy[loadMCIdx];
          s_mcweights[threadIdx.x]=d_mcweights[loadMCIdx];
        }
      else
        {
          s_mc_x[threadIdx.x]=100;
          s_mc_y[threadIdx.x]=0;
          s_mc_sx[threadIdx.x]=100;
          s_mc_sy[threadIdx.x]=10;
          s_mcweights[threadIdx.x]=0.0;
        }
      
      __syncthreads();  
      
      for (int i=0;i<blockDim.x;i++)
        {
          const int l_MCIndex((i));//+threadIdx.x)%blockDim.x);
          const float l_sigmax(s_mc_sx[l_MCIndex]);
          const float l_sigmay(s_mc_sy[l_MCIndex]);
          const float l_minuspt5oversigmaxsqrd(-0.5/l_sigmax/l_sigmax);
          const float l_minuspt5oversigmaysqrd(-0.5/l_sigmay/l_sigmay);
          //const float l_sigmaxsqrd(l_sigmax*l_sigmax);
          //const float l_sigmaysqrd(l_sigmay*l_sigmay);
          //const float l_sigmaxy(l_sigmax*l_sigmay);
          const float l_mc_x(s_mc_x[l_MCIndex]);
          const float l_mc_y(s_mc_y[l_MCIndex]);
          const float l_mcweightOverSigmaxy(s_mcweights[l_MCIndex]/
                                            l_sigmax/l_sigmay);

          for(int l_data_idx=0; l_data_idx<MULTI_DATA; l_data_idx++)
            {
              float l_data_minus_mc_x(l_data_x[l_data_idx]-l_mc_x);
              float l_data_minus_mc_y(l_data_y[l_data_idx]-l_mc_y);
              
              l_pdfatdata[l_data_idx]+=
                
                __expf(
                       
                       l_minuspt5oversigmaxsqrd*
                       l_data_minus_mc_x*l_data_minus_mc_x+

                       l_minuspt5oversigmaysqrd*
                       l_data_minus_mc_y*l_data_minus_mc_y
                       )
                *l_mcweightOverSigmaxy;
            }
        }

      int l_data_idx=0;

      for (int dataIdx=dataChunk;dataIdx<n_data;dataIdx+=n_dataChunks)
        {
          int outIdx(dataIdx+n_data*mcChunk);
          d_datamcweights[outIdx]=l_pdfatdata[l_data_idx];
          l_data_idx++;
        }
    }
  return;
}


__global__ void calcPdfAtData2d(int n_data, int n_things,  
                                double *d_datamcweights, 
                                double *d_pdfatdata)
{
  int thisIdx = (blockIdx.x * blockDim.x + threadIdx.x);
  if(thisIdx<n_data)
    {
      double thisValue(0);
      for(int thingIdx=0; thingIdx<n_things; thingIdx++) 
        {
          int bigIdx(thisIdx+n_data*thingIdx);
          thisValue+=d_datamcweights[bigIdx];
        }
      d_pdfatdata[thisIdx]=thisValue/2.0/CUDART_PI;
    }
  return;
}


extern "C" __host__ double* calcKdeWeights2d(int n_data, int n_mc,
                                             double sigmax, 
                                             double sigmay, 
                                             double h_data_x[], 
                                             double h_data_y[], 
                                             double h_mc_x[],
                                             double h_mc_y[],
                                             double h_mc_sx[],
                                             double h_mc_sy[],
                                             double h_mcweights[])
{
  dim3 block_size;
  block_size.x = BLOCK_SIZE;
  dim3 grid_size;

  grid_size.x = GRID_SIZE;

  int n_mcChunks(n_mc/block_size.x+1);
  int n_mcData(n_data*n_mcChunks);
  std::cout << "n_mc " << n_mc << "     n_data " << n_data 
            << "           n_mcData " << n_mcData << std::endl;

  double *d_mc_x, *d_mc_y,
    *d_mc_sx, *d_mc_sy,
    *d_data_x, *d_data_y,
    *d_mcweights, *d_datamcweights, 
    *d_pdfatdata;

  size_t datasize   = n_data   * sizeof(double);
  size_t mcsize     = n_mc     * sizeof(double);
  size_t mcDatasize = n_mcData * sizeof(double);

  assert(cudaSuccess==cudaMalloc(&d_pdfatdata, datasize));
  double* h_pdfatdata = (double *)malloc(datasize);//n_data);


  assert(cudaSuccess==cudaMalloc(&d_datamcweights, mcDatasize)); 
  //Useful for testing but very slow
  //double* h_datamcweights = (double *)malloc(mcDatasize);
  //mcDatasize);//n_data);
  //for(int i=0;i<n_mcData;i++) h_datamcweights[i]=100000.0;
  //cudaMemcpy(d_datamcweights, h_datamcweights, mcDatasize, cudaMemcpyHostToDevice);
  
  assert(cudaSuccess==cudaMalloc(&d_data_x, datasize));
  cudaMemcpy(d_data_x, h_data_x, datasize, cudaMemcpyHostToDevice);
  assert(cudaSuccess==cudaMalloc(&d_data_y, datasize));
  cudaMemcpy(d_data_y, h_data_y, datasize, cudaMemcpyHostToDevice);
  
  assert(cudaSuccess==cudaMalloc(&d_mcweights, mcsize));
  cudaMemcpy(d_mcweights, h_mcweights, mcsize, cudaMemcpyHostToDevice);  
  assert(cudaSuccess==cudaMalloc(&d_mc_x, mcsize));
  cudaMemcpy(d_mc_x, h_mc_x, mcsize, cudaMemcpyHostToDevice);
  assert(cudaSuccess==cudaMalloc(&d_mc_y, mcsize));
  cudaMemcpy(d_mc_y, h_mc_y, mcsize, cudaMemcpyHostToDevice);
  assert(cudaSuccess==cudaMalloc(&d_mc_sx, mcsize));
  cudaMemcpy(d_mc_sx, h_mc_sx, mcsize, cudaMemcpyHostToDevice);
  assert(cudaSuccess==cudaMalloc(&d_mc_sy, mcsize));
  cudaMemcpy(d_mc_sy, h_mc_sy, mcsize, cudaMemcpyHostToDevice);


  const int n_dataChunks(n_data/MULTI_DATA+1);
  grid_size.x=(n_dataChunks/block_size.x+1)*n_mcChunks;

  std::cout << "Grid size " << grid_size.x << "    block size " 
            << block_size.x << std::endl;

  calcKde2d<<<grid_size, block_size>>>(n_data, n_mc, 
                                       d_data_x, d_data_y,
                                       d_mc_x, d_mc_y,
                                       d_mc_sx, d_mc_sy,
                                       d_mcweights,
                                       d_datamcweights);
  CudaCheckError();

  grid_size.x=n_data/block_size.x+1;
  calcPdfAtData2d<<<grid_size, block_size>>>(n_data, n_mcChunks,
                                             d_datamcweights, d_pdfatdata);
  CudaCheckError();
  
  cudaMemcpy(h_pdfatdata, d_pdfatdata, datasize, cudaMemcpyDeviceToHost);

  //For testing
  /*
  cudaMemcpy(h_datamcweights, d_datamcweights, mcDatasize, 
        cudaMemcpyDeviceToHost);
  int ooooo(0);
  int ooooo2(0);
  for(int i=0;i<n_mcData;i++) 
    {
      //if(h_datamcweights[i]>0&&i/n_data>=1.0)
      //{
      //  std::cout << std::setprecision(23)  
      //            << "ooooo " << h_datamcweights[i] 
      //            << "   " << int(i/n_data) << std::endl;
      //}
      if(h_datamcweights[i]<0) 
        {
          ooooo2++;
        }

      if(h_datamcweights[i]==0) 
        {
          //std::cout << "ooooo " << i << std::endl;
          ooooo++;
        }
    }
  std::cout << "ooooo n " << ooooo <<"/"<< n_mcData << std::endl;
  std::cout << "ooooo2 n " << ooooo2 <<"/"<< n_mcData << std::endl;
  */
  
  cudaFree(d_mc_x);
  cudaFree(d_mc_y);
  cudaFree(d_mc_sx);
  cudaFree(d_mc_sy);
  cudaFree(d_data_x);
  cudaFree(d_data_y);
  cudaFree(d_mcweights);
  cudaFree(d_datamcweights);
  cudaFree(d_pdfatdata);

  //free(h_datamcweights);

  return h_pdfatdata;
}

