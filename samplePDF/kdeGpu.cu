//MaCh3 included
#include "manager/gpuUtils.cu"

// -*- c++ -*-
#include "stdio.h"
#include <iostream>
#include <assert.h>
#include "math_constants.h"

__constant__ double c_data_x[1000];

__global__ void calcKde(int n_data, int n_mc, double sigma, //double *d_data_x, 
                        double *d_mc_x,// double *d_sigmas, 
                        double *d_mcweights,
                        double *d_pdfatdata,
                        double *d_datamcweights)
{
  int theIndex(blockIdx.x * blockDim.x + threadIdx.x);
  int wholeThing(blockDim.x * gridDim.x);
  int thingsForEachData(wholeThing/n_data);
  int dataIdx(theIndex/thingsForEachData);
  if(dataIdx<n_data)
    {
      const double l_data_x(c_data_x[dataIdx]);
      const double  sigmasqrd(sigma*sigma);
      double l_sum(0.0);
      int firstMcIdx(theIndex-dataIdx*thingsForEachData);

      for (int mcIdx=firstMcIdx;mcIdx<n_mc;mcIdx+=thingsForEachData)
        {
          l_sum+= d_mcweights[mcIdx]*exp(-0.5*(l_data_x-d_mc_x[mcIdx])*
                                         (l_data_x-d_mc_x[mcIdx])/sigmasqrd);
        }
      l_sum/=sqrt(2.0*CUDART_PI);//This might be better off elsewhere
      l_sum/=sigma;//but it probably doesn't matter
      d_datamcweights[dataIdx+n_data*firstMcIdx]=l_sum;//shared here?
      //or swap order so coalesced and use shared memory in calcpdfdata
    }
  return;
}

__global__ void calcPdfAtData(int n_data, int n_mc,  double sigma, 
                              double *d_datamcweights, 
                              double *d_pdfatdata)
{
  int thisIdx = (blockIdx.x * blockDim.x + threadIdx.x);
  if (thisIdx < n_data)
    {
      double thisValue(0);
      for(int mcIdx=0; mcIdx<n_mc; mcIdx++) 
        {
          int bigIdx(thisIdx+n_data*mcIdx);
          thisValue+=d_datamcweights[bigIdx];
        }
      d_pdfatdata[thisIdx]=thisValue;
    }
  return;
}


extern "C" __host__ double* calcKdeWeights(int n_data, int n_mc,
                                           double sigma, 
                                           //double norm,
                                           double h_data_x[], 
                                           double h_mc_x[],
                                           //double h_sigmas[],
                                           double h_mcweights[])
{
  if(n_data>1000)
    std::cout << "Bad news: a maximum limit of 1000 data events for any given sample is hardcoded into kdeGpu.cu. This is the price you pay for using lovely fast constant memory."
              << std::endl << "Good news: just change c_data_x[1000] to c_data_x[whatever]. You have 64KB to play with." << std::endl;
  assert(n_data<=1000);

  dim3 block_size;
  block_size.x = 256;//1024;
  dim3 grid_size;
  grid_size.x = 100;
  assert(block_size.x*grid_size.x>(unsigned int)n_data);

  int wholeThing(block_size.x * grid_size.x);
  int thingsForEachData(wholeThing/n_data);
  int n_mcData(n_data*thingsForEachData);//n_mc);

  double *d_mc_x, //*d_sigmas,
    *d_mcweights, *d_datamcweights, 
    *d_pdfatdata;
  
  size_t datasize   = n_data   * sizeof(double);
  size_t mcsize     = n_mc     * sizeof(double);
  size_t mcDatasize = n_mcData * sizeof(double);

  cudaMalloc(&d_pdfatdata, datasize);
  double* h_pdfatdata = (double *)malloc(datasize);//n_data);
  //Initalize d_pdfatdata or BAD THINGS actually I think this is now unnecesary
  for(int i=0;i<n_data;i++) h_pdfatdata[i]=0.0;

  double* h_datamcweights = (double *)malloc(mcDatasize);//n_data);
  
  for(int i=0;i<n_mcData;i++) h_datamcweights[i]=0.0;
  cudaMalloc(&d_datamcweights, mcDatasize);  
  cudaMemcpy(d_datamcweights, h_datamcweights, mcDatasize, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(c_data_x, h_data_x, datasize);

  cudaMalloc(&d_mcweights, mcsize);
  cudaMemcpy(d_mcweights, h_mcweights, mcsize, cudaMemcpyHostToDevice);  
  cudaMalloc(&d_mc_x, mcsize);
  cudaMemcpy(d_mc_x, h_mc_x, mcsize, cudaMemcpyHostToDevice);
  
  calcKde<<<grid_size, block_size>>>(n_data, n_mc, sigma, 
                                     //d_data_x,
                                     d_mc_x, //d_sigmas, 
                                     d_mcweights,d_pdfatdata, 
                                     d_datamcweights);
  grid_size.x = (n_data / block_size.x) + 1;
  calcPdfAtData<<<grid_size, block_size>>>(n_data, thingsForEachData, sigma, 
                                         d_datamcweights, d_pdfatdata);

  cudaMemcpy(h_pdfatdata, d_pdfatdata, datasize, cudaMemcpyDeviceToHost);

  cudaFree(d_mc_x);
  cudaFree(d_mcweights);
  cudaFree(d_datamcweights);
  cudaFree(d_pdfatdata);

  return h_pdfatdata;
}

