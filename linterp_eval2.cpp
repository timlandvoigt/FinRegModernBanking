
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>
#include <cstdarg>
#include <string>
#include <vector>
#include <array>
#include <functional>
#ifdef _CHAR16T
#define CHAR16_T
#endif
#include "mex.h"
#include "linterp.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {

    
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }    

  if (nrhs != 3) {
    mexErrMsgTxt("incorrect number of args. test_linterpC({grid1,grid2,grid3,...},V,{x1,x2,x3,...})");
  }
  
  size_t n_dims = mxGetNumberOfElements(prhs[0]);
  
  assert(n_dims > 0);

  
  // grid vectors
  vector<const double*> grid_list;
  vector<size_t> grid_len_list;
  vector<const double*> xi_list;    
  size_t total_size = 1;
  size_t min_len = mxGetM(prhs[2]);
  double *pts = mxGetPr(prhs[2]);
  for (int i=0; i<n_dims; i++) {				// vectors will be flattened
    mxArray *gridi_m = mxGetCell(prhs[0],i); 
    grid_list.push_back(mxGetPr(gridi_m));
	size_t grid_len = mxGetNumberOfElements(gridi_m);
	grid_len_list.push_back(grid_len);
	total_size *= grid_len;
    xi_list.push_back(pts + i*min_len);
	//min_len = (mxGetNumberOfElements(ptsi_m) < min_len) ? mxGetNumberOfElements(ptsi_m) : min_len;
  }
  
  // F array
  if (total_size != mxGetNumberOfElements(prhs[1])) {
    char pcTemp[1024];
	sprintf(pcTemp, "array sizes do not match. total_size=%d, mxGetNumberOfElements=%d", total_size, mxGetNumberOfElements(prhs[1]));
    mexErrMsgTxt(pcTemp);
  }
  const double *p_F = mxGetPr(prhs[1]);
  
  
  // create result
  plhs[0] = mxCreateDoubleMatrix(1, min_len, mxREAL);
  double *p_result = mxGetPr(plhs[0]);
  
  

  // call interp
  if (n_dims == 1) {
    const int N = 1;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else if (n_dims == 2) {
    const int N = 2;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else if (n_dims == 3) {
    const int N = 3;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else if (n_dims == 4) {
    const int N = 4;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else if (n_dims == 5) {
    const int N = 5;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else if (n_dims == 6) {
    const int N = 6;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else if (n_dims == 7) {
    const int N = 7;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else if (n_dims == 8) {
    const int N = 8;
    InterpMultilinear<N, double, false> interp_obj(grid_list.begin(), grid_len_list.begin(), p_F, p_F + total_size);
    interp_obj.interp_vec(min_len, xi_list.begin(), xi_list.end(), p_result);    
  } else {
    mexErrMsgTxt("dimension not implemented");
  }
}


