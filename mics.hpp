#ifndef MICS_H
#define MICS_H

//#include<gsl/gsl_sf_coupling.h>

#include "RotationalStructure.hpp"
#define EIGEN_USE_MKL_ALL
#include <Eigen/Dense>
using namespace Eigen;

extern"C" {
double thrj_(double *F1, double *F2, double *F3, 
             double *G1, double *G2, double *G3);
}


//double wigner_3j(const int j1, const int j2, const int j3, 
                 //const int m1, const int m2, const int m3, bool useFortran=false);


/* 3j Symbols:  / ja jb jc \
 *              \ ma mb mc /
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
double wigner_3j(const int j1, const int j2, const int j3, 
                 const int m1, const int m2, const int m3) {
//		if (useFortran) {
			double f1 = j1;
			double f2 = j2;
			double f3 = j3;
			double g1 = m1;
			double g2 = m2;
			double g3 = m3;
			return thrj_(&f1, &f2, &f3, &g1, &g2, &g3);
//		}	
//    return gsl_sf_coupling_3j(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3);
};




/*              /  1  1 2 \  /  NA  1 NA_ \  /  NB  1 NB_ \
 *              \ -1 -1 2 /  \ -MA -1 MA_ /  \ -MB -1 MB_ /
 */
double AminusBminus(int NA, int MA,int NB,int MB,int NA_,int MA_, int NB_, int MB_) {
  return wigner_3j(1,1,2,-1,-1,2)*wigner_3j(NA,1,NA_,-MA,-1,MA_)*
         wigner_3j(NB,1,NB_,-MB,-1,MB_);
};


        
/*              /  1  1  2 \  /  NA  1 NA_ \  /  NB  1 NB_ \
 *              \  1  1 -2 /  \ -MA  1 MA_ /  \ -MB  1 MB_ /
 */
double AplusBplus(int NA, int MA,int NB,int MB,int NA_,int MA_, int NB_, int MB_) {
  return wigner_3j(1,1,2,1,1,-2)*wigner_3j(NA,1,NA_,-MA,1,MA_)*
         wigner_3j(NB,1,NB_,-MB,1,MB_);
};

  

/*              /  1  1 2 \  /  NA  1 NA_ \  /  NB  1 NB_ \
 *              \  1 -1 0 /  \ -MA  1 MA_ /  \ -MB -1 MB_ /
 */
double AplusBminus(int NA, int MA,int NB,int MB,int NA_,int MA_, int NB_, int MB_) {
  return wigner_3j(1,1,2,1,-1,0)*wigner_3j(NA,1,NA_,-MA,1,MA_)*
         wigner_3j(NB,1,NB_,-MB,-1,MB_);
};

   

/*              /  1  1 2 \  /  NA  1 NA_ \  /  NB  1 NB_ \
 *              \ -1  1 0 /  \ -MA -1 MA_ /  \ -MB  1 MB_ /
 */
double AminusBplus(int NA, int MA,int NB,int MB,int NA_,int MA_, int NB_, int MB_) {
  return wigner_3j(1,1,2,-1,1,0)*wigner_3j(NA,1,NA_,-MA,-1,MA_)*
         wigner_3j(NB,1,NB_,-MB,1,MB_);
};

  

/*              /  1  1 2 \  /  NA  1 NA_ \  /  NB  1 NB_ \
 *              \  0  0 0 /  \ -MA  0 MA_ /  \ -MB  0 MB_ /
 */
double A0B0(int NA, int MA,int NB,int MB,int NA_,int MA_, int NB_, int MB_) {
  return wigner_3j(1,1,2,0,0,0)*wigner_3j(NA,1,NA_,-MA,0,MA_)*
         wigner_3j(NB,1,NB_,-MB,0,MB_);
};
  
    

/*              /  1  1 2 \  /  NA  1 NA_ \  /  NB  1 NB_ \
 *              \  0 -1 1 /  \ -MA  0 MA_ /  \ -MB -1 MB_ /
 */
double A0Bminus(int NA, int MA,int NB,int MB,int NA_,int MA_, int NB_, int MB_) {
  return wigner_3j(1,1,2,0,-1,1)*wigner_3j(NA,1,NA_,-MA,0,MA_)*
         wigner_3j(NB,1,NB_,-MB,-1,MB_);
};

     

/*              /  1  1  2 \  /  NA  1 NA_ \  /  NB  1 NB_ \
 *              \  0  1 -1 /  \ -MA  0 MA_ /  \ -MB  1 MB_ /
 */
double A0Bplus(int NA, int MA,int NB,int MB,int NA_,int MA_, int NB_, int MB_) {
  return wigner_3j(1,1,2,0,1,-1)*wigner_3j(NA,1,NA_,-MA,0,MA_)*
         wigner_3j(NB,1,NB_,-MB,1,MB_);
};
      

      
/*              /  1  1 2 \  /  NA  1 NA_ \  /  NB  1 NB_ \
 *              \ -1  0 1 /  \ -MA -1 MA_ /  \ -MB  0 MB_ /
 */
double AminusB0(int NA, int MA,int NB,int MB,int NA_,int MA_, int NB_, int MB_) {
  return wigner_3j(1,1,2,-1,0,1)*wigner_3j(NA,1,NA_,-MA,-1,MA_)*
         wigner_3j(NB,1,NB_,-MB,0,MB_);
};

 


/*              /  1  1  2 \  /  NA  1 NA_ \  /  NB  1 NB_ \
 *              \  1  0 -1 /  \ -MA  1 MA_ /  \ -MB  0 MB_ /
 */
double AplusB0(int NA, int MA,int NB,int MB,int NA_,int MA_, int NB_, int MB_) {
  return wigner_3j(1,1,2,1,0,-1)*wigner_3j(NA,1,NA_,-MA,1,MA_)*
         wigner_3j(NB,1,NB_,-MB,0,MB_);
};


#endif
