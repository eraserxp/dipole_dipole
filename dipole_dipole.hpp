// calculate the dipole-dipole interaction between two molecules
#ifndef DD_H
#define DD_H

#include <boost/make_shared.hpp>
#include "mics.hpp"
#include "PolarDiatomicMolecule.hpp"
#include "RotationalStructure.hpp"
#include <vector>
#include <algorithm> //std::sort and unique and distance
#define EIGEN_USE_MKL_ALL
#include <Eigen/Dense>
using namespace Eigen;

//two same molecules in the same DC field
class DipoleDipole {
  private:
    PolarDiatomicMolecule molecule;
    boost::shared_ptr<RotationalStructure>  rs_tmp;
    boost::shared_ptr<RotationalStructure>  rsA_left; // for state |NA, MA>
    boost::shared_ptr<RotationalStructure>  rsB_left; // for state |NB, MB>
    boost::shared_ptr<RotationalStructure>  rsA_right; // for state |NA_, MA_>
    boost::shared_ptr<RotationalStructure>  rsB_right; // for state |NB_, MB_>
    double distance;
    double dipoleA;
    double dipoleB;
    double electricField;
    
  
  public:
    //constructor
    // two molecules are the same and they are in the same DC field
    DipoleDipole(const PolarDiatomicMolecule& moleculeA, 
             const double distance, const double electricField);
    
    // calculate the rotational structures of the two molecules
    void setUp(const int NA, const int MA, const int NB, const int MB,
               const int NA_, const int MA_, const int NB_, const int MB_);
             
    // obtain the interaction between two molecules 
    // <NA,MA| <NB,MB| V |NA_,MA_> |NB_,MB_>
    double getInteraction(const int NA, const int MA, 
                          const int NB, const int MB,
                          const int NA_, const int MA_, 
                          const int NB_, const int MB_);
                          
    // calculate <NA,MA| <NB,MB| V |NA_,MA_> |NB_,MB_> without external fields
    double ddNoField(const int NA, const int MA, 
                     const int NB, const int MB,
                     const int NA_, const int MA_, 
                     const int NB_, const int MB_);
};


/*
 * implementation
 * 
 * 
 */


 
DipoleDipole::DipoleDipole(const PolarDiatomicMolecule& moleculeA, 
             const double distance, const double electricField) {
               
  this->molecule = moleculeA;
  dipoleA = molecule.getDipole();
  dipoleB = dipoleA;
  this->distance = distance;
  this->electricField = electricField;
}


void DipoleDipole::setUp(const int NA, const int MA, 
                         const int NB, const int MB,
                         const int NA_, const int MA_, 
                         const int NB_, const int MB_) {
                           
  if ( (NA-MA>50)||(NB-MB>50)||(NA_-MA_>50)||(NB_-MB_>50) ) {
    std::cout<<"The input rotational levels are out of range!"<<std::endl;
    abort();
  }; 
  
   
  std::vector<int> absolute_m;
  absolute_m.push_back(abs(MA));
  absolute_m.push_back(abs(MB));
  absolute_m.push_back(abs(MA_));
  absolute_m.push_back(abs(MB_));
   
  // remove duplicate elements in absolute_m;
  std::sort(absolute_m.begin(), absolute_m.end());
  std::vector<int>::iterator it = std::unique(absolute_m.begin(),absolute_m.end());
  absolute_m.erase(it, absolute_m.end());
  
      
  for (std::vector<int>::iterator it = absolute_m.begin() ; it != absolute_m.end(); ++it) {
    // there is only one rotational structure for each unique m
    rs_tmp = boost::make_shared<RotationalStructure>(molecule, electricField, *it); 

    if (abs(MA)==*it) 
      {std::cout<< "MA" << std::endl;
      rsA_left = rs_tmp;}
    if (abs(MB)==*it)
      {std::cout<< "MB" << std::endl;
      rsB_left = rs_tmp;}
    if (abs(MA_)==*it) 
      {std::cout<< "MA_" << std::endl;
      rsA_right = rs_tmp;}
    if (abs(MB_)==*it) 
      {std::cout<< "MB_" << std::endl;
      rsB_right = rs_tmp;}
  }; 
 
  
}


double DipoleDipole::getInteraction(const int NA, const int MA, 
                                    const int NB, const int MB,
                                    const int NA_, const int MA_, 
                                    const int NB_, const int MB_){
                                    
  double result = 0.0;
  // find out |NA, MA> belongs to the nth eigenvector
  int nthA_left = NA - abs(MA) + 1;
  int nthB_left = NB - abs(MB) + 1;
  int nthA_right = NA_ -abs(MA_) + 1;
  int nthB_right = NB_ - abs(MB_) + 1;
  
  
  // get the corrsponding eigenvectors
  
  VectorXd vec_A_left = rsA_left->getNthEigenVector(nthA_left);
  
  VectorXd vec_B_left = rsB_left->getNthEigenVector(nthB_left);
  
  VectorXd vec_A_right = rsA_right->getNthEigenVector(nthA_right);
  
  VectorXd vec_B_right = rsB_right->getNthEigenVector(nthB_right);
 
    
  int iStart = abs(MA);
  int iEnd = iStart + rsA_left->getNumStates()-1;

  int jStart = abs(MB);
  int jEnd = jStart + rsB_left->getNumStates()-1;
 
  int kStart = abs(MA_);
  int kEnd = kStart + rsA_right->getNumStates()-1;

  int lStart = abs(MB_);
  int lEnd = lStart + rsB_right->getNumStates()-1;

  
  int i, j, k, l;
  
  #pragma omp parallel for private(i,j,k,l) shared(result)
  for ( i=iStart; i<=iEnd; ++i) {
    for ( j=jStart; j<=jEnd; ++j) {
      for ( k=kStart; k<=kEnd; ++k) {
        for ( l=lStart; l<=lEnd; ++l) {
          #pragma omp atomic
          result+=vec_A_left(i-iStart)*vec_B_left(j-jStart)
                    *vec_A_right(k-kStart)*vec_B_right(l-lStart)
                    *ddNoField(i,MA,j,MB,k,MA_,l,MB_);          
        }
      }
    }
  }
  
  return result;  
}



double DipoleDipole::ddNoField(const int NA, const int MA, 
                               const int NB, const int MB,
                               const int NA_, const int MA_, 
                               const int NB_, const int MB_) {

  double prefactor=-3.0*sqrt(5.0)*pow(-1.0, -(MA+MB))        
                   *( dipoleA*dipoleB/pow(distance,3) )    
            *sqrt( (2.0*NA+1.0)*(2.0*NA_+1.0)*(2.0*NB+1.0)*(2.0*NB_+1.0) )     
            *wigner_3j(NA, 1,NA_,0,0,0)*wigner_3j(NB, 1,NB_,0,0,0);     
  
  return prefactor*( 
                    0.5*AplusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)      
                     + 0.5*AminusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)     
                     - sqrt(1.0/6.0)*AplusBminus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)     
                     - sqrt(1.0/6.0)*AminusBplus(NA,MA,NB,MB,NA_,MA_,NB_,MB_)     
                     - sqrt(1.0/6.0)*A0B0(NA,MA,NB,MB,NA_,MA_,NB_,MB_)  
                     )*pow(-1.0, -(MA+MB));                    
}
                         
#endif
