// test the calculation of rotational structure of LiCs molecule
#ifndef TRS_H
#define TRS_H

#include "unit_conversion.hpp"
#include "PolarDiatomicMolecule.hpp"
#include "RotationalStructure.hpp"
#include "dipole_dipole.hpp"
#include <iomanip>

void test_RotationalStructure(){
  //molecule information for LiCs  
  double rot_const=11.7998e9/2.0;  //unit=Hz 1GHz=10^9 Hz
  rot_const = hz_to_au_freq(rot_const);
  double Dipole=5.529; 
  Dipole = debye_to_au(Dipole);
          
  double electric_field = 687630.0; //1.e5; // V/m
  electric_field = v_m_to_au(electric_field);
  //double theta = pi/2
  double alpha_parallel = 597.0; // a.u.
  double alpha_perpendicular = 262.5; // a.u.
  
  PolarDiatomicMolecule LiCs(Dipole, rot_const,alpha_parallel,alpha_perpendicular);
  RotationalStructure LiCs_rot(LiCs, electric_field, 10);
  std::cout << "The eigenvalues are:"<<std::endl;
  std::cout << LiCs_rot.getEigenValues()<<std::endl;
  
  std::cout<<"The first eigenvector is:"<<std::endl;
  std::cout<< LiCs_rot.getNthEigenVector(10)<<std::endl;
  
  std::cout<<"Check if the eigenvector has been normalized."<<std::endl;
  std::cout<<LiCs_rot.getNthEigenVector(1).dot(LiCs_rot.getNthEigenVector(1))<<std::endl;
                          
}


void test_dipole_dipole() {
    //molecule information for LiCs  
  double rot_const=11.7998e9/2.0;  //unit=Hz 1GHz=10^9 Hz
  rot_const = hz_to_au_freq(rot_const);
  double Dipole=5.529; 
  Dipole = debye_to_au(Dipole);
  double Lattice_Constant=4.e-7; 
  Lattice_Constant = meter_to_au(Lattice_Constant);        
  double electric_field = 687630.0; //0.e5; // V/m
  electric_field = v_m_to_au(electric_field);
  //double theta = pi/2
  double alpha_parallel = 597.0; // a.u.
  double alpha_perpendicular = 262.5; // a.u.
  
  PolarDiatomicMolecule LiCs(Dipole, rot_const,alpha_parallel,alpha_perpendicular);
  std::cout<<"Good here 1"<<std::endl;
  DipoleDipole dd(LiCs, Lattice_Constant,electric_field);
  std::cout<<"Good here 2"<<std::endl;
  dd.setUp(2,0,2,0,2,0,2,0);
  //dd.setUp(1,1,1,1,1,1,1,1);
  std::cout<<"Good here 3"<<std::endl;
  double J = dd.getInteraction(1,0,0,0,0,0,1,0);
  //J = dd.getInteraction(2,-2, 2,2,2,2,2,-2);
  std::cout << "<10|<00| V |00>|10> = " << std::setprecision(13) <<  hartree_to_kHz(J) << " kHz" <<std::endl;
  
  double D = dd.getInteraction(1,0,1,0,1,0,1,0)  
              + dd.getInteraction(0,0,0,0,0,0,0,0)
              -2*dd.getInteraction(1,0,0,0,1,0,0,0);
  std::cout << "D = " << std::setprecision(13) <<  hartree_to_kHz(D) << " kHz" <<std::endl;
  
  //for (int i=0; i<40; ++i) {
    //std::cout << "<"<<i<<",0| " << "<"<<i+1<<",0| " << "V_dd" << " |"<< i+1 << ",0>" << " |"<<i<<",0>  =  "
    //<< hartree_to_kHz(dd.getInteraction(i,0,i+1,0,i+1,0,i,0))<<std::endl;
  //}  

}  

#endif
