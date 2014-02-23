#ifndef UC_H
#define UC_H

#include <math.h>

//Convert Debye to atomic unit
double debye_to_au(const double debye) {
  return debye*0.3934302014076827;
} //end function debye_to_au  




//Convert meter to atomic unit
double meter_to_au(const double meter) {
  return meter*1.889726133921252e10;
} //end function meter_to_au



//Convert from energy atomic unit  to kHz
double hartree_to_kHz(const double hartree) {
  return hartree*6.57968392072144e12;
} //end function hartree_to_kHz

 
 
//Convert electric field to atomic unit
double v_m_to_au(const double electricField) {
  return electricField*1.944690567144141e-12;
} // end function v_m_to_au



/*
 * convert hz to atomic unit of frequency
 */
double hz_to_au_freq(const double hz) {
  return hz*2.418884324306202e-17;
} //end function hz_to_au_freq


/*
 * convert hz to atomic unit of energy
 */
double hz_to_au_energy(const double hz) {

  static double pi = atan(1.0);
  return 2*pi*hz_to_au_freq(hz);
}  

//Convert second to atomic unit
double second_to_au(const double second) {
  return second*4.1341373374e16;
}



// convert w/cm^2 to atomic unit of light intensity
double w_cm2_to_au(const double w_cm2) {
  return w_cm2/3.50944758e16;
}  

#endif
