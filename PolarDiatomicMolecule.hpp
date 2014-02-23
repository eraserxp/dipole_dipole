#ifndef PDM_H
#define PDM_H

class PolarDiatomicMolecule {
  private:
  // everything in atomic unit
    double dipole;
    double rotationalConstant;
    // polarizability
    double alphaParallel;
    double alphaPerpendicular;
    
  public:
    //constructor with no parameters
    PolarDiatomicMolecule();

    //constructor
    PolarDiatomicMolecule(const double dipole, const double rotationalConstant,
                          const double alphaParallel, const double alphaPerpendicular);
    
    // copy constructor
    PolarDiatomicMolecule(const PolarDiatomicMolecule & other);
    
    //assignment operator
    const PolarDiatomicMolecule& operator=(const PolarDiatomicMolecule& other);

    double getDipole();
    
    double getRotationalConstant();
    
    double getAlphaParallel();
    
    double getAlphaPerpendicular();
    
};



/*
 * implementations
 */

//constructor with no parameters
PolarDiatomicMolecule::PolarDiatomicMolecule() {}

//constructor
PolarDiatomicMolecule::PolarDiatomicMolecule(const double dipole, const double rotationalConstant,
                      const double alphaParallel, const double alphaPerpendicular){
  this->dipole = dipole;
  this->rotationalConstant = rotationalConstant;
  this->alphaParallel = alphaParallel;
  this->alphaPerpendicular = alphaPerpendicular;
}

// copy constructor
PolarDiatomicMolecule::PolarDiatomicMolecule(const PolarDiatomicMolecule & other) {
  dipole = other.dipole;
  rotationalConstant = other.rotationalConstant;
  alphaParallel = other.alphaParallel;
  alphaPerpendicular = other.alphaPerpendicular;
}

//every time you need to write your own copy constructor, you need to write your own assignment operator
//assignment operator
const PolarDiatomicMolecule& PolarDiatomicMolecule::operator=(const PolarDiatomicMolecule& other) {
  dipole = other.dipole;
  rotationalConstant = other.rotationalConstant;
  alphaParallel = other.alphaParallel;
  alphaPerpendicular = other.alphaPerpendicular;
  return *this;
}

double PolarDiatomicMolecule::getDipole() {
  return dipole;
}

double PolarDiatomicMolecule::getRotationalConstant() {
  return rotationalConstant;
}

double PolarDiatomicMolecule::getAlphaParallel() {
  return alphaParallel;
}

double PolarDiatomicMolecule::getAlphaPerpendicular() {
  return alphaPerpendicular;
}

#endif
