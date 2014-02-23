#ifndef RS_H
#define RS_H

#include "mics.hpp"
#include "PolarDiatomicMolecule.hpp"
#define EIGEN_USE_MKL_ALL
#include <Eigen/Dense>
using namespace Eigen;

class RotationalStructure {
private:
	PolarDiatomicMolecule molecule;
	double ElectricField;
	// highest rotational levels included in the calculation
	int nmax;
	int rows;
	int columns;
	int absolute_m;
  //the hamiltonian matrix in the bare rotational basis states
	MatrixXd Ham;
	// points to the array of eigen values in atomic unit
	VectorXd eigenValues;
	// eigen vectors
	MatrixXd eigenVectors;


  void formMatrix();
    
  // diagonalize the hamiltonian to get eigenvalues and eigenvectors      
  void solveHam();



public:
  // no-parameter constructor
  RotationalStructure();
  
	//constructor
	RotationalStructure(PolarDiatomicMolecule molecule, const double ElectricField,
			                const int absolute_m, int nmax=50);
  //copy constructor
  RotationalStructure(const RotationalStructure& other);
  
  //assignment operator
  const RotationalStructure& operator =(const RotationalStructure& other);
  
  //return the array of eigen values
  const VectorXd getEigenValues();

  //return the nth eigen vector
  const VectorXd getNthEigenVector(const int n);
  
  //get the number of states included for the calculations
  int getNumStates();

  //destructor
  ~RotationalStructure();

};




/*
 * implementations
 */


void RotationalStructure::formMatrix() {
  double pi = atan(1.0)*4;
  int NStart = absolute_m;
  int NEnd = nmax + absolute_m - 1;
  // off-diagonal part 
  #pragma omp parallel for  
  for (int i=NStart; i<=NEnd; ++i) {
    for (int j=NStart; j<=i; ++j) {
      Ham(i-NStart,j-NStart) = -molecule.getDipole()*pow(-1.0,absolute_m)*ElectricField
                    *sqrt((2.0*i + 1)*(2.0*j + 1))
                    *wigner_3j(i,1,j,-absolute_m,0,absolute_m)
                    *wigner_3j(i,1,j,0,0,0);
      Ham(j-NStart,i-NStart) = Ham(i-NStart,j-NStart);
    }
  }
    
  // diagonal part
  #pragma omp parallel for  
  for (int i=NStart; i<=NEnd; ++i) {
    Ham(i-NStart,i-NStart) += 2*pi*molecule.getRotationalConstant()*(i+1)*i;
  }
  std::cout<<"Form Hamiltonian matrix successfully!"<<std::endl;
}
    
// diagonalize the hamiltonian to get eigenvalues and eigenvectors      
void RotationalStructure::solveHam() {
  formMatrix();  
  SelfAdjointEigenSolver<MatrixXd> eigensolver(Ham);
  if (eigensolver.info()!=Success) {
    abort(); 
  } else {
    std::cout<<"Hamiltonian matrix is diagonalized successfully!"<<std::endl;
  }
  // convert Eigen vector and matrix types to array
  eigenValues = eigensolver.eigenvalues();
  eigenVectors = eigensolver.eigenvectors();
}

//default constructor, do nothing
RotationalStructure::RotationalStructure() {}

//constructor
// nmax default to 50 (see header file)
RotationalStructure::RotationalStructure(PolarDiatomicMolecule molecule, const double ElectricField,
                    const int absolute_m, int nmax) {
  this->molecule = molecule;
  this->ElectricField = ElectricField;
  this->nmax = nmax;
  this->absolute_m = absolute_m;
  rows = nmax;
  columns = nmax;
  Ham.resize(rows,columns);
  //solve the eigen value problem
  solveHam();
}

//copy constructor
RotationalStructure::RotationalStructure(const RotationalStructure& other) {
  molecule = other.molecule;
  ElectricField = other.ElectricField;
  nmax = other.nmax;
  absolute_m = other.absolute_m;
  rows = other.rows;
  columns = other.columns;  
  Ham = other.Ham;
	eigenValues = other.eigenValues;
	eigenVectors = other.eigenVectors;
}

//assignment operator
const RotationalStructure& RotationalStructure::operator =(const RotationalStructure& other){
  molecule = other.molecule;
  ElectricField = other.ElectricField;
  nmax = other.nmax;
  absolute_m = other.absolute_m;
  rows = other.rows;
  columns = other.columns;  
  Ham = other.Ham;
	eigenValues = other.eigenValues;
	eigenVectors = other.eigenVectors;
  return *this;  
}  

//return the pointer that points to the array of eigen values
const VectorXd RotationalStructure::getEigenValues(){
  return eigenValues;
}

//return the pointer that points to the nth eigen vector
const VectorXd RotationalStructure::getNthEigenVector(const int n){
  return eigenVectors.col(n-1);
}

int RotationalStructure::getNumStates() {
  return nmax;
}

RotationalStructure::~RotationalStructure() {
	Ham.resize(0,0);
	eigenValues.resize(0);
	eigenVectors.resize(0,0);
}

#endif
