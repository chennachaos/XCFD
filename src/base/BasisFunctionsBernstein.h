
#ifndef incl_BasisFunctionsBernstein_h
#define incl_BasisFunctionsBernstein_h


/**
  Computes function values and first derivative of
  univariate Bernstein polynomials of a given order p
*/
void Bernstein_BasisFuns1D(int p, double xi, double* N, double* dN_dxi);

/**
  Computes function values of
  univariate Bernstein polynomials of a given order p
*/
void Bernstein_BasisFuns1D(int p, double xi, double* N);


/**
  Computes function values of
  Bernstein polynomials of a given order p for Triangular elements
*/
void BernsteinBasisFunsTria(int p, double xi, double zeta, double* N);

/**
  Computes function values and first derivative of
  Bernstein polynomials of a given order p for Triangular elements
*/
void BernsteinBasisFunsTria(int p, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);


/**
  Computes function values and first derivative of
  Bernstein polynomials of a given order p for Quadrilateral elements
*/
void BernsteinBasisFunsQuad(int p, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);


/**
  Computes function values and first derivative of
  Bernstein polynomials of a given order p for Tetrahedral elements
*/
void BernsteinBasisFunsTetra(int p, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);




/**
  Computes function values of
  Bernstein polynomials of an edge in 2D.
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void BernsteinBasisFunsEdge2D(int p, double* param, double *xNode, double* yNode, double *N, double *normal, double& Jac);

/**
  Computes function values of
  Bernstein polynomials of a face in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void BernsteinBasisFunsFace3D(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac);


/**
  Computes function values of
  Bernstein polynomials of a triangular face in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void BernsteinBasisFunsFaceTria(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac);


/**
  Computes function values of
  Bernstein polynomials of a quadrilateral face in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void BernsteinBasisFunsFaceQuad(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac);


int computeBasisFunctions2D(bool flag, int ETYPE, int degree, double* param, double* xNode, double* yNode, double*  N, double*  dN_dx, double* dN_dy, double&  Jac);


int computeBasisFunctions3D(bool flag, int ETYPE, int degree, double* param, double* xNode, double* yNode, double* zNode, double*  N, double*  dN_dx, double* dN_dy, double* dN_dz, double&  Jac);


#endif
