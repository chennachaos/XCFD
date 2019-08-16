#ifndef incl_LagrangeElement_h
#define incl_LagrangeElement_h

#include "headersBasic.h"
//#include <vector>
//#include <iostream>

using std::vector;
using std::cout;


enum MeshConfiguration
{
  CONFIG_ORIGINAL=0, CONFIG_DEFORMED=1
};


enum ElementShapes
{
  ELEM_SHAPE_TRIA=1, ELEM_SHAPE_QUAD=2, 
  ELEM_SHAPE_TETRA=4, ELEM_SHAPE_PYRAMID=5, ELEM_SHAPE_PENTA=6, ELEM_SHAPE_HEXA=8,
  ELEM_SHAPE_TRIA_BERNSTEIN=100, ELEM_SHAPE_TETRA_BERNSTEIN=400,
  ELEM_SHAPE_QUAD_BERNSTEIN=101, ELEM_SHAPE_HEXA_BERNSTEIN=401
};



class ElementBase
{
  public:

    //member variables
    int  ELEM_TYPE;
    int  ndim, degree, npElem, nlbf, ndof, nsize, nGP;
    int  elmType, matType, secType, subdomId;
    double  elemVol, charlen, dtCrit;

    vector<int>  nodeNums, forAssyVec;

    vector<double>  elemVolGP;

    vector<VectorXd >  Nv, dNvdx, dNvdy, dNvdz, Np, dNpdx, dNpdy, dNpdz;


    //member functions

    ElementBase();

    virtual ~ElementBase();

    int getDimension()
    { return ndim; }

    int getPolynomialDegree()
    { return degree; }

    void  setSubdomainId(int sid)
    {  subdomId = sid; return;  }

    int getSubdomainId()
    {  return  subdomId;  }

    int getNodesPerElement()
    {  return  npElem; }

    int getNdofPerNode()
    { return ndof; }

    int  getNdofPerElement()
    {  return  nsize;  }

    std::vector<int>&  getNodeNumbers()
    {  return  nodeNums; }

    std::vector<int>&  getVectorForAssembly()
    {  return  forAssyVec; }

    double getVolume()
    { return  elemVol;  }

    virtual int calcLoadVector(VectorXd& Flocal)
    { cout << "   'calcLoadVector' is not defined for this element!\n\n"; return 0; }

    virtual double  calcCriticalTimeStep(const double* elemData, const double* timeData, const VectorXd&  veloVec)
    { cout << "   'calcCriticalTimeStep' is not defined for this element!\n\n"; return 0.0; }

    virtual  int calcError(int index)
    { cout << "  'calcError' is not available for this element!\n\n"; return -1; }

    virtual void prepareElemData(const vector<vector<double> >& node_coods)
    { cout << "   'prepareElemData' is not defined for this element!\n\n"; return; }

    virtual int MassMatrices(const vector<vector<double> >& node_coods, const double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2)
    { cout << "   'MassMatrices' is not defined for this element!\n\n"; return -1; }

    virtual double  ResidualIncNavStokesAlgo1(const vector<vector<double> >& node_coods, const double* elemData, const double* timeData, const VectorXd& veloVec, const VectorXd& veloVecPrev, const VectorXd& veloDotVec, const VectorXd& veloDotVecPrev, const VectorXd& presVec, const VectorXd& presVecPrev, VectorXd&  FlocalVelo, VectorXd&  FlocalPres, double timeCur)
    { cout << "   'ResidualIncNavStokesAlgo1' is not defined for this element!\n\n"; return -1; }

    virtual int  ResidualIncNavStokesAlgo2(const vector<vector<double> >& node_coods, const double* elemData, const double* timeData, const VectorXd& veloCur, const VectorXd& veloDotCur, const VectorXd& presCur, VectorXd&  Flocal2)
    { cout << "   'ResidualIncNavStokesAlgo2' is not defined for this element!\n\n"; return -1; }

    virtual double CalculateError(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& veloPrev, const VectorXd& veloDotPrev, const VectorXd& presPrev, double timeCur, int index)
    { cout << "   'CalculateError' is not defined for this element!\n\n"; return 0; }

    virtual int  StiffnessAndResidual(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& solnCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur)
    { cout << "   'StiffnessAndResidual' is not defined for this element!\n\n"; return -1; }

    virtual int toComputeInfSupCondition(const vector<vector<double> >& node_coords, const double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp)
    { cout << "   'toComputeInfSupCondition' is not defined for this element!\n\n"; return -1; }

    virtual int CalculateForces(int side, const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& veloVec, const VectorXd& presVec, VectorXd&  Flocal1)
    { cout << "   'CalculateForces' is not defined for this element!\n\n"; return -1; }

};



#endif //incl_Lagrange_Element_h

