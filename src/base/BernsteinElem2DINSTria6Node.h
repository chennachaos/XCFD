
#ifndef incl_BernsteinElem2DINSTria6Node_h
#define incl_BernsteinElem2DINSTria6Node_h

#include "ElementBase.h"
#include "headersBasic.h"
#include "SolutionData.h"


class  BernsteinElem2DINSTria6Node
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

    BernsteinElem2DINSTria6Node();

    ~BernsteinElem2DINSTria6Node();

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

    int calcLoadVector(VectorXd& Flocal);

    double  calcCriticalTimeStep(const double* elemData, const double* timeData, const VectorXd&  veloVec);

    void prepareElemData(const vector<vector<double> >& node_coods);

    int MassMatrices(const vector<vector<double> >& node_coods, const double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2);

    double  ResidualIncNavStokesAlgo1(const vector<vector<double> >& node_coods, const double* elemData, const double* timeData, const VectorXd& veloVec, const VectorXd& veloVecPrev, const VectorXd& veloDotVec, const VectorXd& veloDotVecPrev, const VectorXd& presVec, const VectorXd& presVecPrev, double* FlocalVelo, double* FlocalPres, double timeCur);

    int  ResidualIncNavStokesAlgo2(const vector<vector<double> >& node_coods, const double* elemData, const double* timeData, const VectorXd& veloCur, const VectorXd& veloDotCur, const VectorXd& presCur, VectorXd&  Flocal2);

    double CalculateError(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& veloPrev, const VectorXd& veloDotPrev, const VectorXd& presPrev, double timeCur, int index);

    int  StiffnessAndResidual(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& solnCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur);

    int toComputeInfSupCondition(const vector<vector<double> >& node_coords, const double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp);

    int CalculateForces(int side, const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& veloVec, const VectorXd& presVec, VectorXd&  Flocal1);
};







#endif


