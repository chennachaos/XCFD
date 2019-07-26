
#ifndef incl_BernsteinElem2DINSTria6Node_h
#define incl_BernsteinElem2DINSTria6Node_h

#include "ElementBase.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "SolutionData.h"


class  BernsteinElem2DINSTria6Node : public ElementBase
{
  public:

    BernsteinElem2DINSTria6Node();

    virtual ~BernsteinElem2DINSTria6Node();

    void prepareElemData(vector<vector<double> >& node_coods);

    virtual double calcCriticalTimeStep(double* elemData, double* timeData, VectorXd&  veloVec);

    virtual int calcLoadVector(VectorXd& Flocal);

    virtual int MassMatrices(vector<vector<double> >& node_coods, double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2);

    //int MassMatrix3nodeDof1(vector<vector<double> >& node_coods, double* elemData, VectorXd&  Mlocal);

    virtual double  ResidualIncNavStokesAlgo1(vector<vector<double> >& node_coods, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& veloVecPrev, VectorXd& veloDotVec, VectorXd& veloDotVecPrev, VectorXd& presVec, VectorXd& presVecPrev, VectorXd&  Flocal1, VectorXd&  Flocal2, double timeCur);

    virtual int  ResidualIncNavStokesAlgo2(vector<vector<double> >& node_coods, double* elemData, double* timeData, VectorXd& veloCur, VectorXd& acceCur, VectorXd& presCur, VectorXd&  Flocal2);

    virtual double CalculateError(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& veloDotPrev, VectorXd& presPrev, double timeCur, int index);

    virtual int  StiffnessAndResidual(vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& solnCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur);

    virtual int toComputeInfSupCondition(vector<vector<double> >& node_coords, double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp);

    virtual int CalculateForces(int side, vector<vector<double> >& node_coords, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& presVec, VectorXd&  Flocal1);
};







#endif

