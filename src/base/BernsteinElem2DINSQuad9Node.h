
#ifndef incl_BernsteinElem2DINSQuad9Node_h
#define incl_BernsteinElem2DINSQuad9Node_h

#include "ElementBase.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "SolutionData.h"


class  BernsteinElem2DINSQuad9Node : public ElementBase
{
  public:

    BernsteinElem2DINSQuad9Node();

    virtual ~BernsteinElem2DINSQuad9Node();

    void prepareElemData(const vector<vector<double> >& node_coods);

    virtual double calcCriticalTimeStep(const double* elemData, const double* timeData, const VectorXd&  veloVec);

    virtual int calcLoadVector(VectorXd& Flocal);

    virtual int MassMatrices(const vector<vector<double> >& node_coods, const double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2);

    virtual double  ResidualIncNavStokesAlgo1(const vector<vector<double> >& node_coods, const double* elemData, const double* timeData, const VectorXd& veloVec, const VectorXd& veloVecPrev, const VectorXd& veloDotVec, const VectorXd& veloDotVecPrev, const VectorXd& presVec, const VectorXd& presVecPrev, VectorXd&  FlocalVelo, VectorXd&  FlocalPres, double timeCur);

    virtual int  ResidualIncNavStokesAlgo2(const vector<vector<double> >& node_coods, const double* elemData, const double* timeData, const VectorXd& veloCur, const VectorXd& veloDotCur, const VectorXd& presCur, VectorXd&  Flocal2);

    virtual double CalculateError(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& veloPrev, const VectorXd& veloDotPrev, const VectorXd& presPrev, double timeCur, int index);

    virtual int  StiffnessAndResidual(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& solnCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur);

    virtual int toComputeInfSupCondition(const vector<vector<double> >& node_coords, const double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp);
};







#endif


