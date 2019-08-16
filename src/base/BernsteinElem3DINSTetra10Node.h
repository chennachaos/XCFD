
#ifndef incl_BernsteinElem3DINSTetra10Node_h
#define incl_BernsteinElem3DINSTetra10Node_h

#include "ElementBase.h"
#include "headersBasic.h"
#include "SolutionData.h"

class  BernsteinElem3DINSTetra10Node : public ElementBase
{
  public:

    BernsteinElem3DINSTetra10Node();

    virtual ~BernsteinElem3DINSTetra10Node();

    void prepareElemData(const vector<vector<double> >& node_coods);

    virtual double calcCriticalTimeStep(const double* elemData, const double* timeData, const VectorXd&  veloVec);

    virtual int calcLoadVector(VectorXd& Flocal);

    virtual int MassMatrices(const vector<vector<double> >& node_coods, const double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2);

    virtual double  ResidualIncNavStokesAlgo1(const vector<vector<double> >& node_coods, const double* elemData, const double* timeData, const VectorXd& veloVec, const VectorXd& veloVecPrev, const VectorXd& veloDotVec, const VectorXd& veloDotVecPrev, const VectorXd& presVec, const VectorXd& presVecPrev, VectorXd&  FlocalVelo, VectorXd&  FlocalPres, double timeCur);

    virtual int  ResidualIncNavStokesAlgo2(const vector<vector<double> >& node_coods, const double* elemData, const double* timeData, const VectorXd& veloCur, const VectorXd& veloDotCur, const VectorXd& presCur, VectorXd&  Flocal2);

    virtual double CalculateError(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& veloPrev, const VectorXd& veloDotPrev, const VectorXd& presPrev, double timeCur, int index);

    virtual int  StiffnessAndResidual(const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& solnCur, MatrixXd& Klocal, VectorXd& Flocal, double timeCur);

    virtual int toComputeInfSupCondition(const vector<vector<double> >& node_coords, const double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp);

    virtual int CalculateForces(int side, const vector<vector<double> >& node_coords, const double* elemData, const double* timeData, const VectorXd& veloVec, const VectorXd& presVec, VectorXd&  Flocal1);
};




#endif


