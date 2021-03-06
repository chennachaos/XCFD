
#ifndef incl_SolutionData_h
#define incl_SolutionData_h

#include "headersBasic.h"


class  SolutionData
{
  public:
    int  nNode_Velo, nNode_Pres, ndim;
    int  nDBC, nFBC;
    int  totalDOFvelo, totalDOFpres;
    double  rho, dt;

  public:

    int  tis;

    VectorXd  pres, presCur, presPrev, presPrev2, presPrev3;
    VectorXd  presDot, presDotPrev, presDotCur, presDiff;
    VectorXd  velo, veloCur, veloDiff, veloPrev, veloPrev2, veloPrev3;
    VectorXd  veloDot, veloDotPrev, acceCur;

    VectorXd  td;

    vector<double>  FluidProps;
    vector<int>  node_map_new_to_old;
    vector<int>  node_map_old_to_new;

    vector<vector<double> > DirichletBCs;


    SolutionData();

    ~SolutionData(){}

    void setTimeIncrementType(int ttt)
    {  tis = ttt; }

    void setSpectralRadius(double ttt)
    {  rho = ttt; }

    void  initialise();

    void  setTimeParam();

    void  timeUpdate();

    void  updateIterStep();

    void  solve();

    void  applyDirichletBCs(double fact);

    void  addNodalForces(double fact);

    void  reset();

    void  perform_Aitken_accelerator();

    void  perform_Aitken_accelerator_velocity();

    void  perform_Aitken_accelerator_pressure();

};





#endif


