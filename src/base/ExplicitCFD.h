
#ifndef incl_ExplicitCFD_CLASS_h
#define incl_ExplicitCFD_CLASS_h


#include <string.h>
#include <vector>
#include <fstream>
#include "ElementBase.h"
#include "BernsteinElem2DINSTria6Node.h"

using std::vector;
using std::cout;
using std::string;

class ElementBase;


#define BERNSTEIN_ELEMENT_TYPE_NAMES {\
                                "BernsteinElem2DPoissonTria6Node", \
                                "BernsteinElem2DSolidTria6Node", \
                                "BernsteinElem2DSolidBbarFbarTria6Node", \
                                "BernsteinElem2DSolidMixedTria6Node", \
                                "BernsteinElem3DPoissonTetra10Node", \
                                "BernsteinElem3DSolidTetra10Node", \
                                "BernsteinElem3DSolidBbarFbarTetra10Node", \
                                "BernsteinElem3DSolidMixedTetra10Node", \
                                "BernsteinElem2DPoissonQuad9Node", \
                                "BernsteinElem3DPoissonHex27Node", NULL};


#define BERNSTEIN_ELEMENT_TYPE_NAMES_FACE {\
                                "BernsteinElem2DEdge3Node", \
                                "BernsteinElem3DFaceTria6Node", \
                                "BernsteinElem3DFaceQuad9Node", NULL};


class ExplicitCFD
{
    //private:

    public:

        int  ndim, ndof, fileCount, nElem;
        int  totalDOF_Velo, totalDOF_Pres, totalDOF;
        int  nNode, nNode_Velo, nNode_Pres, nDBC, nDBC_Velo, nDBC_Pres, nFBC, nOutputFaceLoads;
        int  npElem, npElemVelo, npElemPres;
        int  stepsMax, outputFreq;

        std::string  infilename;
        std::ofstream  fout_convdata;

        double  conv_tol, rhoInf, am, gamm1, gamm2, CFL, timeFinal,  dt;
        double  elemData[50], timeData[50];

        vector<vector<double> >  node_coords;               //!< coordinates of the nodes (or control points)
        vector<vector<int> >     midNodeData;               //!< data for processing the middle nodes
        vector<vector<int> >     outputEdges;               //!< data for computing drag/lift forces
        vector<vector<int> >     elemConn;                  //!< element-node connectivity array

        // Velocity node -> element,nodeidx 
        vector<int>              veloNodeConn;                  //!< velocitynode-element connectivity array // npElemVelo*nElem
        vector<int>              veloNodeConnEind;              //!< velocitynode-element connectivity array helper // nNode

        // Pressure node -> element,nodeidx 
        vector<int>              presNodeConn;                  //!< velocitynode-element connectivity array // npElemPres*nElem
        vector<int>              presNodeConnEind;              //!< velocitynode-element connectivity array helper // nNode



        vector<int>  assyForSolnVelo, assyForSolnPres, OutputNodes;
        vector<int>  pressure_nodes, pressure_nodes_map;


        vector<vector<double> >  DirichletBCsVelo;          //!< Dirichlet BCs for velocity
        vector<vector<double> >  DirichletBCsPres;          //!< Dirichlet BCs for pressure
        vector<vector<double> >  NeumannBCs;                //!< Neumann BCs
        vector<vector<double> >  InitialConds;              //!< Initial conditions
        vector<vector<double> >  OutputData;                //!< data for output
        vector<vector<double> >  nodeForcesData;
        vector<vector<double> >  ElemFaceLoadData;


        vector<BernsteinElem2DINSTria6Node>  elems;
        ElementBase  **elemsFaces;


        VectorXd  soln, solnInit, ForceVectorExternal;
        VectorXd  totalForce, totalMoment, centroid;
        VectorXd  pres, presCur, presPrev, presPrev2, presPrev3;
        VectorXd  presDot, presDotPrev, presDotCur, presDiff;
        VectorXd  velo, veloCur, veloDiff, veloPrev, veloPrev2, veloPrev3;
        VectorXd  veloDot, veloDotPrev, veloDotCur;
        VectorXd  veloApplied, presApplied;
        VectorXd  globalMassVelo, globalMassPres;
        VectorXd  rhsVecVelo, rhsVecVelo2, rhsVecPres, rhsVecPres2;
        VectorXd  veloM, veloMp1, veloDotM, veloDotMp1;
        VectorXd  presM, presMp1, presDotM, presDotMp1;

    public:

        ExplicitCFD();

        ~ExplicitCFD();

        ///////////////////////////////////////////////////////////
        //
        // DATA related member functions
        //
        ///////////////////////////////////////////////////////////

        void  readNodes(string& fname);

        void  readElementConnectivity(string& fname);

        void  readDirichletBCs(string& fname);

        void  readForceBCs(string& fname);

        void  readOutputData(string& fname);

        void  printData(int, int);

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  assignBoundaryConditions();

        void  prepareInputData();

        void  printInfo();

        void  readInputData(string& fname);

        void  readControlParameters();

        void  plotGeom(int, bool, int, bool, int*);

        void  applyExternalForces();

        void  writeNodalData();

        void  writeReadResult(int, string&);

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  setSolver(int, int *parm = NULL, bool cIO = false);

        int   prepareMatrixPattern();

        int   calcMassMatrixForExplicitDynamics();

        bool  converged();

        bool  diverging(double);

        void  setTimeParam();

        void  timeUpdate();

        void  updateIterStep();

        void  reset();

        void  addExternalForces(double loadFact);

        void  applyBoundaryConditions(double timeFact);

        void  computeElementErrors(int);

        void  setInitialConditions();

        int   solveExplicitStep();

        int   solveExplicitStepDTS();

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  postProcess();
        void  postProcessVTKlegacy();
};






#endif






