#ifndef incl_utilities_cfd_h
#define incl_utilities_cfd_h


#include "headersBasic.h"


int  pointInverseTria6node(double* xNode, double* yNode, double* target, double* param);

int  pointInverseTetra10node(double* xNode, double* yNode, double* zNode, double* target, double* param);


inline void setZero(MatrixXd& AA)
{
    int ii, jj;
    for(ii=0;ii<AA.size();ii++)
    {
        for(jj=0;jj<AA[ii].size();jj++)
           AA[ii][jj] = 0.0;
    }

    return;
}




inline void printMatrix(MatrixXd& AA)
{
    int ii, jj;
    printf("\n\n");
    for(ii=0;ii<AA.size();ii++)
    {
        for(jj=0;jj<AA[ii].size();jj++)
           printf("\t%14.8f", AA[ii][jj]);
        printf("\n");
    }
    printf("\n\n");

    return;
}

inline void setZero(VectorXd& AA)
{
    for(int ii=0;ii<AA.size();ii++)
    {
        AA[ii] = 0.0;
    }

    return;
}




inline void printVector(VectorXd& AA)
{
    printf("\n\n");
    for(int ii=0;ii<AA.size();ii++)
        printf("\t%6d\t%12.8f\n", ii, AA[ii]);
    printf("\n\n");

   return;
}



inline void printVector(vector<int>&  vec)
{
    printf("\n\n");
    for(int ii=0;ii<vec.size();ii++)
        printf("\t%6d ", vec[ii]);
    printf("\n\n");

   return;
}


inline void printVector(double* data, int nn)
{
    printf("\n\n");
    for(int ii=0;ii<nn;ii++)
      printf("\t%6d\t%12.8f\n", ii, data[ii]);
    printf("\n\n");

   return;
}


double  norm(vector<double>&  vec);

void  subtract2vecs(const vector<double>&  vec1, const vector<double>&  vec2, vector<double>&  vec);




#endif
