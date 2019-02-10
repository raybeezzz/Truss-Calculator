#include "pch.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>

using namespace std;

//----- PROGRAM CONSTANTS ---------------------------------------------------------------
#define X         0                     // x-direction
#define Y         1                     // y-direction
#define MAXLINE   1024                  // string buffer size

// headers for data file parsing
#define JOINT_COORDINATE_HEADER      "[JOINT COORDINATES]"
#define MEMBER_CONNECTIVITY_HEADER   "[MEMBER JOINT CONNECTIVITY]"
#define REACTIONS_HEADER             "[REACTIONS AT NODES]"
#define EXTERNAL_FORCES_HEADER       "[EXTERNAL FORCES]"
#define FORCE_UNITS_HEADER           "[FORCE UNITS]"

//----- STRUCTURE DEFINITIONS -----------------------------------------------------------

typedef struct MEMBER
{
   int j[2];      // joint indexes
   double F, L;   // internal Force (+ve is tension), member Length
}
MEMBER;

typedef struct JOINT
{
   double p[2];  // joint x,y position.  p[0]=x, p[1]=y
}
JOINT;

typedef struct FORCE  // for external and reaction forces
{
   double F;      // Force
   int j;         // joint index 
   int idir;      // if force is given as X or Y
   double theta;  // if force is given as angle from positive x axis

   FORCE()  // initialize with bad values
   {
      F=1.0e30;
      j=-1, idir=-1;
      theta=1.0e30;
   };
}
FORCE;

typedef struct TRUSS_SIZES
{
   int NJ, NM, NR, NEF;  // number of Joints, Members, Reactions, External Forces
}
TRUSS_SIZES;

//----- FUNCTION PROTOTYPES ----------------------------------------------------------------------------

// parse the file to get number of joints,members,reactions,forces
TRUSS_SIZES getTrussSizes();  

// gets the truss data from a file
void getTrussData(TRUSS_SIZES ts, vector<JOINT> Jf, vector<MEMBER> Mf, vector<FORCE> Ef, vector<FORCE> Rf, string strUnits);

// assemble the matrices and vectors
bool buildSystem();

// function to invert an NxN matrix
bool InverseNxN(double **A, int N, double **Ainv);

// helper function for InverseNxN
void ELGS(double **A, int N, int *indx);
