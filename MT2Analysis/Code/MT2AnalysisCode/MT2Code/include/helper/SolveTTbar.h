#ifndef SolveTTbar_h
#define SolveTTbar_h


//*************************************************************
//
// \class SolveTTbar: Computes a discriminant to distinguish between ttbar and stop direct production
//
//*************************************************************
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

#include "TH2.h"


using namespace std;


// ************************************************************
// class SolveTTbar
//
// The class SolveTTbar contains several functions attempting to find the minimum of the discriminant
//
// - ScanEqn: Scans the plane (qz1, qz2) for a minimum of the discriminant, which is returned
//		The search window and the number of bins can be user defined
//		This function is quite effective in finding the minimum, but is rather slow (many points)
//
// - IterDeriv: Iterates in the plane (qz1, qz2) for a minimum of the discriminant which is returned
//		Uses the analytical form of the derivatives
//		Much faster (fewer steps) than ScanEqn, but fails more often to find a true minimum
//
// - IterEqn: Iterates in the plane (qz1, qz2) for a minimum of the discriminant which is returned
//		Uses derivatives computed analytically
//		Also fast but worse than IterDeriv, would need more work for improvements
//
// ************************************************************


class SolveTTbar {

public:

// Constructor
// constFileName = name of the file of constants (character string)
SolveTTbar ();

// Destructor
~SolveTTbar();

// Compute the discriminator *************************************

// by scanning over the (qz1, qz2) plane
virtual double ScanEqn (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], TH2D * hdiscr);

// To find all solutions by scanning
virtual int ScanFndAll (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[], 
		double discsol[], double qz1sol[], double qz2sol[]);

// by iteration, using analytical derivative (faster)
virtual double IterDeriv (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[]);

// by iteration, using numerically computed derivatives
virtual double IterEqn (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[]);

// Less often needed functions:
// Solves the constraint equations and returns a discriminant of type D2
virtual double SolveEqn_qz12 (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double q1[], double q2[]);

// Solves the constraint equations and returns a discriminant of type D2 and its derivatives
virtual double DerivSolve (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double q1[], double q2[], double derDisc[]);

// Scans the plane (qx1, qy1) for a minimum of the discriminant which is returned
double ScanEqn_qxy (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[], TH2D * hdiscrxy);

// Scans the plane (qx1, qy1) with binary search for a minimum of the discriminant which is returned
double BinSrc_qxy (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[]);

// Computes the solutions for qy for a given qx
int qyfrqx (double p[], double m, double mWt, double qx, double qy[]);

// Solves the constraint equations and returns a discriminant of type D2
// expects qx1 and qy1 in q1[] to be filled, computes qz1 and qz2
double SolveEqn_qxy1 (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double q1[], double q2[]);

// Solves the 4th orderconstraint equations and returns a discriminant of type D2
// returns the number of solutions
int SolveEqn4 (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[]);

// Setting constants *************************************

// Set the debug flag
void SetDebug (int dbg) {
	debug = dbg;
}

// Set the top mass (default = 173)
void Setmt (double mtin) {
	mt = mtin;
}

// Set the W mass (default = 80)
void SetmW (double mWin) {
	mW = mWin;
}

// Set the neutrino mass (default = 0)
void Setmnu (double mnuin) {
	mnu = mnuin;
}

// Set the lower and higher limts on the x-axis window for scan with ScanEqn
void SetXWindow (double xl, double xh) {
	xlo = xl;
	xhi = xh;
}

// Set the lower and higher limts on the y-axis window for scan with ScanEqn
void SetYWindow (double yl, double yh) {
	ylo = yl;
	yhi = yh;
}

// Set the number of bins (or steps) for scan in qz1 and qz2 with ScanEqn or IterDeriv
void SetNscan (int nsc) {
	nscan = nsc;
}

// Set the precison on the discriminant for convergence of the scan in qz1 and qz2 with ScanEqn or IterDeriv
void SetDiscrConv (double dconv) {
	discconverge = dconv;
}

// Set the starting value of qz1 for the search with IterDeriv or IterEqn
void SetStartqz1 (double qz) {
	qz1strt = qz;
}

// Set the starting value of qz2 for the search with IterDeriv or IterEqn
void SetStartqz2 (double qz) {
	qz2strt = qz;
}

// Getting more information *************************************

// getting qz1 at minimum
double GetXMin (void) {
	return Xmin;
}

// getting qz2 at minimum
double GetYMin (void) {
	return Ymin;
}

// getting the status of ending of IterDeriv
// = 0 if converged
// = 1 if ended in looping
// = 2 if ended for too many steps
int GetStatus (void) {
	return status;
}

// getting the number of simple steps of IterDeriv
int GetNSteps (void) {
	return nstps;
}

// getting the number of steps using derivatives of IterDeriv
int GetNDerSteps (void) {
	return nderstps;
}

// getting the number of steps with rotations of IterDeriv
int GetNRotSteps (void) {
	return nrotstps;
}

// getting the number of steps with scanning of IterDeriv
int GetNScanSteps (void) {
	return nscnstps;
}

// getting the number of solutions of the 4th order equation
int GetNsol4 () {
	return nsol4;
}

// Getting the 3-momenta for a solution
int GetMomSol (int i, double pn1[], double pn2[]) {
	if (i > nsol4) return 0;
	pn1[0] = pn1x[i];
	pn1[1] = pn1y[i];
	pn1[2] = pn1z[i];
	pn2[0] = pn2x[i];
	pn2[1] = pn2y[i];
	pn2[2] = pn2z[i];
	return 1;
}

// printing the event contents
virtual void EvtPrint (double pb1[], double pb2[], double pl1[], double pl2[], double pn1[], double pn2[]);

 
private:

// private functions:

// Solves the constraint equations and returns a discriminant of type D1
virtual double SolveEqn1 (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double q1[], double q2[]);

// some kinematics functions
virtual double Minv2 (double p1[], double p2[]);
virtual int InterParabol (double dz1, double dz2, double D1, double D2, double D3, double dz[]);

// solution of a 4th order equation
int solveQuartic(double coef[], double sol[]);

// internal constants
int debug;
double mt, mW, mnu;
double xlo, xhi, ylo, yhi;
int nscan; 
double qz1strt, qz2strt;
double discconverge;
double Xmin, Ymin;
int status, nstps, nderstps, nrotstps, nscnstps;

int nsol4;
double pn1x[4], pn1y[4], pn1z[4];
double pn2x[4], pn2y[4], pn2z[4];




};

#endif


