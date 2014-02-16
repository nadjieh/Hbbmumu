// Compute a discriminant bewteen ttbar and stop direct production

#include "SolveTTbar.h"

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "TH2.h"


//****************************************************************************
// Constructor
//****************************************************************************

SolveTTbar::SolveTTbar () :
debug(0),
mt(173.), mW(80.), mnu(0.),
xlo(-3000.), xhi(3000.), ylo(-3000.), yhi(3000.), nscan(1000), qz1strt(0.), qz2strt(0.), discconverge(0.000001),
Xmin(0.), Ymin(0.),
status(-1), nstps(0), nderstps(0), nrotstps(0), nscnstps(0)
{

 return;
}

//****************************************************************************
// Destructor
//****************************************************************************

SolveTTbar::~SolveTTbar () {

}

	
double SolveTTbar::ScanEqn (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[], TH2D * hdiscr) {
// Scans the plane (qz1, qz2) for a minimum of the discriminant which is returned
//
// input: arrays pi[4] = px, py, pz, E of the visible particles (pbi for b's, pli for leptons)
//		  array  pTm[2] = px, py of the missing transverse momentum
//		  hdiscr = address of histogram if you want to plot the discriminant in the (qz1, qz2) plane

	double qz1, qz2;
	
	double pb1[4], pl1[4], pb2[4], pl2[4], pTm[2];
	for (int i = 0; i < 3; ++i) {
		pb1[i] = pb1in[i];
		pl1[i] = pl1in[i];
		pb2[i] = pb2in[i];
		pl2[i] = pl2in[i];
		if (i < 2) pTm[i] = pTmin[i];
	}
	pb1[3] = sqrt(pb1[0]*pb1[0] + pb1[1]*pb1[1] + pb1[2]*pb1[2]);
	pl1[3] = sqrt(pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]);
	pb2[3] = sqrt(pb2[0]*pb2[0] + pb2[1]*pb2[1] + pb2[2]*pb2[2]);
	pl2[3] = sqrt(pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
	
	double qz1low = xlo;
	double qz1hig = xhi;
	double qz2low = ylo;
	double qz2hig = yhi;
	double dz1 = (qz1hig - qz1low) / (float)nscan;
	double dz2 = (qz2hig - qz2low) / (float)nscan;
	double discmin = 99999.;
	double discminpre = discmin;
	double qz1min = 0., qz2min = 0.;

	double q1[4] = {0., 0., 0., 0.}, q2[4] = {0., 0., 0., 0.};
	
	int itry = 0;
	while(1) {
		itry++;
		for (int i = 0; i < nscan; ++i) {
			qz1 = qz1low + (float)i * dz1;
			q1[2] = qz1;
			for (int j = 0; j < nscan; ++j) {
				qz2 = qz2low + (float)j * dz2;
				q2[2] = qz2;
				double discr = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1, q2);
//				double discr = SolveEqn1 (pb1, pb2, pl1, pl2, pTm, q1, q2);
				if (hdiscr != NULL && itry <= 1) hdiscr->Fill(qz1+0.5*dz1, qz2+0.5*dz2, discr);
				if (discr < discmin) {
					discmin = discr;
					qz1min = qz1;
					qz2min = qz2;
				}
			}
//			cout << " discr = " << discr << endl;
		}
//		cout << " Minimum discr = " << discmin << " for qz1 = " << qz1min << ", qz2 = " << qz2min << endl;
		if (fabs(discmin-discminpre) < discconverge || itry > 100) break;
		discminpre = discmin;
		qz1low = qz1min - 2.*dz1;
		qz1hig = qz1min + 2.*dz1;
		dz1 = (qz1hig - qz1low) / (float)nscan;
		qz2low = qz2min - 2.*dz2;
		qz2hig = qz2min + 2.*dz2;
		dz2 = (qz2hig - qz2low) / (float)nscan;
//		cout << " Scan qz1min = " << qz1low << ", qz1max = " << qz1hig << ", qz2min = " << qz2low << ", qz2max = " << qz2hig
//		  << ", dz1 = " << dz1 << ", dz2 = " << dz2 << endl;
	} 
//	cout << " Minimum discr = " << discmin << " for qz1 = " << qz1min << ", qz2 = " << qz2min << endl;

	Xmin = qz1min;
	Ymin = qz2min;
	return discmin;
	
}

	
int SolveTTbar::ScanFndAll (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[], 
		double discsol[], double qz1sol[], double qz2sol[]) {
// Scans the plane (qz1, qz2) for all the minima of the discriminant 
//
// input: arrays pi[4] = px, py, pz, E of the visible particles (pbi for b's, pli for leptons)
//		  array  pTm[2] = px, py of the missing transverse momentum
//		  it returns the number of solutions and the values of the discriminator, qz1 and qz2 for each

	double qz1, qz2;
	
	double pb1[4], pl1[4], pb2[4], pl2[4], pTm[2];
	for (int i = 0; i < 3; ++i) {
		pb1[i] = pb1in[i];
		pl1[i] = pl1in[i];
		pb2[i] = pb2in[i];
		pl2[i] = pl2in[i];
		if (i < 2) pTm[i] = pTmin[i];
	}
	pb1[3] = sqrt(pb1[0]*pb1[0] + pb1[1]*pb1[1] + pb1[2]*pb1[2]);
	pl1[3] = sqrt(pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]);
	pb2[3] = sqrt(pb2[0]*pb2[0] + pb2[1]*pb2[1] + pb2[2]*pb2[2]);
	pl2[3] = sqrt(pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
	
	double qz1low = xlo;
	double qz1hig = xhi;
	double qz2low = ylo;
	double qz2hig = yhi;
	double dz1 = (qz1hig - qz1low) / (float)nscan;
	double dz2 = (qz2hig - qz2low) / (float)nscan;

	double q1[4] = {0., 0., 0., 0.}, q2[4] = {0., 0., 0., 0.};
	vector<double> disappr, qz1appr, qz2appr;
	int nsols = 0;
	
	// scan and find all candidate solutions
	for (int i = 0; i < nscan; ++i) {
		qz1 = qz1low + (float)i * dz1;
		q1[2] = qz1;
		for (int j = 0; j < nscan; ++j) {
			qz2 = qz2low + (float)j * dz2;
			q2[2] = qz2;
			double discr = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1, q2);
//			double discr = SolveEqn1 (pb1, pb2, pl1, pl2, pTm, q1, q2);
			if (nsols == 0) {
				nsols++;
				disappr.push_back(discr);
				qz1appr.push_back(qz1);
				qz2appr.push_back(qz2);
				continue;
			}
			double tol1 = 5. * dz1;
			double tol2 = 5. * dz2;
			int found = 0;
			for (int k = 0; k < nsols; ++k) {
				if (fabs(qz1-qz1appr[k]) < tol1 || fabs(qz2-qz2appr[k]) < tol2) {
					if (discr < disappr[k]) {
						disappr[k] = discr;
						qz1appr[k] = qz1;
						qz2appr[k] = qz2;
					}
					found = 1;
					break;
				}
			}
			if (!found) {
				nsols++;
				disappr.push_back(discr);
				qz1appr.push_back(qz1);
				qz2appr.push_back(qz2);
				if (nsols < 50) {
//					cout << " Sol = " << nsols << ", discr = " << discr << ", qz1 = " << qz1 << ", qz2 = " << qz2 << endl;
				}
			}
		}
//		cout << " discr = " << discr << endl;
	}
/*	cout << " Nber solutions = " << nsols << endl;
	int nmax = nsols;
	if (nsols > 50) nmax = 50;
	for (int k = 0; k < nmax; ++k) {
		cout << " Sol = " << k << ", discr = " << disappr[k] << ", qz1 = " << qz1appr[k] << ", qz2 = " << qz2appr[k] << endl;
	}
*/	
	// keep the 4 best solutions if too many
	int iflip = 0;
	do {
		iflip = 0;
		for (int i = 0; i < nsols; ++i) {
			for (int j = i+1; j < nsols; ++j) {
				if (disappr[j] < disappr[i]) {
					double discr = disappr[i];
					qz1 = qz1appr[i];
					qz2 = qz2appr[i];
					disappr[i] = disappr[j];
					qz1appr[i] = qz1appr[j];
					qz2appr[i] = qz2appr[j];
					disappr[j] = discr;
					qz1appr[j] = qz1;
					qz2appr[j] = qz2;
					iflip = 1;
				}
			}
		}
	} while(iflip);
	
//	cout << "Kept : " << endl;
	if (nsols > 4) nsols = 4;
	for (int k = 0; k < nsols; ++k) {
		discsol[k] = disappr[k];
		qz1sol[k] = qz1appr[k];
		qz2sol[k] = qz2appr[k];
//		cout << " Sol = " << k << ", discr = " << disappr[k] << ", qz1 = " << qz1appr[k] << ", qz2 = " << qz2appr[k] << endl;
	}
	
	
	// try improve the location of the minima
	for (int k = 0; k < nsols; ++k) {
		double dz1win = 2. * dz1;
		double dz2win = 2. * dz2;
		SetXWindow (qz1sol[k]-dz1win, qz1sol[k]+dz1win);
		SetYWindow (qz2sol[k]-dz2win, qz2sol[k]+dz2win);
		SetNscan (20);
		discsol[k] = ScanEqn (pb1, pb2, pl1, pl2, pTm, NULL);
		qz1sol[k] = GetXMin ();
		qz2sol[k] = GetYMin ();
	}

	return nsols;
}
		
	
double SolveTTbar::IterDeriv (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[]) {
// Iterate in the plane (qz1, qz2) for a minimum of the discriminant which is returned
// uses the analytical form of the derivatives
//
// input: arrays pi[4] = px, py, pz, E of the visible particles (pbi for b's, pli for leptons)
//		  array  pTm[2] = px, py of the missing transverse momentum

	double qz1, qz2;
	
	int useretry = 0;
	int usetryrot = 1;
	int usescan = 1;
	
	double dz1strt = 10.;
	double dz2strt = 10.;
	double dzmax = 100.;
//	double derDmodconv = 0.00001;
	double discmin = 99999.;
	double discpre = 99999.;
	double discConv = discconverge;
//	double discloop = 0.005;
	double discloop = 0.01*discconverge;
	int itrymax = nscan;
	
	double pb1[4], pl1[4], pb2[4], pl2[4], pTm[2];
	for (int i = 0; i < 3; ++i) {
		pb1[i] = pb1in[i];
		pl1[i] = pl1in[i];
		pb2[i] = pb2in[i];
		pl2[i] = pl2in[i];
		if (i < 2) pTm[i] = pTmin[i];
	}
	pb1[3] = sqrt(pb1[0]*pb1[0] + pb1[1]*pb1[1] + pb1[2]*pb1[2]);
	pl1[3] = sqrt(pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]);
	pb2[3] = sqrt(pb2[0]*pb2[0] + pb2[1]*pb2[1] + pb2[2]*pb2[2]);
	pl2[3] = sqrt(pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
	
	double q1[4] = {0., 0., 0., 0.}, q2[4] = {0., 0., 0., 0.};
	double fact1 = 0., fact2 = 0.;
	double derDqz1 = 0., derDqz2 = 0., derDmod = 0., derDqz1min = 0., derDqz2min = 0., derDmodmin = 0.;
	
	qz1 = qz1strt;
	qz2 = qz2strt;
	double dz1 = dz1strt;
	double dz2 = dz2strt;
	double qz1min = 0., qz2min = 0.;
//	double cut = 0.5;
	int itry = 0;
	double discr, derDisc[2];
	int retry = 0;
	int tryrot = 0;
	double dz1try = 0., dz2try = 0.;
	const int nangl = 50;
	double angl = 0., dangl = 6.283185307 / (double)nangl;
	double dzrotmin = 1.;
	double scanwindw = 1000.;
	double scanwindwmin = 1.;
	int ntryscan = 200;
	while(1) {
		itry++;
		if (debug >= 3) {
			cout << endl;
			cout << " itry = " << itry << endl;
		}
		nstps++;
		q1[2] = qz1;
		q2[2] = qz2;
		discr = DerivSolve (pb1, pb2, pl1, pl2, pTm, q1, q2, derDisc);
		derDqz1 = derDisc[0];
		derDqz2 = derDisc[1];
		
		if (fabs(discr) < discConv) {
			if (debug >= 1) {
				cout << " Converged, qz1 = " << qz1 << ", qz2 = " << qz2 << ", discr = " << discr << endl;
			}
			discmin = discr;
			qz1min = qz1;
			qz2min = qz2;
			status = 0;
			break;
		}
		if (fabs(discr-discpre) < discloop) {
			if (debug >= 1) {
				cout << " Looping, qz1 = " << qz1 << ", qz2 = " << qz2 << ", discr = " << discr << ", discpre = " << discpre << endl;
			}
			if (useretry && !retry) {
				retry = 1;
				dz1try = dz1strt;
				dz2try = dz2strt;
				discpre = 99999.;
			} else if (usetryrot && !tryrot) {
				tryrot = 1;
				dz1try = dz1strt;
				discpre = 99999.;
			} else {
				qz1min = qz1;
				qz2min = qz2;
				status = 1;
				break;
			}
		} else if (discr < discpre) {
			discpre = discr;
		}
		
//		derDmod = sqrt(derDqz1*derDqz1 + derDqz2*derDqz2);
		derDmod = 1. / fabs(derDqz1) + 1. / fabs(derDqz2);
		if (debug >= 3) {
//			cout << " itry = " << itry << " qz1 = " << qz1 << ", qz2 = " << qz2 << ", discr = " << discr << ", discmin = " << discmin << endl;
			cout << " discr = " << discr << ", discpre = " << discpre << ", discmin = " << discmin << endl;
		}
//		if (derDmod <= derDmodconv) {
//			break;
//		}
		// else, continue stepping
		if (discr < discmin && !retry) {
			discmin = discr;
			qz1min = qz1;
			qz2min = qz2;
			derDqz1min = derDqz1;
			derDqz2min = derDqz2;
			derDmodmin = derDmod;
			discpre = discmin;
			dz1 = -2. * fabs(dz1) * derDqz1/fabs(derDqz1);
			if ( fabs(dz1) > fabs(discr/derDqz1) ) dz1 = -discr/derDqz1;
			if ( fabs(dz1) > dzmax) dz1 = dzmax * dz1/fabs(dz1);
			dz2 = -2. * fabs(dz2) * derDqz2/fabs(derDqz2);
			if ( fabs(dz2) > fabs(discr/derDqz2) ) dz2 = -discr/derDqz2;
			if ( fabs(dz2) > dzmax) dz2 = dzmax * dz2/fabs(dz2);
			if (debug >= 2) {
				cout << " Accept, step dz1 = " << dz1 << ", dz2 = " << dz2 << endl;
			}
//			fact1 = -derDqz1 / derDmod;
//			fact2 = -derDqz2 / derDmod;
//			qz1 += fact1 * dz1;
//			qz2 += fact2 * dz2;
			qz1 += dz1;
			qz2 += dz2;
		} else if (!retry && !tryrot) {
			discr = discmin;
			derDqz1 = derDqz1min;
			derDqz2 = derDqz2min;
			dz1 = -0.5 * fabs(dz1) * derDqz1/fabs(derDqz1);
			if ( fabs(dz1) > fabs(discr/derDqz1) ) dz1 = -discr/derDqz1;
			dz2 = -0.5 * fabs(dz2) * derDqz2/fabs(derDqz2);
			if ( fabs(dz2) > fabs(discr/derDqz2) ) dz2 = -discr/derDqz2;
			if (debug >= 2) {
				cout << " Reject, step dz1 = " << dz1 << ", dz2 = " << dz2 << endl;
			}
			qz1 = qz1min + dz1;
			qz2 = qz2min + dz2;
		} else if (retry && !tryrot) {
			if (discr < discmin) {
				discmin = discr;
				qz1min = qz1;
				qz2min = qz2;
				derDqz1min = derDqz1;
				derDqz2min = derDqz2;
				derDmodmin = derDmod;
				if (debug >= 2) {
					cout << " Accept, step dz1 = " << dz1 << ", dz2 = " << dz2 << endl;
				}
				fact1 = - 1. / (derDqz1min * derDmodmin);
				fact2 = - 1. / (derDqz2min * derDmodmin);
				dz1try *= 2.;
				dz2try *= 2.;
			} else {
				// orient the step according to 1/derivative
				dz1try *= 0.5;
				dz2try *= 0.5;
				fact1 = - 1. / (derDqz1min * derDmodmin);
				fact2 = - 1. / (derDqz2min * derDmodmin);
				dz1 = fact1 * dz1try;
				dz2 = fact2 * dz2try;
				qz1 = qz1min + dz1;
				qz2 = qz2min + dz2;
				nderstps++;
				if (debug >= 2) {
					cout << " Retry, step dz1 = " << dz1 << ", dz2 = " << dz2 << " derDmod = " << derDmod << endl;
				}
			} 
		} else if (tryrot) {
			// rotate the direction of the step
			double qz1best = 0., qz2best = 0., dz1best = 0., dz2best = 0.;
			dz1try = sqrt(dz1*dz1 + dz2*dz2);
			if (dz1try < dzrotmin) dz1try = dzrotmin;
			for (int i = 0; i < nangl; ++i) {
				angl += dangl;
				dz1 = dz1try * cos (angl);
				dz2 = dz1try * sin (angl);
				qz1 = qz1min + dz1;
				qz2 = qz2min + dz2;
				q1[2] = qz1;
				q2[2] = qz2;
				discr = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1, q2);
				if (discr < discpre) {
					discpre = discr;
					qz1best = qz1;
					qz2best = qz2;
					dz1best = dz1;
					dz2best = dz2;
				}
			}
			nrotstps++;
			if (discpre < discmin) {
				qz1 = qz1best;
				qz2 = qz2best;
				if (debug >= 2) {
					cout << " Rotation improvement found " << discpre << ", qz1 = " << qz1 << ", qz2 = " << qz2 << ", dz = " << dz1try << endl;
				}
				discpre = discmin;
				dz1 = dz1best;
				dz2 = dz2best;
				retry = 0;
			} else if (usescan) {
				// if nothing worked, try a scan
				if (debug >= 2) {
					cout << " No rotation improvement, try a scan around qz1 = " << qz1min << ", qz2 = " << qz2min << endl;
				}
				if (scanwindw > scanwindwmin) {
					SetXWindow (qz1min-scanwindw, qz1min+scanwindw);
					SetYWindow (qz2min-scanwindw, qz2min+scanwindw);
					SetNscan (ntryscan);
					SetDiscrConv (0.000001);
					discr = ScanEqn (pb1, pb2, pl1, pl2, pTm, NULL);
					if (discr < discmin) {
						qz1 = GetXMin ();
						qz2 = GetYMin ();
						if (debug >= 2) {
							cout << " Scan Improved, discr = " << discr << ", qz1 = " << qz1 << ", qz2 = " << qz2 << ", windw = " << scanwindw << endl;
						}
						scanwindw *= 0.5;
//						ntryscan = 0.75 * (float)ntryscan;
					} else {
						if (debug >= 2) {
							cout << " No scan Improvement, discr = " << discr << endl;
						}
						discpre = discmin;
					}
					usescan = 0;
				} else {
					qz1 = qz1min;
					qz2 = qz2min;
				}
				nscnstps++;
			} else {
				qz1 = qz1min;
				qz2 = qz2min;
			}
		}
//		cout << " Minimum discr = " << discmin << " for qz1 = " << qz1min << ", qz2 = " << qz2min << endl;
//		if (fabs(discmin-discminpre) < discconverge) break;
//		discminpre = discmin;
		if (itry > itrymax) {
			if (debug >= 1) {
				cout << " Reached max number of iterations = " << itry << endl;
			}
			qz1min = qz1;
			qz2min = qz2;
			status = 2;
			break;
		}
	} 
//	cout << " Minimum discr = " << discmin << " for qz1 = " << qz1min << ", qz2 = " << qz2min << endl;

	Xmin = qz1min;
	Ymin = qz2min;
	return discmin;
	
}
	
double SolveTTbar::IterEqn (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[]) {
// Iterate in the plane (qz1, qz2) for a minimum of the discriminant which is returned
//
// input: arrays pi[4] = px, py, pz, E of the visible particles (pbi for b's, pli for leptons)
//		  array  pTm[2] = px, py of the missing transverse momentum

// *************** OBSOLETE, NOT TO BE USED *******************

	
	double pb1[4], pl1[4], pb2[4], pl2[4], pTm[2];
	for (int i = 0; i < 3; ++i) {
		pb1[i] = pb1in[i];
		pl1[i] = pl1in[i];
		pb2[i] = pb2in[i];
		pl2[i] = pl2in[i];
		if (i < 2) pTm[i] = pTmin[i];
	}
	pb1[3] = sqrt(pb1[0]*pb1[0] + pb1[1]*pb1[1] + pb1[2]*pb1[2]);
	pl1[3] = sqrt(pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]);
	pb2[3] = sqrt(pb2[0]*pb2[0] + pb2[1]*pb2[1] + pb2[2]*pb2[2]);
	pl2[3] = sqrt(pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
	
	double dz1init = 10.;
	double dz2init = 10.;
	double dzmin = 0.01;
	double disaccept = 1.;
	int itrymax = 1000;
	
	double qz1, qz2;
	qz1 = qz1strt;
	qz2 = qz2strt;
	
	double q1[4][4], q2[4][4];
	
	double qz1min, qz2min;
	qz1min = 0.;
	qz2min = 0.;
	q1[0][2] = qz1min;
	q2[0][2] = qz2min;
	double discmin = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1[0], q2[0]);
	double discr[4];
	double qz1fir = -1., qz1sec = -1., qz2fir = -1., qz2sec = -1.;
	int idirfir = -1, idirsec = -1;
	double discfir = 99999., discsec = 99999.;
	double dz1 = dz1init;
	double dz2 = dz2init;
	double cut = 0.5;
	int trydiag = 0;
	int retry = 0;
	
	int itry = 0;
	while(1) {
		itry++;
		discfir = 99999.;
		discsec = 99999.;
		trydiag = 1;
		q1[0][2] = qz1 + dz1;
		q2[0][2] = qz2;
		discr[0] = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1[0], q2[0]);
		q1[1][2] = qz1;
		q2[1][2] = qz2 + dz2;
		discr[1] = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1[1], q2[1]);
		q1[2][2] = qz1 - dz1;
		q2[2][2] = qz2;
		discr[2] = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1[2], q2[2]);
		q1[3][2] = qz1;
		q2[3][2] = qz2 - dz2;
		discr[3] = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1[3], q2[3]);
		for (int i = 0; i < 4; ++i) {
			if (discr[i] < discfir) {
				discfir = discr[i];
				qz1fir = q1[i][2];
				qz2fir = q2[i][2];
				idirfir = i;
			}
		}
		for (int i = 0; i < 4; ++i) {
			if (i%2 == idirfir%2) continue;
			if (discr[i] < discsec) {
				discsec = discr[i];
				qz1sec = q1[i][2];
				qz2sec = q2[i][2];
				idirsec = i;
			}
		}
		if (debug >= 2) {
			cout << endl;
			cout << " itry = " << itry << ", qz1 = " << qz1 << ", qz2 = " << qz2 << endl;
		}
		if (debug >= 3) {
			cout << " dirfir = " << idirfir << ", discfir = " << discfir
				<< ", dirsec = " << idirsec << ", discsec = " << discsec << ", discmin = " << discmin << endl;
		}
		
		if (discfir <= discmin && discsec <= discmin) {
			discmin = discfir;
			qz1min = qz1;
			qz2min = qz2;
			double discsum = 1./discfir + 1./discsec;
			qz1 = (qz1fir/discfir + qz1sec/discsec) / discsum;
			qz2 = (qz2fir/discfir + qz2sec/discsec) / discsum;
		} else if (discfir <= discmin && discsec > discmin) {
			discmin = discfir;
			qz1min = qz1;
			qz2min = qz2;
			qz1 = qz1fir;
			qz2 = qz2fir;
		} else {
			if (trydiag) {
				// on diagonal, point in the middle between dz1 and dz2
				double q1dia[4][4], q2dia[4][4], discrdia[4];
				itry++;
				trydiag = 0;
				q1dia[0][2] = qz1 + 0.5 * dz1;
				q2dia[0][2] = qz2 + 0.5 * dz2;
				discrdia[0] = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1dia[0], q2dia[0]);
				q1dia[1][2] = qz1 - 0.5 * dz1;
				q2dia[1][2] = qz2 + 0.5 * dz2;
				discrdia[1] = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1dia[1], q2dia[1]);
				q1dia[2][2] = qz1 - 0.5 * dz1;
				q2dia[2][2] = qz2 - 0.5 * dz2;
				discrdia[2] = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1dia[2], q2dia[2]);
				q1dia[3][2] = qz1 + 0.5 * dz1;
				q2dia[3][2] = qz2 - 0.5 * dz2;
				discrdia[3] = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1dia[3], q2dia[3]);
				discfir = 99999.;
				discsec = 99999.;
				for (int i = 0; i < 4; ++i) {
					if (discrdia[i] < discfir) {
						discfir = discrdia[i];
						qz1fir = q1dia[i][2];
						qz2fir = q2dia[i][2];
						idirfir = i;
					}
				}
				for (int i = 0; i < 4; ++i) {
					if (i%2 == idirfir%2) continue;
					if (discrdia[i] < discsec) {
						discsec = discrdia[i];
						qz1sec = q1dia[i][2];
						qz2sec = q2dia[i][2];
						idirsec = i;
					}
				}
				if (debug >= 2) {
					cout << endl;
					cout << " itry = " << itry << ", qz1 = " << qz1 << ", qz2 = " << qz2 << endl;
				}
				if (debug >= 3) {
					cout << " Along the diagonal " << endl;
					cout << " diafir = " << idirfir << ", discfir = " << discfir
						<< ", diasec = " << idirsec << ", discsec = " << discsec << ", discmin = " << discmin << endl;
				}
				if (discfir <= discmin && discsec <= discmin) {
					discmin = discfir;
					qz1min = qz1fir;
					qz2min = qz2fir;
					double discsum = 1./discfir + 1./discsec;
					qz1 = (qz1fir/discfir + qz1sec/discsec) / discsum;
					qz2 = (qz2fir/discfir + qz2sec/discsec) / discsum;
				} else if (discfir <= discmin && discsec > discmin) {
					discmin = discfir;
					qz1min = qz1fir;
					qz2min = qz2fir;
					qz1 = qz1fir;
					qz2 = qz2fir;
				} else {
					double q1int[4], q2int[4];
					double dz[2] = {0., 0.};
					discfir = 99999.;
					for (int i = 0; i < 4; ++i) {
						int j = i + 1;
						if (i == 3) j = 0;
						int flag = InterParabol (dz1, dz2, discr[i], discrdia[i], discr[j], dz);
						if (flag != 0) {
							q1int[2] = qz1 + dz[0];
							q2int[2] = qz2 + dz[1];
							double discint = SolveEqn_qz12 (pb1, pb2, pl1, pl2, pTm, q1int, q2int);
							if (debug >= 3) {
								cout << " Interpolated dir = " << i << ", discr = " << discint << endl;
							}
							if (discint < discfir) {
								discfir = discint;
								qz1fir = q1int[2];
								qz2fir = q2int[2];
								idirfir = i;
							}
						}
					}
					if (discfir < discmin) {
						discmin = discfir;
						qz1min = qz1fir;
						qz2min = qz2fir;
						qz1 = qz1fir;
						qz2 = qz2fir;
						continue;
					}
					if (retry && dz1 <= 0.1*dz1init && dz2 <= 0.1*dz2init) {
						if (debug >= 3) {
							cout << " Give up here ********** " << endl;
						}
						break;
					}
					if (!retry) {
						dz1 = cut * dz1;
						dz2 = cut * dz2;
					} else {
						dz1 -= 10.;
						dz2 -= 10.;
					}
					if (debug >= 3) {
						cout << " reduce step to dz = " << dz1 << endl;
					}
					if (dz1 < dzmin && dz2 < dzmin) {
						if (debug >= 1) {
							cout << " Converged !!!!!!!!!! nsteps = " << itry << ", discmin = " << discmin << endl;
						}
						if (!retry && discmin > disaccept) {
							retry = 1;
							qz1 = qz1min;
							qz2 = qz2min;
							dz1 = 50. * dz1init;
							dz2 = 50. * dz2init;
							if (debug >= 1) {
								cout << " Retry with 50x dz = " << dz1 << endl;
							}
						} else break;
					}
				}
			} else {
				if (retry && dz1 <= 0.1*dz1init && dz2 <= 0.1*dz2init) {
					if (debug >= 3) {
						cout << " Give up here ********** " << endl;
					}
					break;
				}
				dz1 = cut * dz1;
				dz2 = cut * dz2;
				if (debug >= 3) {
					cout << " reduce step to dz = " << dz1 << endl;
				}
			
				if (dz1 < dzmin && dz2 < dzmin) {
					if (debug >= 1) {
						cout << " Converged !!!!!!!!!! nsteps = " << itry << ", discmin = " << discmin << endl;
					}
					if (!retry && discmin > disaccept) {
						retry = 1;
						qz1 = qz1min;
						qz2 = qz2min;
						dz1 = 2. * dz1init;
						dz2 = 2. * dz2init;
						if (debug >= 1) {
							cout << " Retry with double dz = " << dz1 << endl;
						}
					} else break;
				}
			}
		}
		
		if (itry > itrymax) {
			if (debug >= 1) {
				cout << " Reached max number of iterations = " << itry << endl;
			}
			break;
		}
	} 
//	cout << " Minimum discr = " << discmin << " for qz1 = " << qz1min << ", qz2 = " << qz2min << endl;

	Xmin = qz1min;
	Ymin = qz2min;
	return discmin;
	
}


double SolveTTbar::SolveEqn_qz12 (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double q1[], double q2[]) {
// Solves the constraint equations and returns a discriminant of type D2
// expects qz in q[2] to be filled for q1 and q2
// also completes the arrays q1 and q2
	
	
	double pb1mod = pb1[3];
	double pl1mod = pl1[3];
	double pb2mod = pb2[3];
	double pl2mod = pl2[3];
	double qz1 = q1[2];
	double qz2 = q2[2];

	// compute the constants
	double A1 = 0.5 * (mt*mt/pb1mod - mW*mW/pl1mod);
	double A2 = 0.5 * (mt*mt/pb2mod - mW*mW/pl2mod);
	double v1[3], v2[3];
	v1[0] = pb1[0] / pb1mod - pl1[0] / pl1mod;
	v1[1] = pb1[1] / pb1mod - pl1[1] / pl1mod;
	v1[2] = pb1[2] / pb1mod - pl1[2] / pl1mod;
	v2[0] = pb2[0] / pb2mod - pl2[0] / pl2mod;
	v2[1] = pb2[1] / pb2mod - pl2[1] / pl2mod;
	v2[2] = pb2[2] / pb2mod - pl2[2] / pl2mod;
	double Ml1b1 = Minv2 (pl1, pb1);
	double Ml2b2 = Minv2 (pl2, pb2);
	double c1 = (Ml1b1*Ml1b1 + mW*mW) / (2.*pb1mod) - mnu*mnu/(2.*pl1mod) - A1;
	double c2 = (Ml2b2*Ml2b2 + mW*mW) / (2.*pb2mod) - mnu*mnu/(2.*pl2mod) - A2 - v2[0]*pTm[0] - v2[1]*pTm[1];
	
	// Now compute qx1 and qy1
	double num1 = v1[2]*qz1 - c1;
	double num2 = v2[2]*qz2 - c2;
	double denom = -v1[0]*v2[1] + v1[1]*v2[0];
	q1[0] =  (v2[1]*num1 + v1[1]*num2) / denom;
	q1[1] = -(v2[0]*num1 + v1[0]*num2) / denom;
	q1[2] = qz1;
	q2[0] = pTm[0] - q1[0]; 
	q2[1] = pTm[1] - q1[1];
	q2[2] = qz2;
	if (debug >= 4) {
		cout << " q1 px, py, pz = " << q1[0] << ", " << q1[1] << ", " << q1[2] << endl;
		cout << " q2 px, py, pz = " << q2[0] << ", " << q2[1] << ", " << q2[2] << endl;
	}
	
	double Eq1 = sqrt(q1[0]*q1[0] + q1[1]*q1[1] + q1[2]*q1[2] + mnu*mnu);
	double Eq2 = sqrt(q2[0]*q2[0] + q2[1]*q2[1] + q2[2]*q2[2] + mnu*mnu);
	q1[3] = Eq1;
	q2[3] = Eq2;
	
	double mW1 = sqrt(2.*(pl1mod*Eq1 - pl1[0]*q1[0] - pl1[1]*q1[1] - pl1[2]*q1[2]) + mnu*mnu);
	double mW2 = sqrt(2.*(pl2mod*Eq2 - pl2[0]*q2[0] - pl2[1]*q2[1] - pl2[2]*q2[2]) + mnu*mnu);
//	double dW = ( (mW1-mW)*(mW1-mW) + (mW2-mW)*(mW2-mW) ) / (mW*mW);
	double dW = (mW1-mW)*(mW1-mW) + (mW2-mW)*(mW2-mW);
	double mt1 = sqrt(2.*(pb1mod*Eq1 - pb1[0]*q1[0] - pb1[1]*q1[1] - pb1[2]*q1[2]) + Ml1b1*Ml1b1 + mW*mW);
	double mt2 = sqrt(2.*(pb2mod*Eq2 - pb2[0]*q2[0] - pb2[1]*q2[1] - pb2[2]*q2[2]) + Ml2b2*Ml2b2 + mW*mW);
//	double dt = ( (mt1-mt)*(mt1-mt) + (mt2-mt)*(mt2-mt) ) / (mt*mt);
	double dt = (mt1-mt)*(mt1-mt) + (mt2-mt)*(mt2-mt);
	double discr = 0.5 * sqrt(dW + dt);
//	cout << " mW1 = " << mW1 << ", mW2 = " << mW2 << ", Disc = " << dW << endl;
//	cout << " mt1 = " << mt1 << ", mt2 = " << mt2 << ", Disc = " << dt << endl;

	return discr;
	
}


double SolveTTbar::SolveEqn1 (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double q1[], double q2[]) {
// Solves the constraint equations and returns a discriminant of type D1
// expects qz in q[2] to be filled for q1 and q2
// also completes the arrays q1 and q2

// *************** OBSOLETE, NOT TO BE USED *******************
	
	
	double pb1mod = pb1[3];
	double pl1mod = pl1[3];
	double pb2mod = pb2[3];
	double pl2mod = pl2[3];
	double qz1 = q1[2];
	double qz2 = q2[2];

	// compute the constants
	double A1 = 0.5 * (mt*mt/pb1mod - mW*mW/pl1mod);
	double A2 = 0.5 * (mt*mt/pb2mod - mW*mW/pl2mod);
	double v1[3], v2[3];
	v1[0] = pb1[0] / pb1mod - pl1[0] / pl1mod;
	v1[1] = pb1[1] / pb1mod - pl1[1] / pl1mod;
	v1[2] = pb1[2] / pb1mod - pl1[2] / pl1mod;
	v2[0] = pb2[0] / pb2mod - pl2[0] / pl2mod;
	v2[1] = pb2[1] / pb2mod - pl2[1] / pl2mod;
	v2[2] = pb2[2] / pb2mod - pl2[2] / pl2mod;
	double Ml1b1 = Minv2 (pl1, pb1);
	double Ml2b2 = Minv2 (pl2, pb2);
	double c1 = (Ml1b1*Ml1b1 + mW*mW) / (2.*pb1mod) - A1;
	double c2 = (Ml2b2*Ml2b2 + mW*mW) / (2.*pb2mod) - A2 - v2[0]*pTm[0] - v2[1]*pTm[1];
	
	// Now compute qx1 and qy1
	double num1 = v1[2]*qz1 - c1;
	double num2 = v2[2]*qz2 - c2;
	double denom = -v1[0]*v2[1] + v1[1]*v2[0];
	q1[0] =  (v2[1]*num1 + v1[1]*num2) / denom;
	q1[1] = -(v2[0]*num1 + v1[0]*num2) / denom;
	q1[2] = qz1;
	q2[0] = pTm[0] - q1[0];
	q2[1] = pTm[1] - q1[1];
	q2[2] = qz2;
//	cout << " q1 px, py, pz = " << q1[0] << ", " << q1[1] << ", " << q1[2] << endl;
//	cout << " q2 px, py, pz = " << q2[0] << ", " << q2[1] << ", " << q2[2] << endl;
	
	double q1mod = sqrt(q1[0]*q1[0] + q1[1]*q1[1] + q1[2]*q1[2]);
	double q2mod = sqrt(q2[0]*q2[0] + q2[1]*q2[1] + q2[2]*q2[2]);
	q1[3] = q1mod;
	q2[3] = q2mod;
	
	double mW1 = sqrt(2.*(pl1mod*q1mod - pl1[0]*q1[0] - pl1[1]*q1[1] - pl1[2]*q1[2]));
	double mW2 = sqrt(2.*(pl2mod*q2mod - pl2[0]*q2[0] - pl2[1]*q2[1] - pl2[2]*q2[2]));
	double dW = fabs(mW1-mW) + fabs(mW2-mW);
	double mt1 = sqrt(2.*(pb1mod*q1mod - pb1[0]*q1[0] - pb1[1]*q1[1] - pb1[2]*q1[2]) + Ml1b1*Ml1b1 + mW*mW);
	double mt2 = sqrt(2.*(pb2mod*q2mod - pb2[0]*q2[0] - pb2[1]*q2[1] - pb2[2]*q2[2]) + Ml2b2*Ml2b2 + mW*mW);
	double dt = fabs(mt1-mt) + fabs(mt2-mt);
	double discr = 0.25 * (dW + dt);
//	cout << " mW1 = " << mW1 << ", mW2 = " << mW2 << ", Disc = " << dW << endl;
//	cout << " mt1 = " << mt1 << ", mt2 = " << mt2 << ", Disc = " << dt << endl;

	return discr;
	
}


double SolveTTbar::DerivSolve (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double q1[], double q2[], double derDisc[]) {
// Solves the constraint equations and returns a discriminant and its derivatives
// expects qz in q[2] to be filled for q1 and q2 and completes the arrays q1 and q2
// derivatives of the discrimunant are returned in derDisc[0] = dD/dqz1, [1] = dD/dqz2
	
	double pb1mod = pb1[3];
	double pl1mod = pl1[3];
	double pb2mod = pb2[3];
	double pl2mod = pl2[3];
	double qz1 = q1[2];
	double qz2 = q2[2];

	// compute the constants
	double A1 = 0.5 * (mt*mt/pb1mod - mW*mW/pl1mod);
	double A2 = 0.5 * (mt*mt/pb2mod - mW*mW/pl2mod);
	double v1[3], v2[3];
	v1[0] = pb1[0] / pb1[3] - pl1[0] / pl1[3];
	v1[1] = pb1[1] / pb1[3] - pl1[1] / pl1[3];
	v1[2] = pb1[2] / pb1[3] - pl1[2] / pl1[3];
	v2[0] = pb2[0] / pb2[3] - pl2[0] / pl2[3];
	v2[1] = pb2[1] / pb2[3] - pl2[1] / pl2[3];
	v2[2] = pb2[2] / pb2[3] - pl2[2] / pl2[3];
	double Ml1b1 = Minv2 (pl1, pb1);
	double Ml2b2 = Minv2 (pl2, pb2);
	double c1 = (Ml1b1*Ml1b1 + mW*mW) / (2.*pb1mod) - mnu*mnu/(2.*pl1mod) - A1;
	double c2 = (Ml2b2*Ml2b2 + mW*mW) / (2.*pb2mod) - mnu*mnu/(2.*pl2mod) - A2 - v2[0]*pTm[0] - v2[1]*pTm[1];
	
//	double fact1 = 0., fact2 = 0.;
	double derDqz1 = 0., derDqz2 = 0.;
	
//	double discr;

//	discr = SolveEqn (pb1, pb2, pl1, pl2, pTm, q1, q2);
	
	// Now compute qx1 and qy1
	double num1 = v1[2]*qz1 - c1;
	double num2 = v2[2]*qz2 - c2;
	double denom = -v1[0]*v2[1] + v1[1]*v2[0];
	q1[0] =  (v2[1]*num1 + v1[1]*num2) / denom;
	q1[1] = -(v2[0]*num1 + v1[0]*num2) / denom;
	q1[2] = qz1;
	q2[0] = pTm[0] - q1[0];
	q2[1] = pTm[1] - q1[1];
	q2[2] = qz2;
//	cout << " q1 px, py, pz = " << q1[0] << ", " << q1[1] << ", " << q1[2] << endl;
//	cout << " q2 px, py, pz = " << q2[0] << ", " << q2[1] << ", " << q2[2] << endl;
	
	double Eq1 = sqrt(q1[0]*q1[0] + q1[1]*q1[1] + q1[2]*q1[2] + mnu*mnu);
	double Eq2 = sqrt(q2[0]*q2[0] + q2[1]*q2[1] + q2[2]*q2[2] + mnu*mnu);
	q1[3] = Eq1;
	q2[3] = Eq2;
	
	double mW1 = sqrt(2.*(pl1mod*Eq1 - pl1[0]*q1[0] - pl1[1]*q1[1] - pl1[2]*q1[2]) + mnu*mnu);
	double mW2 = sqrt(2.*(pl2mod*Eq2 - pl2[0]*q2[0] - pl2[1]*q2[1] - pl2[2]*q2[2]) + mnu*mnu);
	double mt1 = sqrt(2.*(pb1mod*Eq1 - pb1[0]*q1[0] - pb1[1]*q1[1] - pb1[2]*q1[2]) + Ml1b1*Ml1b1 + mW*mW);
	double mt2 = sqrt(2.*(pb2mod*Eq2 - pb2[0]*q2[0] - pb2[1]*q2[1] - pb2[2]*q2[2]) + Ml2b2*Ml2b2 + mW*mW);
	double derqx1qz1 =  v2[1] * v1[2] / denom;
	double derqy1qz1 = -v2[0] * v1[2] / denom;
	double derqT1dz1 = q1[0]*derqx1qz1 + q1[1]*derqy1qz1;
	double derqT2dz1 = -(pTm[0]-q1[0])*derqx1qz1 - (pTm[1]-q1[1])*derqy1qz1;
	double dermW1sqqz1 = 2. * ( pl1[3]*(qz1+derqT1dz1)/q1[3] - (pl1[2]+pl1[0]*derqx1qz1+pl1[1]*derqy1qz1) );
	double dermW1qz1 = dermW1sqqz1 / (2.*mW1);
	double dermt1sqqz1 = 2. * ( pb1[3]*(qz1+derqT1dz1)/q1[3] - (pb1[2]+pb1[0]*derqx1qz1+pb1[1]*derqy1qz1) ) + dermW1sqqz1;
	double dermt1qz1 = dermt1sqqz1 / (2.*mt1);
	double dermW2sqqz1 = 2. * ( pl2[3]*derqT2dz1/q2[3] + (pl2[0]*derqx1qz1+pl2[1]*derqy1qz1) );
	double dermW2qz1 = dermW2sqqz1 / (2.*mW2);
	double dermt2sqqz1 = 2. * ( pb2[3]*derqT2dz1/q2[3] + (pb2[0]*derqx1qz1+pb2[1]*derqy1qz1) ) + dermW2sqqz1;
	double dermt2qz1 = dermt2sqqz1 / (2.*mt2);
//		if (mW1 > mW) dermW1qz1 = -dermW1qz1;
//		if (mt1 > mt) dermt1qz1 = -dermt1qz1;
//		if (mW2 > mW) dermW2qz1 = -dermW2qz1;
//		if (mt2 > mt) dermt2qz1 = -dermt2qz1;
//		double derDqz1 = dermW1qz1 + dermW2qz1 + dermt1qz1 + dermt2qz1;
//		double D = 0.25 * sqrt( ((mW1-mW)*(mW1-mW) + (mW2-mW)*(mW2-mW))/(mW*mW) + ((mt1-mt)*(mt1-mt) + (mt2-mt)*(mt2-mt))/(mt*mt) );
//		double derDqz1 = 0.0625 * (((mW1-mW)*dermW1qz1 + (mW2-mW)*dermW2qz1)/(mW*mW) + ((mt1-mt)*dermt1qz1 + (mt2-mt)*dermt2qz1)/(mt*mt)) / D;
	double D = 0.5 * sqrt( ((mW1-mW)*(mW1-mW) + (mW2-mW)*(mW2-mW)) + ((mt1-mt)*(mt1-mt) + (mt2-mt)*(mt2-mt)) );
	derDqz1 = 0.25 * ( (mW1-mW)*dermW1qz1 + (mW2-mW)*dermW2qz1 + (mt1-mt)*dermt1qz1 + (mt2-mt)*dermt2qz1 ) / D;
	double derqx2qz2 = -v1[1] * v2[2] / denom;
	double derqy2qz2 =  v1[0] * v2[2] / denom;
	double derqT2dz2 = q2[0]*derqx2qz2 + q2[1]*derqy2qz2;
	double derqT1dz2 = -(pTm[0]-q2[0])*derqx2qz2 - (pTm[1]-q2[1])*derqy2qz2;
	double dermW2sqqz2 = 2. * ( pl2[3]*(qz2+derqT2dz2)/q2[3] - (pl2[2]+pl2[0]*derqx2qz2+pl2[1]*derqy2qz2) );
	double dermW2qz2 = dermW2sqqz2 / (2.*mW2);
	double dermt2sqqz2 = 2. * ( pb2[3]*(qz2+derqT2dz2)/q2[3] - (pb2[2]+pb2[0]*derqx2qz2+pb2[1]*derqy2qz2) ) + dermW2sqqz2;
	double dermt2qz2 = dermt2sqqz2 / (2.*mt2);
	double dermW1sqqz2 = 2. * ( pl1[3]*derqT1dz2/q1[3] + (pl1[0]*derqx2qz2+pl1[1]*derqy2qz2) );
	double dermW1qz2 = dermW1sqqz2 / (2.*mW1);
	double dermt1sqqz2 = 2. * ( pb1[3]*derqT1dz2/q1[3] + (pb1[0]*derqx2qz2+pb1[1]*derqy2qz2) + dermW1sqqz2);
	double dermt1qz2 = dermt1sqqz2 / (2.*mt1);
//		if (mW1 > mW) dermW1qz2 = -dermW1qz2;
//		if (mt1 > mt) dermt1qz2 = -dermt1qz2;
//		if (mW2 > mW) dermW2qz2 = -dermW2qz2;
//		if (mt2 > mt) dermt2qz2 = -dermt2qz2;
//		double derDqz2 = dermW2qz2 + dermW1qz2 + dermt2qz2 + dermt1qz2;
//		double derDqz2 = 0.25 * (((mW2-mW)*dermW2qz2 + (mW1-mW)*dermW1qz2)/(mW*mW) + ((mt2-mt)*dermt2qz2 + (mt1-mt)*dermt1qz2)/(mt*mt)) / D;
	derDqz2 = 0.25 * (((mW2-mW)*dermW2qz2 + (mW1-mW)*dermW1qz2) + ((mt2-mt)*dermt2qz2 + (mt1-mt)*dermt1qz2)) / D;
	if (debug >= 3) {
		cout << " qz1 = " << qz1 << ", qz2 = " << qz2 << endl;
		cout << " mW1 = " << mW1 << ", mW2 = " << mW2 << ", mt1 = " << mt1 << ", mt2 = " << mt2 << endl;
		cout << " dermW1qz1 = " << dermW1qz1 << ", dermW2qz1 = " << dermW2qz1 << ", dermt1qz1 = " << dermt1qz1 << ", dermt2qz1 = " << dermt2qz1 << endl;
		cout << " dermW1qz2 = " << dermW1qz2 << ", dermW2qz2 = " << dermW2qz2 << ", dermt1qz2 = " << dermt1qz2 << ", dermt2qz2 = " << dermt2qz2 << endl;
		cout << " D = " << D << ", derDqz1 = " << derDqz1 << ", derDqz2 = " << derDqz2 << endl;
	}
	
	derDisc[0] = derDqz1;
	derDisc[1] = derDqz2;
	return D;
}

	
double SolveTTbar::ScanEqn_qxy (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[], TH2D * hdiscrxy) {
// Scans the plane (qx1, qy1) for a minimum of the discriminant which is returned
//
// input: arrays pi[4] = px, py, pz, E of the visible particles (pbi for b's, pli for leptons)
//		  array  pTm[2] = px, py of the missing transverse momentum
//		  hdiscrxy = address of histogram if you want to plot the discriminant in the (qx1, qy1) plane

	double pb1[4], pl1[4], pb2[4], pl2[4], pTm[2];
	for (int i = 0; i < 3; ++i) {
		pb1[i] = pb1in[i];
		pl1[i] = pl1in[i];
		pb2[i] = pb2in[i];
		pl2[i] = pl2in[i];
		if (i < 2) pTm[i] = pTmin[i];
	}
	pb1[3] = sqrt(pb1[0]*pb1[0] + pb1[1]*pb1[1] + pb1[2]*pb1[2]);
	pl1[3] = sqrt(pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]);
	pb2[3] = sqrt(pb2[0]*pb2[0] + pb2[1]*pb2[1] + pb2[2]*pb2[2]);
	pl2[3] = sqrt(pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
	
	double qx1, qy1;
	
	double qx1low = xlo;
	double qx1hig = xhi;
	double dx1 = (qx1hig - qx1low) / (float)nscan;
	double dy1 = dx1;
	double discmin = 99999.;
	double discminpre = discmin;
	double qx1min = 0., qy1min = 0.;

	double q1[4] = {0., 0., 0., 0.}, q2[4] = {0., 0., 0., 0.};
	
	int itry = 0;
	while(1) {
		itry++;
		nstps++;
		for (int i = 0; i < nscan; ++i) {
			qx1 = qx1low + (float)i * dx1;
			q1[0] = qx1;
			double qy[2];
			int ifl = qyfrqx (pl1, 0., mW, qx1, qy);
			if (!ifl) continue;
			double qy1low = qy[0];
			double qy1hig = qy[1];
			int nscy = (qy1hig - qy1low) / dy1 + 1.;
			for (int j = 0; j < nscy; ++j) {
				qy1 = qy1low + (float)j * dy1;
				q1[1] = qy1;
				double discr = SolveEqn_qxy1 (pb1, pb2, pl1, pl2, pTm, q1, q2);
				if (debug >= 4) {
					cout << "  qx1 = " << qx1 << " qy1 = " << qy1 << " discr = " << discr << endl;
				}
				if (hdiscrxy != NULL && itry <= 1) hdiscrxy->Fill(qx1+0.5*dx1, qy1+0.5*dy1, discr);
				if (discr < discmin) {
					discmin = discr;
					qx1min = qx1;
					qy1min = qy1;
				}
			}
			if (debug >= 3) {
				cout << " itry = " << itry << " ibin = " << i << " qx1 = " << qx1 << " min discr = " << discmin << endl;
			}
		}
//		cout << " Minimum discr = " << discmin << " for qx1 = " << qx1min << ", qy1 = " << qy1min << endl;
//		if (fabs(discmin-discminpre) < discconverge || hdiscrxy != NULL) break;
		if (fabs(discmin-discminpre) < discconverge || itry > 100) break;
//		break;
		if (debug >= 2) {
			cout << " itry = " << itry << " qx1 = " << qx1min << " min discr = " << discmin << " discminpre = " << discminpre << endl;
		}
		nscnstps++;
		discminpre = discmin;
		qx1low = qx1min - 2.*dx1;
		qx1hig = qx1min + 2.*dx1;
		dx1 = (qx1hig - qx1low) / (float)nscan;
//		cout << " Scan qz1min = " << qz1low << ", qz1max = " << qz1hig << ", qz2min = " << qz2low << ", qz2max = " << qz2hig
//		  << ", dz1 = " << dz1 << ", dz2 = " << dz2 << endl;
	} 
	if (debug >= 1) {
		cout << " Minimum discr = " << discmin << " for qx1 = " << qx1min << ", qy1 = " << qy1min << endl;
	}

	if (fabs(discmin-discminpre) < discconverge) status = 0;
	if (itry > 100) status = 2;
	Xmin = qx1min;
	Ymin = qy1min;
	return discmin;
	
}

	
double SolveTTbar::BinSrc_qxy (double pb1in[], double pb2in[], double pl1in[], double pl2in[], double pTmin[]) {
// Scans the plane (qx1, qy1) with binary search for a minimum of the discriminant which is returned
//
// input: arrays pi[4] = px, py, pz, E of the visible particles (pbi for b's, pli for leptons)
//		  array  pTm[2] = px, py of the missing transverse momentum

	double qylowmin = -6000.;
	double qyhigmax = 6000.;
	
	double pb1[4], pl1[4], pb2[4], pl2[4], pTm[2];
	for (int i = 0; i < 3; ++i) {
		pb1[i] = pb1in[i];
		pl1[i] = pl1in[i];
		pb2[i] = pb2in[i];
		pl2[i] = pl2in[i];
		if (i < 2) pTm[i] = pTmin[i];
	}
	pb1[3] = sqrt(pb1[0]*pb1[0] + pb1[1]*pb1[1] + pb1[2]*pb1[2]);
	pl1[3] = sqrt(pl1[0]*pl1[0] + pl1[1]*pl1[1] + pl1[2]*pl1[2]);
	pb2[3] = sqrt(pb2[0]*pb2[0] + pb2[1]*pb2[1] + pb2[2]*pb2[2]);
	pl2[3] = sqrt(pl2[0]*pl2[0] + pl2[1]*pl2[1] + pl2[2]*pl2[2]);
	
	double qx1low = xlo;
	double qx1hig = xhi;
	double dx1 = (qx1hig - qx1low) / (float)nscan;
	double dy1 = dx1;
	double discmin = 99999.;
	double discminpre = discmin;
	double qx1min = 0., qy1min = 0.;

	double q1[4] = {0., 0., 0., 0.}, q2[4] = {0., 0., 0., 0.};
	
	int itry = 0;
	while(1) {
		itry++;
		int ifnd = 0;
		for (int i = 0; i < nscan; ++i) {
//			int i = 0;
			double qx1 = qx1low + (float)i * dx1;
//			double qx1 = 63.1534;
			q1[0] = qx1;
			double qy[2];
			int ifl = qyfrqx (pl1, 0., mW, qx1, qy);
			if (!ifl) {
				if (debug >= 2) {
					cout << " ******** ibin = " << i << " qx1 = " << qx1 << " No limits for qy1 " << endl;
				}
				if (ifnd) {
					if (debug >= 2) {
						cout << " stop scan here " << endl;
					}
					break;
				}
				continue;
			}	
			ifnd = 1;
			if (qy[0] < qylowmin) qy[0] = qylowmin;
			if (qy[1] > qyhigmax) qy[1] = qyhigmax;
			double qy1low = qy[0];
			q1[1] = qy1low;
			double disclow = SolveEqn_qxy1 (pb1, pb2, pl1, pl2, pTm, q1, q2);
			double qy1hig = qy[1];
			q1[1] = qy1hig;
			double dischig = SolveEqn_qxy1 (pb1, pb2, pl1, pl2, pTm, q1, q2);
			double qy1mid = 0.5 * (qy1hig + qy1low);
			q1[1] = qy1mid;
			double discmid = SolveEqn_qxy1 (pb1, pb2, pl1, pl2, pTm, q1, q2);
			if (debug >= 3) {
				cout << " qx1 = " << qx1 << endl;
				cout << " qy1low = " << qy1low << " disclow = " << disclow;
				cout << " qy1hig = " << qy1hig << " dischig = " << dischig;
				cout << " qy1mid = " << qy1mid << " discmid = " << discmid << endl;
			}
			double discbst, qy1bst;
			if (disclow < dischig) {
				discbst = disclow;
				qy1bst = qy1low;
			} else {
				discbst = dischig;
				qy1bst = qy1hig;
			}
			double discpre = 99999.;
			if (discmid > disclow && discmid > dischig) {
				if (debug >= 3) {
					cout << " Crazy: disc low = " << disclow << " mid = " << discmid << " hig = " << dischig << endl;
				}
			}
			int nscy = (qy1hig - qy1low) / dy1 + 1.;
			for (int j = 0; j < nscy; ++j) {
				nstps++;
				double qy1up = 0.5 * (qy1mid + qy1hig);
				q1[1] = qy1up;
				double discup = SolveEqn_qxy1 (pb1, pb2, pl1, pl2, pTm, q1, q2);
				double qy1dw = 0.5 * (qy1mid + qy1low);
				q1[1] = qy1dw;
				double discdw = SolveEqn_qxy1 (pb1, pb2, pl1, pl2, pTm, q1, q2);
				if (debug >= 3) {
					cout << " y-step = " << j << "  qx1 = " << qx1 << " qy1mid = " << qy1mid << " discdw = " << discdw << " discup = " << discup
						<< " qy1dw = " << qy1dw << " qy1up = " << qy1up << endl;
				}
				if (discup < discdw) {
					qy1low = qy1mid;
					disclow = discmid;
					qy1mid = qy1up;
					discmid = discup;
				} else {
					qy1hig = qy1mid;
					dischig = discmid;
					qy1mid = qy1dw;
					discmid = discdw;
				}
				discbst = discmid;
				qy1bst = qy1mid;
				if (debug >= 3) {
					cout << " discbst = " << discbst << " discpre = " << discpre << endl;
				}
				if (fabs(discpre-discbst) < 0.005) break;
				discpre = discbst;
			} // end loop over bins in qy1
			if (debug >= 2) {
				cout << " nsteps = " << nstps << endl;
			}
			
			if (discbst < discmin) {
				discmin = discbst;
				qx1min = qx1;
				qy1min = qy1bst;
			}
			if (debug >= 2) {
				cout << " itry = " << itry << " ibin = " << i << " qx1 = " << qx1 << " min discr = " << discmin << " for qy1 = " << qy1min
				<< " discminpre = " << discminpre << endl;
			}
		} // end loop over bins in qx1
		
//		cout << " Minimum discr = " << discmin << " for qx1 = " << qx1min << ", qy1 = " << qy1min << endl;
		nscnstps++;
		if (discmin < discconverge || fabs(discmin-discminpre) < discconverge || itry > 20) {
			if (itry <= 20) status = 0;
			else status = 2;
			break;
		}
//		break;
		if (debug >= 2) {
			cout << " itry = " << itry << " qx1 = " << qx1min << " min discr = " << discmin << " discminpre = " << discminpre << endl;
		}
		discminpre = discmin;
		qx1low = qx1min - 5.*dx1;
		qx1hig = qx1min + 5.*dx1;
		dx1 = (qx1hig - qx1low) / (float)nscan;
		if (debug >= 3) {
			cout << " Reduced window qx1low = " << qx1low << ", qx1hig = " << qx1hig << endl;
		}
	} 
	if (debug > 0) {
		cout << " Minimum discr = " << discmin << " for qx1 = " << qx1min << ", qy1 = " << qy1min << endl;
	}

	if (fabs(discmin-discminpre) < discconverge) status = 0;
	if (itry > 100) status = 2;
	Xmin = qx1min;
	Ymin = qy1min;
	return discmin;
	
}


int SolveTTbar::qyfrqx (double p[], double m, double mWt, double qx, double qy[]) {
// Computes the solutions for qy for a given qx using MT
// returns the 2 solutions in qy[2]	
	
	qy[0] = 0.;
	qy[1] = 0.;
	double ETsq = p[0]*p[0] + p[1]*p[1] + m*m;
	double mdifsq = mWt*mWt - m*m;
	double a = ETsq - p[1]*p[1];
	double b = -(2.*p[0]*qx + mdifsq) * p[1];
	double c = (ETsq - p[0]*p[0])*qx*qx - mdifsq*( p[0]*qx + 0.25*mdifsq );
	double dis = b*b - 4.*a*c;
//	cout << a << " " << b << " " << c << " " << dis << endl;
	if (dis < 0.) return 0;
	double sqrtdis = sqrt(dis);
	qy[0] = (-b - sqrtdis) / (2.*a);
	qy[1] = (-b + sqrtdis) / (2.*a);
//	double test = 2.*qx*sqrt(ETsq) - 2.*p[0]*qx;
//	cout << " test = " << sqrt(test) << endl;

	return 1;
}


double SolveTTbar::SolveEqn_qxy1 (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[], double q1[], double q2[]) {
// Solves the constraint equations and returns a discriminant of type D2
// expects qx1 and qy1 in q1[] to be filled, computes qz1 and qz2
// also completes the arrays q1 and q2
	
	
	double pb1mod = pb1[3];
	double pl1mod = pl1[3];
	double pb2mod = pb2[3];
	double pl2mod = pl2[3];

	// compute the constants
	double A1 = 0.5 * (mt*mt/pb1mod - mW*mW/pl1mod);
	double A2 = 0.5 * (mt*mt/pb2mod - mW*mW/pl2mod);
	double v1[3], v2[3];
	v1[0] = pb1[0] / pb1mod - pl1[0] / pl1mod;
	v1[1] = pb1[1] / pb1mod - pl1[1] / pl1mod;
	v1[2] = pb1[2] / pb1mod - pl1[2] / pl1mod;
	v2[0] = pb2[0] / pb2mod - pl2[0] / pl2mod;
	v2[1] = pb2[1] / pb2mod - pl2[1] / pl2mod;
	v2[2] = pb2[2] / pb2mod - pl2[2] / pl2mod;
	double Ml1b1 = Minv2 (pl1, pb1);
	double Ml2b2 = Minv2 (pl2, pb2);
	double c1 = (Ml1b1*Ml1b1 + mW*mW) / (2.*pb1mod) - A1;
	double c2 = (Ml2b2*Ml2b2 + mW*mW) / (2.*pb2mod) - A2 - v2[0]*pTm[0] - v2[1]*pTm[1];
	
	// Now compute qz1 and qz2
	double qz1 = (-v1[0]*q1[0] - v1[1]*q1[1] + c1) / v1[2];
	double qz2 = ( v2[0]*q1[0] + v2[1]*q1[1] + c2) / v2[2];
	q1[2] = qz1;
	q2[0] = pTm[0] - q1[0];
	q2[1] = pTm[1] - q1[1];
	q2[2] = qz2;
	if (debug >= 4) {
		cout << " q1 px, py, pz = " << q1[0] << ", " << q1[1] << ", " << q1[2] << endl;
		cout << " q2 px, py, pz = " << q2[0] << ", " << q2[1] << ", " << q2[2] << endl;
	}
	
	double q1mod = sqrt(q1[0]*q1[0] + q1[1]*q1[1] + q1[2]*q1[2]);
	double q2mod = sqrt(q2[0]*q2[0] + q2[1]*q2[1] + q2[2]*q2[2]);
	q1[3] = q1mod;
	q2[3] = q2mod;
	
	double mW1 = sqrt(2.*(pl1mod*q1mod - pl1[0]*q1[0] - pl1[1]*q1[1] - pl1[2]*q1[2]));
	double mW2 = sqrt(2.*(pl2mod*q2mod - pl2[0]*q2[0] - pl2[1]*q2[1] - pl2[2]*q2[2]));
//	double dW = ( (mW1-mW)*(mW1-mW) + (mW2-mW)*(mW2-mW) ) / (mW*mW);
	double dW = (mW1-mW)*(mW1-mW) + (mW2-mW)*(mW2-mW);
	double mt1 = sqrt(2.*(pb1mod*q1mod - pb1[0]*q1[0] - pb1[1]*q1[1] - pb1[2]*q1[2]) + Ml1b1*Ml1b1 + mW*mW);
	double mt2 = sqrt(2.*(pb2mod*q2mod - pb2[0]*q2[0] - pb2[1]*q2[1] - pb2[2]*q2[2]) + Ml2b2*Ml2b2 + mW*mW);
//	double dt = ( (mt1-mt)*(mt1-mt) + (mt2-mt)*(mt2-mt) ) / (mt*mt);
	double dt = (mt1-mt)*(mt1-mt) + (mt2-mt)*(mt2-mt);
	double discr = 0.5 * sqrt(dW + dt);
//	cout << " mW1 = " << mW1 << ", mW2 = " << mW2 << ", Disc = " << dW << endl;
//	cout << " mt1 = " << mt1 << ", mt2 = " << mt2 << ", Disc = " << dt << endl;

	return discr;
	
}


int SolveTTbar::SolveEqn4 (double pb1[], double pb2[], double pl1[], double pl2[], double pTm[]) {
// Solves the 4th order constraint equations and returns a discriminant of type D2
// returns the number of solutions
// Implements the Sonnenschein solution
		
	double Eb1 = pb1[3];
	double El1 = pl1[3];
	double Eb2 = pb2[3];
	double El2 = pl2[3];
	
	// use approximation that masses are zero for b, l, nu
	double mb = 0., ml = 0., mn = 0.;
	double mWsq = mW*mW;
	double mtsq = mt*mt;
	double mbsq = mb*mb;
	double mlsq = ml*ml;
	double mnsq = mn*mn;
	double mWln = mWsq-mlsq-mnsq;
	double mtbln = mtsq-mbsq-mlsq-mnsq;
	
	// Compute the coefficients ai
	double a1 = (Eb1+El1)*mWln - El1*mtbln + 2.*Eb1*El1*El1 - 2.*El1*(pb1[0]*pl1[0]+pb1[1]*pl1[1]+pb1[2]*pl1[2]);
	double a2 = 2.*(Eb1*pl1[0]-El1*pb1[0]);
	double a3 = 2.*(Eb1*pl1[1]-El1*pb1[1]);
	double a4 = 2.*(Eb1*pl1[2]-El1*pb1[2]);
	double a14 = a1 / a4;
	double a24 = a2 / a4;
	double a34 = a3 / a4;
	// Compute the coefficients bi
	double b1 = (Eb2+El2)*mWln - El2*mtbln + 2.*Eb2*El2*El2 - 2.*El2*(pb2[0]*pl2[0]+pb2[1]*pl2[1]+pb2[2]*pl2[2]);
	double b2 = 2.*(Eb2*pl2[0]-El2*pb2[0]);
	double b3 = 2.*(Eb2*pl2[1]-El2*pb2[1]);
	double b4 = 2.*(Eb2*pl2[2]-El2*pb2[2]);
	double b14 = b1 / b4;
	double b24 = b2 / b4;
	double b34 = b3 / b4;
	// Compute the coefficients cij
	double c22 = mWln*mWln - 4.*(El1*El1-pl1[2]*pl1[2])*a14*a14 - 4.*mWln*pl1[2]*a14;
	double c21 = 4.*mWln*(pl1[0]-pl1[2]*a24) - 8.*(El1*El1-pl1[2]*pl1[2])*a1*a2/(a4*a4) - 8.*pl1[0]*pl1[2]*a14;
	double c20 = -4.*(El1*El1-pl1[0]*pl1[0]) - 4.*(El1*El1-pl1[2]*pl1[2])*a24*a24 - 8.*pl1[0]*pl1[2]*a24;
	double c11 = 4.*mWln*(pl1[1]-pl1[2]*a34) - 8.*(El1*El1-pl1[2]*pl1[2])*a1*a3/(a4*a4) - 8.*pl1[1]*pl1[2]*a14;
	double c10 = - 8.*(El1*El1-pl1[2]*pl1[2])*a2*a3/(a4*a4) + 8.*pl1[0]*pl1[1] - 8.*pl1[0]*pl1[2]*a34 - 8.*pl1[1]*pl1[2]*a24;
	double c00 = -4.*(El1*El1-pl1[1]*pl1[1]) - 4.*(El1*El1-pl1[2]*pl1[2])*a34*a34 - 8.*pl1[1]*pl1[2]*a34;
	// Compute the coefficients dij prime
	double dp22 = mWln*mWln - 4.*(El2*El2-pl2[2]*pl2[2])*b14*b14 - 4.*mWln*pl2[2]*b14;
	double dp21 = 4.*mWln*(pl2[0]-pl2[2]*b24) - 8.*(El2*El2-pl2[2]*pl2[2])*b1*b2/(b4*b4) - 8.*pl2[0]*pl2[2]*b14;
	double dp20 = -4.*(El2*El2-pl2[0]*pl2[0]) - 4.*(El2*El2-pl2[2]*pl2[2])*b24*b24 - 8.*pl2[0]*pl2[2]*b24;
	double dp11 = 4.*mWln*(pl2[1]-pl2[2]*b34) - 8.*(El2*El2-pl2[2]*pl2[2])*b1*b3/(b4*b4) - 8.*pl2[1]*pl2[2]*b14;
	double dp10 = - 8.*(El2*El2-pl2[2]*pl2[2])*b2*b3/(b4*b4) + 8.*pl2[0]*pl2[1] - 8.*pl2[0]*pl2[2]*b34 - 8.*pl2[1]*pl2[2]*b24;
	double dp00 = -4.*(El2*El2-pl2[1]*pl2[1]) - 4.*(El2*El2-pl2[2]*pl2[2])*b34*b34 - 8.*pl2[1]*pl2[2]*b34;
	// Compute the coefficients dij
	double d22 = dp22 + pTm[0]*pTm[0]*dp20 + pTm[1]*pTm[1]*dp00 + pTm[0]*pTm[1]*dp10 + pTm[0]*dp21 + pTm[1]*dp11;
	double d21 = -dp21 - 2.*pTm[0]*dp20 - pTm[1]*dp10;
	double d20 = dp20;
	double d11 = -dp11 - 2.*pTm[1]*dp00 - pTm[0]*dp10;
	double d10 = dp10;
	double d00 = dp00;
	// Compute the coefficients hi
	double h4 = c00*c00*d22*d22 + c11*d22*(c11*d00-c00*d11) + c00*c22*(d11*d11-2.*d00*d22) + c22*d00*(c22*d00-c11*d11);
	double h3 = c00*d21*(2.*c00*d22-c11*d11) + c00*d11*(2.*c22*d10+c21*d11) + c22*d00*(2.*c21*d00-c11*d10) - c00*d22*(c11*d10+c10*d11)
			- 2.*c00*d00*(c22*d21+c21*d22) - d00*d11*(c11*c21+c10*c22) + c11*d00*(c11*d21+2.*c10*d22);
	double h2 = c00*c00*(2.*d22*d20+d21*d21) - c00*d21*(c11*d10+c10*d11) + c11*d20*(c11*d00-c00*d11) + c00*d10*(c22*d10-c10*d22)
			+ c00*d11*(2.*c21*d10+c20*d11) + d00*d00*(2.*c22*c20+c21*c21) - 2.*c00*d00*(c22*d20+c21*d21+c20*d22)
			+ c10*d00*(2.*c11*d21+c10*d22) - d00*d10*(c11*c21+c10*c22) - d00*d11*(c11*c20+c10*c21);
	double h1 = c00*d21*(2.*c00*d20-c10*d10) - c00*d20*(c11*d10+c10*d11) + c00*d10*(c21*d10+2.*c20*d11) - 2.*c00*d00*(c21*d20+c20*d21) 
			+ c10*d00*(2.*c11*d20+c10*d21) + c20*d00*(2.*c21*d00-c10*d11) - d00*d10*(c11*c20+c10*c21);
	double h0 = c00*c00*d20*d20 + c10*d20*(c10*d00-c00*d10) + c20*d10*(c00*d10-c10*d00) + c20*d00*(c20*d00-2.*c00*d20) ;

	double coef[5];
	coef[0] = h4;
	coef[1] = h3;
	coef[2] = h2;
	coef[3] = h1;
	coef[4] = h0;
	
	// solve the quartic equation
	nsol4 = solveQuartic(coef, pn1x);
	
	for (int i = 0; i < nsol4; ++i) {
		double c0 = c00;
		double c1 = c11 + c10*pn1x[i];
		double c2 = c22 + c21*pn1x[i] + c20*pn1x[i]*pn1x[i];
		double d0 = d00;
		double d1 = d11 + d10*pn1x[i];
		double d2 = d22 + d21*pn1x[i] + d20*pn1x[i]*pn1x[i];
		pn1y[i] = (c0*d2 - c2*d0) / (c1*d0 - c0*d1);
		pn1z[i] = -(a1 + a2*pn1x[i] + a3*pn1y[i]) / a4;
		pn2x[i] = pTm[0] - pn1x[i];
		pn2y[i] = pTm[1] - pn1y[i];
		pn2z[i] = -(b1 + b2*pn2x[i] + b3*pn2y[i]) / b4;
	}
	if (debug >= 2) {
		cout << " Numberr of solutions = " << nsol4 << endl;
		for (int i = 0; i < nsol4; ++i) {
			cout << " " << i << " pn1x =" << pn1x[i] << " pn1y = " << pn1y[i] << " pn1z = " << pn1z[i]
			  << " pn2x =" << pn2x[i] << " pn2y = " << pn2y[i] << " pn2z = " << pn2z[i] << endl;
		}
	}
	
	return nsol4;
	
}


void SolveTTbar::EvtPrint (double pb1[], double pb2[], double pl1[], double pl2[], double pn1[], double pn2[]) {


	double pTm[2];
	pTm[0] = -(pb1[0] + pb2[0] + pl1[0] + pl2[0]);
	pTm[1] = -(pb1[1] + pb2[1] + pl1[1] + pl2[1]);

	cout << " Event print: " << endl;
	cout << " b1 px, py, pz, pmod = " << pb1[0] << ", " << pb1[1] << ", " << pb1[2] << ", " << pb1[3] << endl;
	cout << " l1 px, py, pz, pmod = " << pl1[0] << ", " << pl1[1] << ", " << pl1[2] << ", " << pl1[3] << endl;
	cout << " b2 px, py, pz, pmod = " << pb2[0] << ", " << pb2[1] << ", " << pb2[2] << ", " << pb2[3] << endl;
	cout << " l2 px, py, pz, pmod = " << pl2[0] << ", " << pl2[1] << ", " << pl2[2] << ", " << pl2[3] << endl;
	cout << " pTmiss x, y =   " << pTm[0] << ", " << pTm[1] << endl;
	cout << " Generated neutrinos: " << endl;
	cout << " n1 px, py, pz, pmod = " << pn1[0] << ", " << pn1[1] << ", " << pn1[2] << ", " << pn1[3] << endl;
	cout << " n2 px, py, pz, pmod = " << pn2[0] << ", " << pn2[1] << ", " << pn2[2] << ", " << pn2[3] << endl;
	
	double pTsum[] = {0., 0.};
	for (int i = 0; i < 2; ++i) {
		pTsum[i] = pb1[i] + pb2[i] + pl1[i] + pl2[i] + pn1[i] + pn2[i];
	}
	cout << " pT sum x = " << pTsum[0] << ", y = " << pTsum[1] << endl;

	return;

}


// double Minv2 (vector<double> p1, vector<double> p2) {
double SolveTTbar::Minv2 (double p1[], double p2[]) {
	
	double minvsq = ((p1[3]+p2[3])*(p1[3]+p2[3])  - (p1[0]+p2[0])*(p1[0]+p2[0])
					 - (p1[1]+p2[1])*(p1[1]+p2[1]) - (p1[2]+p2[2])*(p1[2]+p2[2]) );
	if (minvsq < 0.) minvsq = 0.;
	return sqrt(minvsq);
	
}


// double Minv2 (vector<double> p1, vector<double> p2) {
int SolveTTbar::InterParabol (double dz1, double dz2, double D1, double D2, double D3, double dz[]) {
// finds the lowest point of the parabola a(x-x0)^2 + c = D, where x is the length along the joining line
// returns 0 if no solution
	
	if (D2 > D1 && D2 > D3) return 0;
	
	double x3 = sqrt(dz1*dz1 + dz2*dz2);
	double x2 = 0.5 * x3;
	
	double a = (D1 * (x2-x3) + D2*x3 - D3*x2) / (x2*x3*(x2-x3));
	double x0 = ((D1-D2)*x3*x3 - (D1-D3)*x2*x2) / (2.*a*x2*x3*(x3-x2));
	dz[0] = dz1 * (1. - x0 / x3);
	dz[1] = dz2 * x0 / x3;
	
//	double c = D1 - a * x0*x0;
//	cout << " Interp dz1 = " << dz1 << ", dz2 = " << dz2 << ", D1 = " << D1 << ", D2 = " << D2 << ", D3 = " << D3 << endl;
//	cout << "  a = " << a << ", x0 = " << x0 << ", discr min = " << c
//		<< ", dz[0] = " << dz[0] << ", dz[1] = " << dz[1] << endl;

	return 1;

}


int SolveTTbar::solveQuartic(double coef[], double sol[]) {
  // computes the solutions of a quartic equation
  //  uses the Ferrari solution, as presented in wikipedia
  // input: coefficients of
  //   coef[0] + coef[1]*x + coef[2]*x^2 + coef[3]*x^3 + coef[4]* x^4 = 0
  // output: nsol   = number of real solutions (0 to 4)
  //         sol[4] = values of the solutions, in increasing order

  int nsol = 0;
  if (coef[4] == 0.) {
    return nsol;
  }
  sol[0] = 0.;

  // compute coefficients for a depressed quartic
  // x^4 + a x^2 + b x + c = 0
  double a = - 3.*coef[3]*coef[3] / (8.*coef[4]*coef[4]) + coef[2] / coef[4];
  double b = coef[3]*coef[3]*coef[3] / (8.*coef[4]*coef[4]*coef[4])
           - coef[3]*coef[2] / (2.*coef[4]*coef[4]) + coef[1] / coef[4];
  double c = - 3.*coef[3]*coef[3]*coef[3]*coef[3]
               / (256.*coef[4]*coef[4]*coef[4]*coef[4])
           + coef[3]*coef[3]*coef[2] / (16.*coef[4]*coef[4]*coef[4])
           - coef[3]*coef[1] / (4.*coef[4]*coef[4]) + coef[0] / coef[4];
  if (debug > 3){
    cout << " Solve the quartic equation: " << endl;
    cout << " a = " << a << " b = " << b << " c = " << c << endl;
  }
  int n = 0;
  double term1, term2, term3;

  // if the equation is biquartic, solve imediately
  if (b == 0.) {
    double root = a*a - 4.*c;
    if (debug > 3){
      cout << " root = " << root << endl;
    }
    if (root >= 0.) {
      root = sqrt(root);
      term1 = - coef[3] / (4.*coef[4]);
      term2 = 0.5 * (-a + root);
      term3 = 0.5 * (-a - root);
      if (term2 >= 0.) {
        term2 = sqrt(term2);
        sol[n]   = term1 - term2;
        sol[n+1] = term1 + term2;
        n = n + 2;
      }
      if (term3 >= 0.) {
        term3 = sqrt(term3);
        sol[n]   = term1 - term3;
        sol[n+1] = term1 + term3;
        n = n + 2;
      }
    }
    nsol = n;

  // else if not biquartic
  } else {
    // compute the solution of the depressed cubic equation
    double p = - a*a / 12. - c;
    double q = - a*a*a / 108. + a*c / 3. - b*b / 8.;
    double qsq = q*q;
    double p3rd = p*p*p;
    double root = 0.25*qsq + p3rd / 27.;
    double prec = 0.001;
    if (fabs(root) < prec) { root = 0.;}
    if (debug > 3){
      cout << " p = " << p << " q = " << q << " root = " << root << endl;
    }
	  double u = 0.; //, u1, u2;
    if (root > 0.) {
      double rp, rm;
//      if (fabs(p3rd/qsq) > 1.E-7) {
        double sqrtroot = sqrt(root);
        rp = - 0.5*q + sqrtroot;
        rm = - 0.5*q - sqrtroot;
/*      } else if (q > 0.) {
        rp = 1./27. * p3rd / q;
        rm = -q;
      } else {
        rp = -q;
        rm = -1./27. * p3rd / q;
      }
*/
      if (debug > 3){
        cout << " rp,m = " << rp << " " << rm << endl;
      }
      double sgn = rp < 0. ? -1. : 1.;
      u = sgn * exp( log(fabs(rp)) / 3. );
      sgn = rm < 0. ? -1. : 1.;
      u = u + sgn * exp( log(fabs(rm)) / 3. );
    } else if (root <= 0.) {
      double sgn = root < 0. ? -1. : 1.;
      double modul = sqrt(0.25*q*q + sgn*root);
      double argum = 0.;
      if (modul != 0.) {
        argum = acos(-0.5*q/modul);
        u = 2. * exp(log(modul) / 3.) * cos(argum/3.);
      } else {
        u = 0.;
      }
      if (debug > 3){
        cout << " modul = " << modul << " argum = " << argum << endl;
      }
    }
    if (debug > 3){
      cout << " u = " << u << endl;
      double v = u;
      double eqn = v*v*v + p*v + q;
      cout << " check depressed cubic eqn = " << eqn << endl;
    }
    // transform to the solution of the non-depressed cubic equation
    double y;
    if (u != 0.) {
      y = - 5.*a / 6. + u;
//      y = - 5.*a / 6. + u - p / (3.*u);
      if (debug > 3){
        cout << " u = " << u << " y = " << y << endl;
        double eqn = y*y*y + 2.5*a*y*y + (2.*a*a-c)*y
                   + 0.5*a*a*a-0.5*a*c-b*b/8;
        cout << " check non-depressed cubic eqn = " << eqn << endl;
      }
    } else {
      double sgn = q < 0. ? -1. : 1.;
      y = - 5.*a / 6. - sgn * exp( log(fabs(q)) / 3. );
      if (debug > 3){
        cout << " y = " << y << endl;
        double eqn = y*y*y + 2.5*a*y*y + (2.*a*a-c)*y
                   + 0.5*a*a*a-0.5*a*c-b*b/8;
        cout << " check non-depressed cubic eqn = " << eqn << endl;
      }
    }
    // compute the solution of the quartic equation
    double w = a + 2.*y;
    if (debug > 3){
      cout << " wsq = " << w << endl;
    }
    if (w >= 0.) {
      w = sqrt(w);
      term1 = - coef[3] / (4.*coef[4]);
      term2 = - (3.*a + 2.*y + 2.*b/w);
      term3 = - (3.*a + 2.*y - 2.*b/w);
      if (debug > 3){
        cout << " term1 = " << term1 << " term2 = " << term2
             << " term3 = " << term3 << endl;
      }
//      double prec = 0.01;
      if (term2 >= 0.) {
        term2 = sqrt(term2);
        sol[n]   = term1 + 0.5*(w - term2);
        sol[n+1] = term1 + 0.5*(w + term2);
        n = n + 2;
      }
      if (term3 >= 0.) {
        term3 = sqrt(term3);
        sol[n]   = term1 - 0.5*(w - term3);
        sol[n+1] = term1 - 0.5*(w + term3);
        n = n + 2;
      } 
    }
    nsol = n;
  }
  if (debug > 3){
    // check quartic equation
    for (int i = 0; i < nsol; ++i) {
      double b0 = sol[i];
      double eqn = (((coef[4]*b0 + coef[3])*b0
                      + coef[2])*b0 + coef[1])*b0 + coef[0];
      cout << " check quartic: sol " << i << " eqn = " << eqn << endl;
    }
  }

  // now reorder the solutions according to increasing values
  bool changed;
  do {
    changed = false;
    for (int i = 0; i < nsol-1; ++i) {
      if (sol[i] > sol[i+1]) {
        double save = sol[i];
        sol[i] = sol[i+1];
        sol[i+1] = save;
        changed = true;
      }
    }
  } while (changed);

  return nsol;
}

