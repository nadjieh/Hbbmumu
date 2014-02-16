// Search for top in events
#include "helper/TopSearch.h"

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
// #include <iostream.h>

using namespace std;


//****************************************************************************
// Constructor
//***************************************************************************/

TopSearch::TopSearch () : 
debug(1), inAcc(1), withbtag(1), recoW1j(1), recoWlept(1), ckWoverlap(1), rescaleW(1),
chisqWmax(7.), chisqTmax(2.), chisqTWlepmax(0.5), massWjbmin(40.), massWjbmax(130.), chisqW1b(1.), chisqWlept(0.),
dpbypnom(1.), massW(80.4), widthW2j(2.1), widthW1j(10.), massT(172.), widthT(10.)
{
//	int nevt, nlum, nrun, nreadmax = 1;
	
//   if (debug >= 1) {
//     // initialize function for final statistics
//     ttbarStats (0);
//   }
	

 return;
}

//****************************************************************************
// Destructor
//***************************************************************************/

TopSearch::~TopSearch(){
//   //write out final statistics
//     if (debug >= 1) {
//       ttbarStats (2);
//     }
}

//***************************************************************************/
int TopSearch::GoSearch () {

	if (npart <= 0) return 0;
	
	nbJets = GetnbJets ();

	// collect the W candidates
	
	npW = GetWcand ();
	
	// remove overlapping W candidates
	//if ckWoverlap==0, only scale is applied
	
	//	if (ckWoverlap) 
	npW = RemWoverlaps ();
	
	// then look for Top candidates
	
	npT = GetTcand ();
	
	// remove overlapping Top candidates
	
	npT = RemToverlaps ();
	
	// clean up the Ws
	
	CleanWs ();
	
	// give the final interpretation of the event
	
	EvtFinal ();
		
	// collect data for final statistics
// 	if (debug >= 1) {
// 		ttbarStats (1);
// 	}
		

	return 1;
}


//***************************************************************************/
int TopSearch::GetnbJets () {
  int nBJets = 0;
  for (int i = 0; i < npart; ++i) {
    if (inAcc && fabs (evt[i][7]) > 2.4) continue;
    if (evt[i][0] > 1) continue;
    
    if(evt[i][9] > 0){
      if (debug >= 3)
	cout<<" Jet ["<<i<<"] is a bTagged Jet"<<endl;
      nBJets++;
    }
  }
  if (debug >= 3)
    cout<<" This event has "<<nBJets<<" bTagged Jets"<<endl;

  return nBJets;
}

int TopSearch::GetWcand () {

	int npW = 0;

	// loop over the particles and collect the good W candidates
	if (debug >= 3) {
		cout << " Looking for W candidates (with chisq < " << chisqWmax << ")" << endl;
	}
	// first, check if a single jet mass is compatible with the W
	if (recoW1j) {
		int ilablW1j = -999;
		int ilablW1b = -990;
		for (int i = 0; i < npart; ++i) {
			if (inAcc && fabs (evt[i][7]) > 2.4) continue;
			if (evt[i][0] > 1) continue;
			double chi2W = GetChisqW (i, -1);
			int ilabl = 0;
			// either because 2 jets coalesced
			if (chi2W <= chisqWmax) {
				ilabl = ilablW1j;
			// or because the b from top coalesced with a leg from W
			} else if (evt[i][9] > 0 && evt[i][5] > massWjbmin && evt[i][5] < massWjbmax) {
				chi2W = chisqW1b;
				ilabl = ilablW1b;
			} else continue;
			double maW = evt[i][5];
			if (debug >= 3) {
				cout << "  all W ( " << i << " ), M = " << maW << " chisq = " << chi2W;
				if (ilabl == ilablW1b) cout << ", b-tagged";
				cout << endl;
			}
			mW[npW] = maW;
			chisqW[npW] = chi2W;
			ipW[npW][0] = i;
			ipW[npW][1] = ilabl;
			ilabl++;
			npW++;
		}
	}
	// then look for dijet Ws
	for (int i = 0; i < npart-1; ++i) {
		if (inAcc && fabs (evt[i][7]) > 2.4) continue;
		if (evt[i][0] > 1) continue;
		for (int j = i+1; j < npart; ++j) {
			if (inAcc && fabs (evt[j][7]) > 2.4) continue;
			if (evt[j][0] > 1) continue;
			double chi2W = GetChisqW (i, j);
			double maW = InvMass2 (evt[i], evt[j]);
//			cout << " try W " << i << " " << j << ", M = " << maW << " chisq = " << chisqW << endl;
			if (chi2W > chisqWmax) continue;
//			if (mW < mWmin || massW > mWmax) continue;
			if (debug >= 3) {
			  cout << "  all W ( " << i << ", " << j << " ), M = " << maW << " chisq = " << chi2W << endl;
			}
			mW[npW] = maW;
			chisqW[npW] = chi2W;
			ipW[npW][0] = i;
			ipW[npW][1] = j;
			npW++;
		}
	}
	// also look for leptonic Ws
	//compute the MHT vector as a particle (as only scalar MET is given)
	int nlept = 0;
	if (recoWlept) {
		GetMHT(evt[npart]);
//		ilabl = -909;
		for (int i = 0; i < npart; ++i) {
			if (inAcc && fabs (evt[i][7]) > 2.4) continue;
			if (evt[i][0] != 2 && evt[i][0] != 3) continue;
			if (evt[i+1][0] == 2 || evt[i+1][0] == 3) {
				if (debug >= 3) {
					cout << "  Too many leptons, no W leptonic solution" << endl;
				}
				npart -= 3;
				break;
			}
			npart++;
			nlept++;
			int nsol = pzmissfrW (evt[i], evt[npart-1], evt[npart-1], evt[npart]);
			if (nsol == 0) {
			  if (debug >= 3) {
			    cout << "  No solution for leptonic W" << endl;
			  }
			  npart--;
			  continue;
			}
			if (debug >= 3) {
				cout << "  all W ( " << npart-1 << ", " << i << " ), M = " << massW << " leptonic W created " << endl;
				cout << "  all W ( " << npart << ", " << i << " ), M = " << massW << " leptonic W created " << endl;
			}
			mW[npW] = massW;
			chisqW[npW] = chisqWlept;
			ipW[npW][0] = npart-1;
//			ipW[npW][1] = ilabl;
			ipW[npW][1] = i;
			npW++;
			mW[npW] = massW;
			chisqW[npW] = chisqWlept;
			ipW[npW][0] = npart;
//			ipW[npW][1] = ilabl;
			ipW[npW][1] = i;
			npW++;
//			ilabl++;
			npart++;
		}
	}
	
	// order the candidates in increasing chisq
	int changed;
	do {
	  changed = false;
	  for (int i = 0; i < npW-1; ++i) {
	    if (ipW[i][0] >= 1000) continue;
	    for (int j = i+1; j < npW; ++j) {
	      if (chisqW[i] > chisqW[j]) {
		int i1 = ipW[i][0];
		int i2 = ipW[i][1];
		double mass = mW[i];
		double chisq = chisqW[i];
		ipW[i][0] = ipW[j][0];
		ipW[i][1] = ipW[j][1];
		mW[i] = mW[j];
		chisqW[i] = chisqW[j];
		ipW[j][0] = i1;
		ipW[j][1] = i2;
		mW[j] = mass;
		chisqW[j] = chisq;
		changed = true;
	      }
	    }
	  }
	} while (changed);

	return npW;
}
	


//***************************************************************************/
int TopSearch::RemWoverlaps () {
  if(ckWoverlap){
	// tag the overlapping solutions
	for (int i = 0; i < npW-1; ++i) {
		if (ipW[i][0] >= 1000) continue;
//		cout << " " << ipW[i][0] << " " << ipW[i][1] << endl;
		for (int j = i+1; j < npW; ++j) {
//			if (ipW[j][0] < 0) continue;
//			cout << " " << ipW[i][0] << " " << ipW[i][1] << " " << ipW[j][0] << " " << ipW[j][1] << endl;
			if (ipW[i][0] == ipW[j][0] || ipW[i][0] == ipW[j][1] || ipW[i][1] == ipW[j][0]
			    || (ipW[i][1] == ipW[j][1] && ipW[i][1] >= 0 && evt[ipW[i][1]][0] <= 1) ) {
				if (debug >=3) {
					cout << "  overlapping: ( " << ipW[i][0]%1000;
					if (ipW[i][1] >= 0) cout << ", " << ipW[i][1]%1000;
					cout << " ) with ( " << ipW[j][0]%1000;
					if (ipW[j][1] >= 0) cout << ", " << ipW[j][1]%1000;
					cout << " ) -> reject ( " << ipW[j][0]%1000;
					if (ipW[j][1] >= 0) cout << ", " << ipW[j][1]%1000;
					cout << " )" << endl;
				}
				ipW[j][0] += 1000;
				if (ipW[j][1] >= 0) ipW[j][1] += 1000;
			}
		}
	}
}
	// squeeze out the rejected overlaps
	int ilast = 0;
	for (int i = 0; i < npW; ++i) {
		if (ipW[i][0] >= 1000) continue;
		int i0 = ipW[i][0];
		int i1 = ipW[i][1];
		ipW[ilast][0] = i0;
		ipW[ilast][1] = i1;
		mW[ilast] = mW[i];
		double scalfact = 1.;
		if (rescaleW > 0 && (ipW[ilast][1] >= 0 || ipW[ilast][1] == -999) ) scalfact = massW / mW[ilast];

		pW[ilast][1] = scalfact * evt[i0][1];

		if (i1 >= 0) pW[ilast][1] += scalfact * evt[i1][1];


		pW[ilast][2] = scalfact * evt[i0][2];
		if (i1 >= 0) pW[ilast][2] += scalfact * evt[i1][2];
		pW[ilast][3] = scalfact * evt[i0][3];
		if (i1 >= 0) pW[ilast][3] += scalfact * evt[i1][3];
		pW[ilast][4] = scalfact * evt[i0][4];
		if (i1 >= 0) pW[ilast][4] += scalfact * evt[i1][4];
//		ipW[ilast][0] = ipW[i][0];
//		ipW[ilast][1] = ipW[i][1];
		chisqW[ilast] = chisqW[i];
		if (debug >= 3) {
			cout << "  kept W ( " << ipW[ilast][0];
			if (ipW[ilast][1] >= 0) cout << ", " << ipW[ilast][1];
			cout << " ) M " << mW[ilast]
				<< " px " << pW[ilast][1] << " py " << pW[ilast][2] << " pz " << pW[ilast][3] << " E " << pW[ilast][4];
			if (ipW[ilast][1] < -900 && ipW[ilast][1] > -910) cout << " ambig ";
			
			cout<< endl;
		}
		ilast++;
	}
	npW = ilast;

	return npW;
}


//***************************************************************************/
int TopSearch::GetTcand () {

	int npT = 0;

	// then look for top candidates
	if (debug >= 3) {
		cout << " Looking for T candidates (accepted if chisq < " << chisqTmax << ")" << endl;
	}
	npT = 0;
	for (int j = 0; j < npW; ++j) {
	  for (int i = 0; i < npart; ++i) {
	    if (inAcc && fabs (evt[i][7]) > 2.4) continue;
			if (withbtag == 2 && !(evt[i][9] || ipW[j][1] == -990) ) continue;
			if (evt[i][0] > 1) continue;
//			if (!ileft[i]) continue;
			if (i == ipW[j][0] || i == ipW[j][1]) continue;
			double mass = InvMass2 (evt[i], pW[j]);
			// cout<<" E1, x1, y1, z1 "<<evt[i][4]<<", "<<evt[i][1]<<", "<<evt[i][2]<<", "<<evt[i][3]<<endl;
			// cout<<" E2, x2, y2, z2 "<<pW[j][4]<<", "<<pW[j][1]<<", "<<pW[j][2]<<", "<<pW[j][3]<<" mass: "<<mass<<endl;

			double chisq = GetChisqT (pW[j], ipW[j][0], ipW[j][1], i);
			double chisqtot = chisqW[j] + chisq;
			if (debug >= 3) {
				cout << "  all T " << "( " << ipW[j][0];
				if (ipW[j][1] >= 0) cout << ", " << ipW[j][1];
				cout << " ) " << i << ", M " << mass << " chisq = " << chisq << " chisq tot = " << chisqtot;
				if (chisqtot <= chisqTmax) cout << " accepted " << endl;
				else cout << endl;
			}
			if (evt[ipW[j][0]][0] == 1 && chisqtot > chisqTmax) continue;
			else if (ipW[j][1] >= 0 && evt[ipW[j][1]][0] > 1 && chisqtot > chisqTWlepmax) continue;
			ipT[npT][0] = j;
			ipT[npT][1] = i;
			ipT[npT][2] = 0;
			
			if(evt[i][9] > 0)
			  ipT[npT][2]++;
			if(evt[ipW[j][0]][0] == 1 && evt[ipW[j][0]][9] > 0 )
			  ipT[npT][2]++;
			if(evt[ipW[j][1]][0] == 1 && evt[ipW[j][1]][9] > 0 )
			  ipT[npT][2]++;

			mT[npT] = mass;
			chisqT[npT] = chisqtot;
			npT++;
		}
	}
	
	if(withbtag == 3){
	  for (int i = 0; i < npT; ++i){
	    if((nbJets - ipT[i][2]) < 1){
	      ipT[i][0] += 1000;
	      ipT[i][1] += 1000;
	    }
	  }
	}
	
	// order the candidates in increasing chisq
	int changed;
	do {
		changed = false;
		for (int i = 0; i < npT-1; ++i) {
			if (ipT[i][0] >= 1000) continue;
			for (int j = i+1; j < npT; ++j) {
				int isbtagi = 0;
				if (evt[ipT[i][1]][9] > 0) isbtagi = 1;
				else if (ipW[ipT[i][0]][1] == -990) isbtagi = 1;
				int isbtagj = 0;
				if (evt[ipT[j][1]][9] > 0) isbtagj = 1;
				else if (ipW[ipT[j][0]][1] == -990) isbtagj = 1;
				int jbest = 0;
				// choose the  b-tagged (if requested) or best chisq
				if ( (withbtag >= 1 && isbtagj && !isbtagi) ) jbest = 1;
				else if (withbtag == 0 && chisqT[i] > chisqT[j]) jbest = 1;
				else if (chisqT[i] > chisqT[j] && !isbtagj && !isbtagi) jbest = 1;
				else if (chisqT[i] > chisqT[j] && isbtagj && isbtagi) jbest = 1;
				if (jbest) {
					int i1 = ipT[i][0];
					int i2 = ipT[i][1];
					int i3 = ipT[i][2];
					double mass = mT[i];
					double chisq = chisqT[i];
					ipT[i][0] = ipT[j][0];
					ipT[i][1] = ipT[j][1];
					ipT[i][2] = ipT[j][2];
					mT[i] = mT[j];
					chisqT[i] = chisqT[j];
					ipT[j][0] = i1;
					ipT[j][1] = i2;
					ipT[j][2] = i3;
					mT[j] = mass;
					chisqT[j] = chisq;
					changed = true;
				}
			}
		}
	} while (changed);
	
	if (debug >= 3) {
		for (int i = 0; i < npT; ++i) {
			cout << "  ordered T ( " << ipW[ipT[i][0]][0];
			if (ipW[ipT[i][0]][1] >= 0) cout << ", " << ipW[ipT[i][0]][1];
			cout << " ) " << ipT[i][1] << " chisq " << chisqT[i];
			if (evt[ipT[i][1]][9] || ipW[ipT[i][0]][1] == -990) cout << ", with "<<ipT[i][2]<<" b-tagged ";
			cout << endl;
		}
	}
	
	return npT;
}


//***************************************************************************/
int TopSearch::RemToverlaps () {	

	// tag the overlapping solutions
	for (int i = 0; i < npT-1; ++i) {
		if (ipT[i][0] >= 1000) continue;
		for (int j = i+1; j < npT; ++j) {
			if (ipT[j][0] >= 1000) continue;
			int iovl = 0;
//			if (ipT[i][0] == ipT[j][0] || ipT[i][1] == ipT[j][1]
			if (ipT[i][1] == ipT[j][1]) iovl = 1;
			else if (ipT[i][1] == ipW[ipT[j][0]][0] || ipT[i][1] == ipW[ipT[j][0]][1]) iovl = 1;
			else if (ipT[j][1] == ipW[ipT[i][0]][0] || ipT[j][1] == ipW[ipT[i][0]][1]) iovl = 1;
			else if (ipW[ipT[i][0]][0] == ipW[ipT[j][0]][0] || ipW[ipT[i][0]][1] == ipW[ipT[j][0]][1]) iovl = 1;
			else if (ipW[ipT[i][0]][0] == ipW[ipT[j][0]][1] || ipW[ipT[i][0]][1] == ipW[ipT[j][0]][0]) iovl = 1;
			if (iovl) {
				if (debug >= 3) {
//					cout << " " << ipT[i][0] << " " << ipT[j][0] << endl;
					cout << "  overlapping: ( " << ipW[ipT[i][0]][0];
					if (ipW[ipT[i][0]][1] >= 0) cout << ", " << ipW[ipT[i][0]][1];
					cout << " ) " << ipT[i][1] << " with ( " << ipW[ipT[j][0]][0];
					if (ipW[ipT[j][0]][1] >= 0) cout << ", " << ipW[ipT[j][0]][1];
					cout << " ) " << ipT[j][1];
					cout << " -> reject ( " << ipW[ipT[j][0]][0];
					if (ipW[ipT[j][0]][1] >= 0) cout << ", " << ipW[ipT[j][0]][1];
					cout << " ) " << ipT[j][1] << endl;
				}
//				if (evt[ipT[j][1]][9] > 0 && evt[ipT[i][1]][9] <= 0) {
//					cout << " wrong reject" << endl;
//				}
				ipT[j][0] += 1000;
				ipT[j][1] += 1000;
			}
		}
	}
	// squeeze out the rejected overlaps
	int ilast = 0;
	for (int i = 0; i < npT; ++i) {
		if (ipT[i][1] >= 1000) continue;
		ipT[ilast][0] = ipT[i][0];
		ipT[ilast][1] = ipT[i][1];
		ipT[ilast][2] = ipT[i][2];
		mT[ilast] = mT[i];
		chisqT[ilast] = chisqT[i];
		if (debug >= 3) {
			cout << "  kept T ( " << ipW[ipT[ilast][0]][0];
			if (ipW[ipT[ilast][0]][1] >= 0) cout << ", " << ipW[ipT[ilast][0]][1];
			cout << " ) " << ipT[ilast][1] << endl;
		}
		ilast++;
	}
	npT = ilast;	
	
	return npT;
}


//***************************************************************************/
void TopSearch::CleanWs () {
// tags the Ws (+1000) which are in accepted tops or duplicates from leptonic decay
// and tags the leftover jets

	// set all jets as leftover
	for (int i = 0; i < npart; ++i) {
		ileft[i] = 1;
	}
	
	// remove the jets used in top from the list
	for (int i = 0; i < npT; ++i) {
//		if (ipT[i][0] >= 1000) cout << " top 1000" << endl;
//		if (ipW[ipT[i][0]][0] >= 1000) cout << " W 1000" << endl;
		int iW1 = ipW[ipT[i][0]][0];
		int iW2 = ipW[ipT[i][0]][1];
		ileft[iW1] = 0;
		if (iW2 >= 0) ileft[iW2] = 0;
		ileft[ipT[i][1]] = 0;
	}
	
	if (debug >= 3) {
		cout << " Cleaning for lonely Ws:" << endl;
	}
	if ( npW <= 0) return;	
	// reject Ws with particles included in accepted tops
	for (int i = 0; i < npW; ++i) {
		if (ipW[i][0] >= 1000) continue;
		for (int j = 0; j < npT; ++j) {
//			if (ipT[j][0] >= 1000) continue;
			// if the W is from the top, flag it with +1000
			if (ipW[i][0] == ipW[ipT[j][0]][0]%1000 && (ipW[i][1] < 0 || ipW[i][1] == ipW[ipT[j][0]][1]) ) {
				ipW[i][0] += 1000;
				if (debug >= 3) {
					cout << "  reject W ( " << ipW[i][0]%1000;
					if (ipW[i][1] >= 0) cout << ", " << ipW[i][1];
					cout << " ) is from Top" << endl;
				}
				break;
			}
			// if the W is overlapping with a W from top, flag it with +2000
			else if(ipW[i][0] == ipT[j][1] || ipW[i][1] == ipT[j][1]
				|| ipW[i][0] == ipW[ipT[j][0]][0]%1000 || ipW[i][1] == ipW[ipT[j][0]][1]
				|| ipW[i][0] == ipW[ipT[j][0]][1] || ipW[i][1] == ipW[ipT[j][0]][0]%1000 ) {
//				cout << " " << ipW[i][0] << " " << ipW[i][1] << " " << ipW[ipT[j][0]][0] << " " << ipW[ipT[j][0]][1] << endl;
				ipW[i][0] += 2000;
//				cout << " " << ipW[i][0] << " " << ipW[i][1] << " " << ipW[ipT[j][0]][0] << " " << ipW[ipT[j][0]][1] << endl;
				if (debug >= 3) {
					cout << "  remove W ( " << ipW[i][0]%1000;
					if (ipW[i][1] >= 0) cout  << ", " << ipW[i][1];
					cout << " ) " << endl;
				}
				break;
			}
		}
	}
	
	// also reject lonely Ws with b-tag, flag them with +2000
	for (int i = 0; i < npW; ++i) {
		if (ipW[i][0] >= 1000) continue;
		if (ipW[i][1] == -990) {
			ipW[i][0] += 2000;
			if (debug >= 3) {
					cout << "  remove W ( " << ipW[i][0]%1000;
					if (ipW[i][1] >= 0) cout  << ", " << ipW[i][1];
					cout << " ) b-tag" << endl;
			}
		}
	}
	
	// also reject ambiguous W from leptonic decay, flag them with +2000
	int remove = 0;
//	for (int i = npW-2; i < npW; ++i) {
	if (npW > 1) {
		for (int i = 0; i < 2; ++i) {
			if (ipW[i][0] >= 1000) continue;
			if (i == 0 && ipW[i][1] == ipW[i+1][1] && ipW[i+1][0] >= 1000) remove = 1;
			if (i == 1 && ipW[i][1] == ipW[i-1][1]) remove = 1;
			if (remove) {
				ipW[i][0] += 2000;
				if (debug >= 3) {
					cout << "  remove W ( " << ipW[i][0]%1000;
					if (ipW[i][1] >= 0) cout  << ", " << ipW[i][1];
					cout << " ) amb sol" << endl;
				}
			}
		}
	}
	// squeeze out the rejected Ws with flag = +2000
	int ilast = 0;
	for (int i = 0; i < npW; ++i) {
		if (ipW[i][0] >= 2000) continue;
		for (int j = 0; j < npT; ++j) {
			if (ipT[j][0] == i) ipT[j][0] = ilast;
		}
		int i0 = ipW[i][0];
		int i1 = ipW[i][1];
		ipW[ilast][0] = i0;
		ipW[ilast][1] = i1;
		mW[ilast] = mW[i];
		pW[ilast][1] = pW[i][1];
		pW[ilast][2] = pW[i][2];
		pW[ilast][3] = pW[i][3];
		pW[ilast][4] = pW[i][4];
		chisqW[ilast] = chisqW[i];
		ilast++;
	}
	npW = ilast;
	

	//Saeid:: This part flags the jets in the lonely W's committed for the time being.
// 	// remove the jets used in top from the list
// 	for (int i = 0; i < npW; ++i) {
// 		if (ipW[i][0] >= 1000 && ipW[i][0] < 2000) continue;
// 		ileft[ipW[i][0]] = 0;
// 		if (ipW[i][1] >= 0) ileft[ipW[i][1]] = 0;
// 	}
	
	// consistency check of W tagging (only to make sure things are right)
	int nWfrT = 0, nWflg = 0;
	for (int i = 0; i < npT; ++i) {
		nWfrT++;
//		cout << "  final T ( " << ipW[ipT[i][0]][0];
//		if (ipW[ipT[i][0]][1] >= 0) cout << ", " << ipW[ipT[i][0]][1];
//		cout << " ) " << ipT[i][1] << endl;
	}
	for (int i = 0; i < npW; ++i) {
		if (ipW[i][0] >= 1000) nWflg++;
	}
	if (nWfrT != nWflg) {
		cout << " Inconsistent W tagging: from T = " << nWfrT << " tagged +1000 = " << nWflg << endl;
		for (int i = 0; i < npW; ++i) {
			if (ipW[i][0] < 1000) continue;
			int fnd = -1;
			for (int j = 0; j < npT; ++j) {
				if (ipT[j][0] == i) break;
				fnd = i;
			}
			if (fnd >= 0) cout << "  W found " << fnd << " " << ipW[i][0] << " " << ipW[i][1] << endl;
		}
	}

	
	return;
}

//***************************************************************************/
void TopSearch::EvtFinal () {

	// give the final interpretation of the event
	// list the tops
	if (debug >= 2) {
		cout << endl;
		cout << " Final event interpretation: " << endl;
		cout << " Number of Top candidates = " << npT << endl;
	}
	int ipb = -1, ip1 = -1, ip2 = -1;
	int ihemib, ihemi1, ihemi2;
	int nT = npT;
	if (npT > 2) nT = 2;
	for (int i = 0; i < npT; ++i) {
		if (debug >= 2) {
			cout << "  top: ( " << ipW[ipT[i][0]][0]%1000;
			if (ipW[ipT[i][0]][1] >= 0) cout << ", " << ipW[ipT[i][0]][1];
			cout << " ) " << ipT[i][1]
				<< ", M_t = " << mT[i] << " M_W = " << mW[ipT[i][0]] << " chisq = " << chisqT[i] << endl;
		}
		
		// handle the b-tags
		// and flag the particles according to their origin: b or W
		ipb = ipT[i][1];
		ip1 = ipW[ipT[i][0]][0]%1000;
		if (evt[ipb][9]) {
			evt[ipb][9] += 10; // correct b from top
			if (debug >= 2) {
				cout << "   with correct b-tag" << endl;
			}
		}
		if (ipW[ipT[i][0]][1] == -990 && !evt[ipb][9]) {
			evt[ip1][9] += 10; // correct b from top
			if (debug >= 2) {
				cout << "   with correct b-tag" << endl;
			}
		}
		if (ipW[ipT[i][0]][1] == -990 && evt[ipb][9]) {
			evt[ip1][9] += 20; // correct b from top
			if (debug >= 2) {
				cout << "   with 2 correct b-tags" << endl;
			}
		}
		if (ipW[ipT[i][0]][1] >= 0) {
			if (evt[ip1][9]) evt[ip1][9] += 30; // b from W2j from top
			ip2 = ipW[ipT[i][0]][1];
			if (evt[ip2][9]) evt[ip2][9] += 30;
		} else if (ipW[ipT[i][0]][1] == -999) {
			if (evt[ip1][9]) evt[ip1][9] += 40; // b from W1j from top
		}

		// handle the hemi association
		/*ihemib = IsInHemi (ipb);
		ihemi1 = IsInHemi (ip1);
		if (ipW[ipT[i][0]][1] >= 0) ihemi2 = IsInHemi (ip2);
		else if (ipW[ipT[i][0]][1] == -990) ihemi2 = ihemib;
		else ihemi2 = ihemi1;
		if(ihemib == ihemi1 && ihemi1 == ihemi2) {
			if (debug >= 2) {
				cout << "   all in same hemi" << endl;
			}
			nhTOK[nT]++;
		} else if(ihemib != ihemi1 && ihemi1 == ihemi2) {
			if (debug >= 2) {
				cout << "   b hemi split from W" << endl;
			}
			nhbinT[nT]++;
		} else if(ihemi1 != ihemi2) {
			if (debug >= 2) {
				cout << "   W split in 2 hemis" << endl;
			}
			nhWinT[nT]++;
			}*/
	}
	if (debug >= 2) {
		cout << " Lonely Ws : " << endl;
	}
	// use only lonely Ws, i.e. not flagged by +1000 as from top
	for (int i = 0; i < npW; ++i) {
		if (ipW[i][0] >= 1000) continue;
		ip1 = ipW[i][0];
		ip2 = ipW[i][1];
//		if (i < npW-1 && ip2 == ipW[i+1][1] && ipW[i+1][0] >= 1000) continue;
//		if (i == npW-1 && ip2 == ipW[i-1][1]) continue;
		if (debug >= 2) {
			cout << "  W: ( " << ipW[i][0];
			if (ipW[i][1] >= 0) cout << ", " << ipW[i][1];
			cout << " ), M = " << mW[i] << " chisq = " << chisqW[i] << endl;
		}
		// handle the b-tags
		if (ip2 >= 0) {
			if (evt[ip1][9]) evt[ip1][9] += 50;		// b from lonely W2j
			if (evt[ip2][9]) evt[ip2][9] += 50;
		} else {
			if (evt[ip1][9]) evt[ip1][9] += 60;		// b from lonely W1j
		}
		// handle the hemi association
		/*	ihemi1 = IsInHemi (ip1);
		if (ipW[i][1] >= 0) {ihemi2 = IsInHemi (ip2);}
		else {ihemi2 = ihemi1;}
		if(ihemi1 == ihemi2) {
			if (debug >= 2) {
				cout << "   all in same hemi" << endl;
			}
			nhWOK[nT]++;
		} else if (ipW[i][1] >= 0) {
			if (debug >= 2) {
				cout << "   W split in 2 hemis" << endl;
			}
			nhW[nT]++;
		}*/
	}
	if (debug >= 2) {
		cout << " Lonely jets/leptons : ";
	}
	for (int i = 0; i < npart; ++i) {
		if (!ileft[i]) continue;
//		if (evt[i][5] < 50.) continue;
		if (debug >= 2) {
			cout << " " << i;
		}
		if (evt[i][9]) evt[i][9] += 70;  // b from lonely jet
	}
	if (debug >= 2) {
		cout << endl;
		cout << " B-tagged Jets : " << endl;
	}
	int nb = 0;
	for (int i = 0; i < npart; ++i) {
		if (evt[i][9] == 0) continue;
		if (debug >= 2) {
			cout << "  jet " << i << " is b-tagged";
			if (evt[i][9] == 11) cout << ", is b from top";
			if (evt[i][9] == 21) cout << ", is 2 b from top";
			if (evt[i][9] == 31) cout << ", is b in W2j from top";
			if (evt[i][9] == 41) cout << ", is b in W1j from top";
			if (evt[i][9] == 51) cout << ", is b in lonely W2j";
			if (evt[i][9] == 61) cout << ", is b in lonely W1j";
			if (evt[i][9] == 71) cout << ", is b in lonely jet";
			cout << endl;
		}
		if (evt[i][9] > 0 && evt[i][9] < 11) {
			cout << "  lonely b-tag, jet " << i << ", tag = " << evt[i][9];
		}
		nb++;
	}
	if (nb > 4) nb = 4;
	nbevt[nb]++;
	if (debug >= 2) {
		cout << endl;
	}
	

	return;
}


//***************************************************************************/
int TopSearch::ttbarStats (int ifl) {

//	int nbT[] = {0, 0, 0};
//	int nbW[] = {0, 0, 0};
//	int nbjet[] = {0, 0, 0};
	
	char buff[100];
	
//	cout << " " << ifl << " ntop = " << ntop[0] << " " << ntop[1] << " " << ntop[2] << endl;
	if (ifl == 0) {
		// initialize counters for the final statistics
		nread = 0;
		for (int i = 0; i < 3; ++i) {
			ntop[i] = 0;
			naddW[i] = 0;
			nWinT[i] = 0;
			nWlone[i] = 0;
			nbtot[i] = 0;
			nbinT[i] = 0;
			nbinTbW1b[i] = 0;
			nbinW1jT[i] = 0;
			nbinW2jT[i] = 0;
			nbinW1j[i] = 0;
			nbinW2j[i] = 0;
			nbinJ[i] = 0;
			nhTOK[i] = 0;
			nhbinT[i] = 0;
			nhWinT[i] = 0;
			nhWOK[i] = 0;
			nhW[i] = 0;
			njlone[i] = 0;
			ntop1[i] = 0;
			ntop2[i] = 0;
			nbevt[i] = 0;
		}
		ntop1[3] = 0;
		ntop2[3] = 0;
		ntop2[4] = 0;
		ntop2[5] = 0;
		ntop2[6] = 0;
		ntop2[7] = 0;
		ntop2[8] = 0;
		ntop2[9] = 0;
		nbevt[3] = 0;
		nbevt[4] = 0;
		
		cout << endl;
		if(debug > 1){
		cout << "*************************************************" << endl;
		cout << "**   Looking for top quarks   *******************" << endl;
		cout << "*************************************************" << endl;
		cout << endl;
		cout << " Default values: " << endl;
		cout << " Debug level = " << debug << endl;
		if (inAcc) cout << " Particles only if within eta acceptance" << endl;
		if (withbtag == 0) cout << " Tops ordered by chisquared only" << endl;
		else if (withbtag == 1) cout << " Correct b-tag for the top preferred over best chisq" << endl;
		else if (withbtag == 2) cout << " Correct b-tag required for the top" << endl;
		else if (withbtag == 3) cout << " At least 1 bJet out of the top" << endl;
		if (rescaleW == 0) cout << " No rescaling of W 4-vector done" << endl;
		else if (rescaleW == 1) cout << " W 4-vector rescaled to the correct W mass" << endl;
		if (recoW1j) {
			cout << " Ws reconstructed from single heavy jets" << endl;
			cout << "   mass range for W1b = " << massWjbmin << " to " << massWjbmax << " GeV" << endl;
			cout << "   chisquared penalty for W1b = " << chisqW1b << endl;
		}
		else cout << " No Ws reconstructed from single heavy jets" << endl;
		if (recoWlept) {
			cout << " Ws reconstructed from single leptons" << endl;
			cout << "   chisquared penalty for Wlept = " << chisqWlept << endl;
		}
		else cout << " No Ws reconstructed from single leptons" << endl;
		if (ckWoverlap) cout << " W overlaps are cleaned" << endl;
		else cout << " W overlaps are not cleaned" << endl;
		
		cout << " Chisq cuts: W = " << chisqWmax << ", top = " << chisqTmax << endl;
		cout << " Momentum resolution constant = " << dpbypnom << endl;
		cout << " W mass = " << massW << ", width 2j = " << widthW2j << ", width 1j = " << widthW1j << " GeV" << endl;
		cout << " T mass = " << massT << ", width = " << widthT << " GeV" << endl;
		cout << endl;}
	
	} else if (ifl == 1) {
		// count the top events and classify the tops
		nread++;
		int nT = npT;
		if (npT > 2) nT = 2;
		ntop[nT]++;
		int iTclass[] = {0, 0, 0};
		for (int i = 0; i < nT; ++i) {
			if (ipW[ipT[i][0]][1] == -999) iTclass[i] = 1;				// 1 jet W
			else if (ipW[ipT[i][0]][1] == -990) iTclass[i] = 2;			// 1 b-jet W
			else if (evt[ipW[ipT[i][0]][1]][0] == 1) iTclass[i] = 3;	// 2 jets W
			else if (evt[ipW[ipT[i][0]][1]][0] > 1) iTclass[i] = 4;		// lept W
		}
		if (npT == 1) {
			ntop1[iTclass[0]-1]++;
		} else if (npT >= 2) {
			if (iTclass[0] == 1 && iTclass[1] == 1) ntop2[0]++;
			else if (iTclass[0] == 1 && iTclass[1] == 2) ntop2[1]++;
			else if (iTclass[0] == 2 && iTclass[1] == 1) ntop2[1]++;
			else if (iTclass[0] == 1 && iTclass[1] == 3) ntop2[2]++;
			else if (iTclass[0] == 3 && iTclass[1] == 1) ntop2[2]++;
			else if (iTclass[0] == 1 && iTclass[1] == 4) ntop2[3]++;
			else if (iTclass[0] == 4 && iTclass[1] == 1) ntop2[3]++;
			else if (iTclass[0] == 2 && iTclass[1] == 2) ntop2[4]++;
			else if (iTclass[0] == 2 && iTclass[1] == 3) ntop2[5]++;
			else if (iTclass[0] == 3 && iTclass[1] == 2) ntop2[5]++;
			else if (iTclass[0] == 2 && iTclass[1] == 4) ntop2[6]++;
			else if (iTclass[0] == 4 && iTclass[1] == 2) ntop2[6]++;
			else if (iTclass[0] == 3 && iTclass[1] == 3) ntop2[7]++;
			else if (iTclass[0] == 3 && iTclass[1] == 4) ntop2[8]++;
			else if (iTclass[0] == 4 && iTclass[1] == 3) ntop2[8]++;
			else if (iTclass[0] == 4 && iTclass[1] == 4) ntop2[9]++;
		}
		// count the events with additional W
		for (int i = 0; i < npW; ++i) {
			if (ipW[i][0] >= 1000) continue;
			naddW[nT]++;
			break;
		}
		// count W jet multiplicities if from Top and lonely
		for (int i = 0; i < npW; ++i) {
			if (ipW[i][0] >= 1000) {
				if (ipW[i][1] == -999) {nWinT[0]++;}
				else if (ipW[i][1] == -990) {nWinT[1]++;}
				else {nWinT[2]++;}
			} else {
				if (ipW[i][1] == -999) {nWlone[0]++;}
				else if (ipW[i][1] == -990) {nWlone[1]++;}
				else {nWlone[2]++;}
			}
		}
		// count the lonely jets
		for (int i = 0; i < npart; ++i) {
//			if (evt[i][5] < 50.) continue;
			if (ileft[i]) njlone[nT]++;
		}
		// count the b-tagged jets
		for (int i = 0; i < npart; ++i) {
			if (evt[i][9] == 0) continue;
			nbtot[nT]++;
			if (evt[i][9] == 11) nbinT[nT]++;
			if (evt[i][9] == 21) nbinTbW1b[nT]++;
			if (evt[i][9] == 31) nbinW2jT[nT]++;
			if (evt[i][9] == 41) nbinW1jT[nT]++;
			if (evt[i][9] == 51) nbinW2j[nT]++;
			if (evt[i][9] == 61) nbinW1j[nT]++;
			if (evt[i][9] == 71) nbinJ[nT]++;
		}
	
	} else if (ifl == 2) {
		cout << " Final statistics" << endl;
		cout << " ________________" << endl;
		cout << endl;
		cout << " Total number of events processed = " << nread << endl;
		cout << endl;
		cout << " Note: events with >2 tops are included in 2t" << endl;
		cout << endl;

		cout << " Distribution of top and W events: " << endl;
		cout <<       "               0t     1t     2t" << endl;
//		cout << " ntop = " << ntop[0] << " " << ntop[1] << " " << ntop[2] << endl;
		sprintf(buff, " events      %.4d   %.4d   %.4d", ntop[0], ntop[1], ntop[2]);
		cout << buff << endl;
		sprintf(buff, " with add Ws %.4d   %.4d   %.4d", naddW[0], naddW[1], naddW[2]);
		cout << buff << endl;
		cout << endl;

		cout << " Categories of Ws for reconstructed tops: " << endl;
		cout <<       "       1j     1b     2j   lept" << endl;
		sprintf(buff, " 1t  %.4d   %.4d   %.4d   %.4d", ntop1[0], ntop1[1], ntop1[2], ntop1[3]);
		cout << buff << endl;
		cout << endl;
		cout <<       "    1j-1j  1j-1b  1j-2j 1j-lep  1b-1b  1b-2j 1b-lep  2j-2j 2j-lep lep-lep" << endl;
		sprintf(buff, " 2t  %.4d   %.4d   %.4d   %.4d   %.4d   %.4d   %.4d   %.4d   %.4d   %.4d",
			    ntop2[0], ntop2[1], ntop2[2], ntop2[3], ntop2[4], ntop2[5], ntop2[6], ntop2[7], ntop2[8], ntop2[9]);
		cout << buff << endl;
		cout << endl;

		cout << " Top reconstruction efficiencies (if 2 had tops produced): " << endl;
		cout << "  from N_2t events p = " << sqrt((float)ntop[2] / (float)nread) << endl;
		double rt = sqrt(1.-2.*(float)ntop[1] / (float)nread);
		double pup = 0.5 * (1. + rt);
		double pdw = 0.5 * (1. - rt);
		cout << "  from N_1t events(p or 1-p) = " << pup << " " << pdw << endl;
		cout << "  from N_0t events p = " << sqrt(1. - (float)ntop[0] / (float)nread) << endl;
		cout << "  from N_1t + 2 N_2t events p = " << (float)(ntop[1]+2*ntop[2]) / (float)(2*nread) << endl;
		cout << " Top reconstruction efficiencies (if 1 had top produced): " << endl;
		cout << "  from N_1t + N_2t events p = " << (float)(ntop[1]+ntop[2]) / (float)(nread) << endl;
		cout << "  from N_2t events fake = " << (float)ntop[2] / (float)(nread) << endl;
		cout << " Top reconstruction efficiencies (if 1 had + 1 lept top produced): " << endl;
		int N10 = ntop1[0] + ntop1[1] + ntop1[2] + ntop2[0] + ntop2[1] + ntop2[2] + ntop2[4] + ntop2[5] + ntop2[7];
		int N01 = ntop1[3] + ntop2[9];
		int N11 = ntop2[3] + ntop2[6] + ntop2[8];
		cout << "  from N_10 + N_11 events p_had = " << (float)(N10 + N11) / (float)(nread) << endl;
		cout << "  from N_01 + N_11 events q_lep = " << (float)(N01 + N11) / (float)(nread) << endl;
		cout << endl;

		cout << " Distribution of W decay jet multiplicities: " << endl;
		cout <<       "               1j     1b     2j" << endl;
		sprintf(buff, " W from Top   %.3d    %.3d    %.3d", nWinT[0], nWinT[1], nWinT[2]);
		cout << buff << endl;
		sprintf(buff, " W lonely     %.3d    %.3d    %.3d", nWlone[0], nWlone[1], nWlone[2]);
		cout << buff << endl;
		cout << endl;
		
		cout << " Distribution of lonely jets events: " << endl;
		cout <<       "               0t     1t     2t" << endl;
		sprintf(buff, " #jets        %.3d    %.3d    %.3d", njlone[0], njlone[1], njlone[2]);
		cout << buff << endl;
		double fjevt[5];
		for (int i = 0; i < 3; ++i) {
			fjevt[i] = (float)njlone[i] / ntop[i];
		}
		sprintf(buff, " jets/evt    %0.3f  %0.3f  %0.3f", fjevt[0], fjevt[1], fjevt[2]);
		cout << buff << endl;
		cout << endl;
		
		cout << " Distribution of b-tagged jets per event: " << endl;
		cout <<       "              0b     1b     2b     3b   >=4b" << endl;
		sprintf(buff, " b tags/evt  %.3d    %.3d    %.3d    %.3d    %.3d", nbevt[0], nbevt[1], nbevt[2], nbevt[3], nbevt[4]);
		cout << buff << endl;
		double ntotb = nbevt[0] + nbevt[1] + nbevt[2] + nbevt[3] + nbevt[4];
		double fbevt[5];
		for (int i = 0; i < 5; ++i) {
			fbevt[i] = (float)nbevt[i] / ntotb;
		}
		sprintf(buff, " fraction  %0.3f  %0.3f  %0.3f  %0.3f  %0.3f", fbevt[0], fbevt[1], fbevt[2], fbevt[3], fbevt[4]);
		cout << buff << endl;
		cout << endl;
		
		cout << " Association of b-tagged jets: " << endl;
		cout <<       "                         0t     1t     2t" << endl;
		sprintf(buff, " b total                %.3d    %.3d    %.3d", nbtot[0], nbtot[1], nbtot[2]);
		cout << buff << endl;
		sprintf(buff, " b from t (b or W1b)    %.3d    %.3d    %.3d", nbinT[0], nbinT[1], nbinT[2]);
		cout << buff << endl;
		sprintf(buff, " b from t (b and W1b)   %.3d    %.3d    %.3d", nbinTbW1b[0], nbinTbW1b[1], nbinTbW1b[2]);
		cout << buff << endl;
		sprintf(buff, " b from W2j in t        %.3d    %.3d    %.3d", nbinW2jT[0], nbinW2jT[1], nbinW2jT[2]);
		cout << buff << endl;
		sprintf(buff, " b from W1j in t        %.3d    %.3d    %.3d", nbinW1jT[0], nbinW1jT[1], nbinW1jT[2]);
		cout << buff << endl;
		sprintf(buff, " b from lonely W2j      %.3d    %.3d    %.3d", nbinW2j[0], nbinW2j[1], nbinW2j[2]);
		cout << buff << endl;
		sprintf(buff, " b from lonely W1j      %.3d    %.3d    %.3d", nbinW1j[0], nbinW1j[1], nbinW1j[2]);
		cout << buff << endl;
		sprintf(buff, " b from lonely jet      %.3d    %.3d    %.3d", nbinJ[0], nbinJ[1], nbinJ[2]);
		cout << buff << endl;
		cout << endl;
		
		cout << " Hemisphere association: " << endl;
		cout <<       "                      0t     1t     2t" << endl;
		cout << " For top: " << endl;
		sprintf(buff, "  b and W same hemi  %.3d    %.3d    %.3d", nhTOK[0], nhTOK[1], nhTOK[2]);
		cout << buff << endl;
		sprintf(buff, "  b split from W     %.3d    %.3d    %.3d", nhbinT[0], nhbinT[1], nhbinT[2]);
		cout << buff << endl;
		sprintf(buff, "  W split in 2       %.3d    %.3d    %.3d", nhWinT[0], nhWinT[1], nhWinT[2]);
		cout << buff << endl;
		cout << " For lonely Ws: (counts #Ws)" << endl;
		sprintf(buff, "  both in same hemi  %.3d    %.3d    %.3d", nhWOK[0], nhWOK[1], nhWOK[2]);
		cout << buff << endl;
		sprintf(buff, "  W split in 2       %.3d    %.3d    %.3d", nhW[0], nhW[1], nhW[2]);
		cout << buff << endl;
		cout << endl;
	}

	return 0;
}


//***************************************************************************/
int TopSearch::IsInHemi (int ind) {

	int ihemi = -1;

	int isrc = (int)ind;

	for (int ih = 0; ih < 2; ++ih) {

		for (int i = 0; i < npsjet[ih]; ++i) {

			if (isrc != ipsjet[ih][i]) continue;
			ihemi = ih;
			break;
		}
	}
	
	return ihemi;
}


//***************************************************************************/
double TopSearch::InvMass2 (double pvect1[], double pvect2[])	 {

	 double msq =
    (pvect1[4]+pvect2[4])*(pvect1[4]+pvect2[4])
  - (pvect1[1]+pvect2[1])*(pvect1[1]+pvect2[1])
  - (pvect1[2]+pvect2[2])*(pvect1[2]+pvect2[2])
  - (pvect1[3]+pvect2[3])*(pvect1[3]+pvect2[3]);
 double mass = sqrt(fabs(msq));
 if (msq < 0.) mass = -mass;

 return mass;
}


//***************************************************************************/
double TopSearch::InvMass3 (double pvect1[], double pvect2[], double pvect3[]) {

 double msq =
    (pvect1[4]+pvect2[4]+pvect3[4])*(pvect1[4]+pvect2[4]+pvect3[4])
  - (pvect1[1]+pvect2[1]+pvect3[1])*(pvect1[1]+pvect2[1]+pvect3[1])
  - (pvect1[2]+pvect2[2]+pvect3[2])*(pvect1[2]+pvect2[2]+pvect3[2])
  - (pvect1[3]+pvect2[3]+pvect3[3])*(pvect1[3]+pvect2[3]+pvect3[3]);
 double mass = sqrt(fabs(msq));
 if (msq < 0.) mass = -mass;

 return mass;
}

//***************************************************************************/
double TopSearch::GetChisqW (int i, int j) {
// convention: j = -1 if single particle for W

	double chisqW = -1.;
	double mass, dmass;
	
	if (j < 0) {
		double pmom = sqrt(evt[i][1]*evt[i][1] + evt[i][2]*evt[i][2] + evt[i][3]*evt[i][3]);
		double dpbyp1 = dpbypnom / sqrt(pmom);
		mass = evt[i][5];
		double fact = (evt[i][4] - pmom) * pmom / mass;
		dmass = sqrt(fact*fact*dpbyp1*dpbyp1 + widthW1j*widthW1j);
	} else {
		mass = InvMass2 (evt[i], evt[j]);
		double dpbyp1 = dpbypnom / sqrt(evt[i][4]);
		double dpbyp2 = dpbypnom / sqrt(evt[j][4]);
		dmass = sqrt(0.25*(dpbyp1*dpbyp1 + dpbyp2*dpbyp2) * mass*mass + widthW2j*widthW2j);
	}
	chisqW = (mass-massW)*(mass-massW) / (dmass*dmass);
	
	return chisqW;
}

//***************************************************************************/
double TopSearch::GetChisqT (double pWcor[], int i, int j, int k) {
// pWcor[5] W 4-vector, possibly rescaled to the correct W mass
// i = particle 0 from W (quark or MET)
// j = particle 1 from W (quark or nothing or lepton)
// k = additional jet
// convention: j = -999 to -910 if single particle for W

	double chisqT = -1.;
	double mass, dmass;
	
	// for W from 1j
	if (j < -910) {
		double pmom1 = sqrt(evt[i][1]*evt[i][1] + evt[i][2]*evt[i][2] + evt[i][3]*evt[i][3]);
		mass = InvMass2 (pWcor, evt[k]);
		double fact = 2. * (evt[i][4] - pmom1) * pmom1;
		double mdiff = mass*mass - evt[i][5]*evt[i][5];
		double term1 = (fact + mdiff) * (fact + mdiff);
		double term2 = mdiff * mdiff;
		double dpbyp1 = dpbypnom / sqrt(pmom1);
		double dpbyp3 = dpbypnom / sqrt(evt[k][4]);
		dmass = sqrt(0.25*(term1*dpbyp1*dpbyp1 + term2*dpbyp3*dpbyp3)/(mass*mass) + widthT*widthT);
//		cout << " " << i << " " << j << " " << k << " M = " << mass << " dM = " << dmass << endl;
	} 
	
	// for W leptonic
	else if (evt[j][0] > 1) {
//		cout << " " << i << " " << j << " " << k << endl;
		mass = InvMass3 (evt[i], evt[j], evt[k]);
		double dpbyp1 = dpbypnom * sqrt(HT) / sqrt(evt[i][4]);
		double dpbyp3 = dpbypnom / sqrt(evt[k][4]);
//		cout << " E " << evt[i][4] << endl;
//		cout << " dpbyp " << dpbyp1 << " " << dpbyp3 << endl;
		dmass = sqrt(0.25*(dpbyp1*dpbyp1 + dpbyp3*dpbyp3) * mass*mass + widthT*widthT);
//		cout << " mass " << mass << " " << dmass << endl;
	}
	
	// for W from 2j
	 else {
		mass = InvMass2 (pWcor, evt[k]);
		double mass12 = InvMass2 (evt[i], evt[j]);
		double mass23 = InvMass2 (evt[j], evt[k]);
		double mass31 = InvMass2 (evt[i], evt[k]);
		double dpbyp1 = dpbypnom / sqrt(evt[i][4]);
		double dpbyp2 = dpbypnom / sqrt(evt[j][4]);
		double dpbyp3 = dpbypnom / sqrt(evt[k][4]);
		double term1 = dpbyp1*dpbyp1 * (mass12*mass12 + mass31*mass31) * (mass12*mass12 + mass31*mass31);
		double term2 = dpbyp2*dpbyp2 * (mass12*mass12 + mass23*mass23) * (mass12*mass12 + mass23*mass23);
		double term3 = dpbyp3*dpbyp3 * (mass23*mass23 + mass31*mass31) * (mass23*mass23 + mass31*mass31);
		dmass = sqrt(0.25*(term1 + term2 + term3) / (mass*mass) + widthT*widthT);
//		dmass = sqrt(0.25*(dpbyp1*dpbyp1 + dpbyp2*dpbyp2 + dpbyp3*dpbyp3) * mass + widthT*widthT);
	}
	chisqT = (mass-massT)*(mass-massT) / (dmass*dmass);
	
	return chisqT;
}

//***************************************************************************/
void TopSearch::GetMHT(double MHT[]) {
// returns MHT into array MHT[]

	for (int i = 0; i < ndat; ++i) {
		MHT[i] = 0.;
	}
	
	MHT[0] = 10;
	for (int i = 0; i < npart; ++i) {
		MHT[1] -= evt[i][1];
		MHT[2] -= evt[i][2];
	}
	MHT[6] = sqrt(MHT[1]*MHT[1] + MHT[2]*MHT[2]);
	MHT[8] = atan2(MHT[1], MHT[2]);
	MHT[9] = 0.;
//	cout << " MHT = " << MHT[6] << endl;
	
	return;

}

//***************************************************************************/
int TopSearch::pzmissfrW (double lept[], double MET[], double MET1[], double MET2[]) {
// computes the solutions of pzmiss from W mass constraint
// returns the 2 solutions (if any) in MET1[] and MET2[]
	
	double pi = 3.141592654;
	int nsol;
	
	double pTl = sqrt(lept[1]*lept[1] + lept[2]*lept[2]);
	double pzl = lept[3];
	double pl = sqrt(pTl*pTl + pzl*pzl);
	double pTm = sqrt(MET[1]*MET[1] + MET[2]*MET[2]);
	double pTlpTm = lept[1]*MET[1] + lept[2]*MET[2];
	
	double A = massW*massW + 2.*pTlpTm;
	double B = A*A - 4.*pl*pl*pTm*pTm;
	
	double disc = A*A*pzl*pzl + (pl*pl-pzl*pzl)*B;
//	cout << " discr = " << disc << endl;
	if (disc < 0.) return 0;
	nsol = 2;
	double num1 = A * pzl;
	double num2 = sqrt(disc);
	double denom = pl*pl - pzl*pzl;
	if (denom <= 0.) return 0;
//	cout << " num1 = " << num1 << ", num2 = " << num2 << ", denom = " << denom << endl;
	double sol1 = 0.5 * (num1 - num2) / denom;
	double sol2 = 0.5 * (num1 + num2) / denom;
//	cout << " sol1 = " << sol1 << ", sol2 = " << sol2 << endl;
	for (int i = 0; i < ndat; ++i) {
		MET1[0] = MET[0];
		MET2[0] = MET[0];
		MET1[1] = MET[1];
		MET2[1] = MET[1];
		MET1[2] = MET[2];
		MET2[2] = MET[2];
		MET1[3] = sol1;
		MET2[3] = sol2;
		MET1[5] = massW;
		MET2[5] = massW;
		double pmom1 = sqrt(MET1[1]*MET1[1] + MET1[2]*MET1[2] + MET1[3]*MET1[3]);
		double pmom2 = sqrt(MET2[1]*MET2[1] + MET2[2]*MET2[2] + MET2[3]*MET2[3]);
		MET1[4] = sqrt(pmom1*pmom1 + MET1[5]*MET1[5]);
		MET2[4] = sqrt(pmom2*pmom2 + MET2[5]*MET2[5]);
		MET1[6] = MET[6];
		MET2[6] = MET[6];
		double thet1 = asin(MET1[6] / pmom1);
		if (thet1 < 0.) thet1 += pi;
		double thet2 = asin(MET2[6] / pmom2);
		if (thet2 < 0.) thet2 += pi;
		MET1[7] = -log(tan(0.5*thet1));
		MET2[7] = -log(tan(0.5*thet2));
		MET1[8] = MET[8];
		MET2[8] = MET[8];
		MET1[9] = MET[9];
		MET2[9] = MET[9];
	}
	
	return nsol;

}



