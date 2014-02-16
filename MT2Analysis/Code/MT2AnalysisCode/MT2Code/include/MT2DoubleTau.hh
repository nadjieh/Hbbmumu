#ifndef MT2DoubleTau_hh
#define MT2DoubleTau_hh


class MT2DoubleTau: public TObject {
public:
        MT2DoubleTau():tau0Ind(-1),tau1Ind(-1),chargeSum(-1),MT2(-1.),MT2Imbalanced(-1.),METImbalanced(-1.),METImbalancedPhi(-10.),
		       hasNoVetoElec(false),hasNoVetoMu(false),signalDoubleTau(false),isTau0NonIso(false),isTau1NonIso(false),QCDCategory(-1){};
        ~MT2DoubleTau(){};
        //Reset values to default
        void Reset(){
                tau0Ind = -1;
                tau1Ind = -1;
                chargeSum = -1;
                MT2 = -1.;
                MT2Imbalanced = -1.;
		METImbalanced = -1.;
		METImbalancedPhi = -10.;
                hasNoVetoElec = false;
                hasNoVetoMu = false;
		signalDoubleTau = false;
                isTau0NonIso = false;
                isTau1NonIso = false;
		QCDCategory = -1;
        };

        //Setters
        void SetTauIndex0(int input){ tau0Ind = input; };
        void SetTauIndex1(int input){ tau1Ind = input; };
        void SetSumCharge(int input){ chargeSum = input; };
        void SetMT2(float input){ MT2 = input; };
        void SetMT2Imbalanced(float input){ MT2Imbalanced = input; };
	void SetMETImbalanced( float input){ METImbalanced = input;};
	void SetMETImbalancedPhi( float input){ METImbalancedPhi = input;};
        void SetElecVeto(bool input){ hasNoVetoElec = input; };
        void SetMuVeto(bool input){ hasNoVetoMu = input; };
	void SetBeingSignal(bool input) { signalDoubleTau = input;};
        void SetTau0NonIso(bool input){ isTau0NonIso = input; };
        void SetTau1NonIso(bool input){ isTau1NonIso = input; };
	void SetQCDCategory( int input){ QCDCategory = input;}
        //Getters
        int GetTauIndex0(){ return tau0Ind; };
        int GetTauIndex1(){ return tau1Ind; };
        int GetSumCharge(){ return chargeSum; };
        float GetMT2(){ return MT2; };
        float GetMT2Imbalanced(){ return MT2Imbalanced; };
	float GetMETImbalanced(){ return METImbalanced;};
	float GetMETImbalancedPhi(){ return METImbalancedPhi;};
        bool HasNoVetoElec(){ return hasNoVetoElec; };
        bool HasNoVetoMu(){ return hasNoVetoMu; };
        bool IsTau0NonIso(){ return isTau0NonIso; };
        bool IsTau1NonIso(){ return isTau1NonIso; };

        //Event Selection Methods (only focus on leptons)
        bool isSignalDoubleTau(){
                /*bool ret = (tau0Ind != -1);
                ret = ret && (tau1Ind != -1);
                ret = ret && (chargeSum == 0);
                ret = ret && hasNoVetoElec;
                ret = ret && hasNoVetoMu;
                return ret;*/
		return signalDoubleTau;
        }

        int GetQCDCategory(){
                /*
		 * Not suitable for QCD: -1
                 * SS-Iso: 1
                 * SS-anti-Iso: 2
                 * OS-anti-Iso: 3
                */
		return QCDCategory;
                /*int ret = -1;
		if(tau0Ind == -1 || tau1Ind == -1) //QCD event should have valid tau indecies
					return ret;
		if(!hasNoVetoElec || !hasNoVetoMu) //QCD event is vetoed against additional e/mu
					return ret;
                if(chargeSum != 0 && !isTau0NonIso && !isTau1NonIso)
                        ret = 1;
                if(chargeSum != 0 && isTau0NonIso && isTau1NonIso)
                        ret = 2;
                if(chargeSum == 0 && !isTau0NonIso && !isTau1NonIso)
                        ret = 3;
                return ret;*/
        }

        void printObject(){
                cout<<"TauTau quantities in object: "<<endl;
                cout<<"\ttau0Ind: "<<tau0Ind<<endl;
                cout<<"\ttau1Ind: "<<tau1Ind<<endl;
                cout<<"\tchargeSum: "<<chargeSum<<endl;
                cout<<"\tMT2: "<<MT2<<endl;
                cout<<"\tMT2Imbalanced: "<<MT2Imbalanced<<endl;
                cout<<"\tMETImbalanced: "<<METImbalanced<<endl;
                cout<<"\tMETImbalancedPhi: "<<METImbalancedPhi<<endl;
                cout<<"\thasNoVetoElec: "<<hasNoVetoElec<<endl;
                cout<<"\thasNoVetoMu: "<<hasNoVetoMu<<endl;
                cout<<"\tsignalDoubleTau: "<<signalDoubleTau<<endl;
                cout<<"\tQCDCategory: "<<QCDCategory<<endl;
                cout<<"\t\tisTau0NonIso: "<<isTau0NonIso<<endl;
                cout<<"\t\tisTau1NonIso: "<<isTau1NonIso<<endl;
        }

private:
        //General properties
        int tau0Ind;
        int tau1Ind;
        int chargeSum;
        float MT2;
        float MT2Imbalanced;
	float METImbalanced;
        float METImbalancedPhi;
        bool hasNoVetoElec;
        bool hasNoVetoMu;
	//Check for being Signal
	bool signalDoubleTau;
        //Tau propertis for QCD backgrund estmation
        bool isTau0NonIso;
        bool isTau1NonIso;
	int  QCDCategory;
        ClassDef(MT2DoubleTau, 1)

};
#endif /*MT2DoubleTau_hh*/
                         
