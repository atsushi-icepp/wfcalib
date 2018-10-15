Double_t timemin = -400;
Double_t timemax = 700;
Double_t lambda;
string   strLambda="Lambda";
Double_t alpha;
string   strAlpha="Alpha";
Double_t ScintDecay;
string   strScintDecay="ScintDecay";
string   LightSource;
string   strLightSource="LightSource";
Double_t LEDWidth;
string   strLEDWidth  ="LEDWidth";
Double_t SPwidth;
string   strSPwidth="SinglePulseWidth";
Double_t APtimeconstant;
string   strAPtimeconstant="AfterPulseTimeConstant";
Double_t Gain=1;
string   strGain="Gain";
Double_t RangeMin;
string   strRangeMin="RangeMin";
Double_t RangeMax;
string   strRangeMax="RangeMax";
Int_t    Nstep;
string   strNstep="Nstep";
Int_t    Ndiff;
string   strNdiff="Ndiff";
Int_t    Nevent;
string   strNevent  ="Nevent";
Double_t    DNFreq;
string   strDNFreq  ="DarkNoiseFrequency";
Int_t    Nbins;
string   strNbins  ="Nbins";
Double_t BaselineStart;
string   strBLstart  ="BLstart";
Double_t BaselineEnd;
string   strBLend  ="BLend";
Double_t IntStart;
string   strIntStart  ="IntStart";
Double_t IntEnd;
string   strIntEnd  ="IntEnd";
// Int_t    Nevent;
std::vector<Double_t> noiselist;
string   strNoiseLevel  ="NoiseLevel";
Double_t PixelNoise;
string   strPixelNoise ="PixelNoise";
TString safeName(TString name);
Double_t funcsctime(Double_t *x, Double_t *par);
Double_t funcledtime(Double_t *x, Double_t *par);
Double_t funcsingle(Double_t *x, Double_t *par);
Double_t APtiming(Double_t *x, Double_t *par);
Double_t Borel(Double_t,Int_t); // caluculate the P(n) of borel distribution
Int_t Borel_gen(Double_t); // generate the random variable obeying Borel
TF1 *gFSingle = new TF1("fsingle", funcsingle, timemin, timemax, 3);
TF1 *gFSctime = new TF1("fsctime", funcsctime, timemin, timemax, 1);
TF1 *gFLEDtime = new TF1("fledtime", funcledtime, timemin, timemax, 1);
TF1 *gAPtime = new TF1("APtime", APtiming, 0, timemax, 1);

TString safeName(TString name){
	TObject* old = gROOT->FindObject(name);
	if(old) delete old;
	return name;
}

Double_t funcsctime(Double_t *x, Double_t *par) {
	//par 0: decay time const
	if (x[0] < 0) return 0;
	else return TMath::Exp(-x[0]/par[0]);
}

Double_t funcledtime(Double_t *x, Double_t *par) {
	//par 0: decay time const
	if (x[0] < 0) return 0;
	else if(x[0] > par[0]) return 0;
	else return 1;
}

Double_t funcsingle(Double_t *x, Double_t *par) {
	//par 0: decay time const
	//par 1: start time of pulse
	//par 2: scale of the single pulse
	//Area normalized
	if (x[0] < par[1]) return 0;
	else if (x[0] > par[1] + par[0] * 4) return 0;
	else return par[2]*TMath::Exp(-(x[0] - par[1])/par[0]) / (par[0]*(1-TMath::Exp(-4)));
}

Double_t APtiming(Double_t *x, Double_t *par) {
	//par 0: time const for AP timing distribution (assuming exponential dist)
	if (x[0] < 0) return 0;
	else return TMath::Exp(-x[0]/par[0]);
}

Double_t Borel(Double_t lambda, Int_t n){
	Double_t exp = TMath::Exp(-lambda*n);
	Double_t power = 1.;
	for (Int_t i=1;i<n;++i){
		power *= lambda*n;
	}
	Double_t numerator = TMath::Factorial(n);
	return exp*power/numerator;
}

Int_t Borel_gen(Double_t lambda){
	Int_t return_rn = 0;
	Double_t rand_uni = gRandom->Uniform();
	while(rand_uni > 0){
		return_rn++;
		rand_uni -= Borel(lambda,return_rn);
	}
	return return_rn;
}


class Waveform{
protected:
	Int_t           fNPoints;
	Double_t        fBinsize;
	Double_t        *fNoiseAmp;
	Double_t        *fAmplitude;
	Double_t        *fTime;
	Double_t        fTimeMin=timemin;      //Time of the first (0) point used for fix bin
	Double_t        fTimeMax=timemax;
	Double_t        fDNFreq=0;
	Double_t        fNoiseLevel=0;
	Int_t           fNDN;
public:
	// Double_t variance;
	Waveform(Int_t npoints,Double_t timemin,Double_t timemax);
	~Waveform();
	void SetAmplitude(Double_t* amplitude);
	void SetAmplitudeAt(Int_t ibin,Double_t amp);
	void Differentiate(Int_t ndiff);
	void MakeEvent(Int_t npe);
	void MakeDarkNoise();
	void MakeElectricNoise();
	void SetNoiseLevel(Double_t noiselevel){fNoiseLevel=noiselevel;};
	void SetDarkNoiseFrequency(Double_t DNfreq){fDNFreq=DNfreq;};
	void Draw();
	Double_t GetAmplitude(Int_t ibin);
	Double_t GetTimeAt(Int_t ibin);
	Double_t GetTotalVariance(Double_t IntStart,Double_t IntEnd);
	Double_t GetSignalVariance(Double_t IntStart,Double_t IntEnd);
	Double_t GetChargeIntegration(Double_t IntStart,Double_t IntEnd);
	Double_t GetPulseHeight(Double_t IntStart,Double_t IntEnd);
	Double_t GetBaseLineVariance(Double_t BaselineStart,Double_t BaselineEnd);
	Double_t GetRMS();
	Double_t GetNoiseConst();
	Double_t GetInterference();
	Int_t    GetNDN(){return fNDN;};
	Int_t    FindPoint(Double_t x);
	Int_t    GetRegionNpoints(Double_t Start,Double_t End);
	// Int_t    GetSignalNpoints(Double_t IntStart,Double_t IntEnd);

};

Waveform::Waveform(Int_t npoints,Double_t timemin,Double_t timemax){
	fNPoints   = npoints;
	fTime      = new Double_t[npoints];
	fAmplitude = new Double_t[npoints];
	fNoiseAmp  = new Double_t[npoints];
	// fTimeMin   = timemin;
	// fTimeMax   = timemax;
	memset(fAmplitude, 0, sizeof(Double_t)*npoints);
	memset(fNoiseAmp,  0, sizeof(Double_t)*npoints);
	for (int ipnt = 0; ipnt < npoints; ipnt++) {
		fTime[ipnt]=timemin+(timemax-timemin)*ipnt/fNPoints;
	}
	fBinsize= (fTimeMax-fTimeMin)/(Double_t)fNPoints;
	// std::cout<<"fbinsize : "<<fBinsize<<" fNPoints: "<<fNPoints<<std::endl;
}

Waveform::~Waveform(){
	// fNPoints   = npoints;
	delete[] fTime;
	delete[] fAmplitude;
	delete[] fNoiseAmp;
}

void Waveform::SetAmplitudeAt(Int_t ibin,Double_t amp){
	if (ibin>0&&ibin<fNPoints) {
		fAmplitude[ibin]=amp;
	}
}

Double_t Waveform::GetAmplitude(Int_t ibin){
	// if (ibin>0&&ibin<fNPoints) {
	return   fAmplitude[ibin];
	// }
}

Double_t Waveform::GetTimeAt(Int_t ibin){
	// if (ibin>0&&ibin<fNPoints) {
	return   fTime[ibin];
	// }
}

void Waveform::Draw(){
	TCanvas* canvastemp= new TCanvas(safeName("canvastemp"),safeName("canvastemp"),600,600);
	TGraph* grwf= new TGraph();
	for (Int_t ipoint = 0; ipoint < fNPoints; ipoint++) {
		grwf->SetPoint(ipoint,fTime[ipoint],fAmplitude[ipoint]);
	}
	grwf->Draw("ap");

}

void Waveform::Differentiate(Int_t ndiff){

	if (ndiff>0) {
		Double_t* tempWF= new Double_t[fNPoints];
		Double_t* diffWF= new Double_t[fNPoints];
		Double_t* tempNWF= new Double_t[fNPoints];
		Double_t* diffNWF= new Double_t[fNPoints];
		memcpy(tempWF,fAmplitude,fNPoints*sizeof(Double_t) );
		memcpy(tempNWF,fNoiseAmp,fNPoints*sizeof(Double_t) );

		for(int idiff=0;idiff<ndiff;idiff++){

			for (int iBin = 1; iBin < fNPoints-1; iBin++) {
				diffWF[iBin]=(tempWF[iBin+1]-tempWF[iBin-1]);
				diffNWF[iBin]=(tempNWF[iBin+1]-tempNWF[iBin-1]);
			}
			diffWF[0]=0;
			diffWF[fNPoints-1]=0;
			diffNWF[0]=0;
			diffNWF[fNPoints-1]=0;
			memcpy(tempWF,diffWF,fNPoints*sizeof(Double_t) );
			memcpy(tempNWF,diffNWF,fNPoints*sizeof(Double_t) );
		}
		memcpy(fAmplitude,tempWF,fNPoints*sizeof(Double_t) );
		memcpy(fNoiseAmp,tempNWF,fNPoints*sizeof(Double_t) );
	}
}

void Waveform::MakeElectricNoise(){
	for (int ipnt = 0; ipnt < fNPoints; ipnt++) {
		Double_t noise=gRandom->Gaus(0,fNoiseLevel);
		// fAmplitude[ipnt]+=noise;
		fNoiseAmp[ipnt]=noise;
	}
}

void Waveform::SetAmplitude(Double_t* amplitude)
{
	// Copy amplitude
	if (amplitude != fAmplitude)
	memcpy(fAmplitude, amplitude, sizeof(Double_t)*fNPoints);
}

void Waveform::MakeEvent(Int_t npe){ // define npe as initial number of photon
	gRandom->SetSeed(0);
	for (size_t inpe = 0; inpe < npe; inpe++) {
		Double_t pulsetime;
		if (LightSource=="Scint") {
			pulsetime = gFSctime->GetRandom();
		}else if(LightSource=="LED"){
			pulsetime = gFLEDtime->GetRandom();
		}

		Int_t CT_num = Borel_gen(lambda);
                Double_t fluctuate_pix = gRandom->Gaus(1.,PixelNoise/TMath::Sqrt(CT_num));
		// Double_t pulsetime=gRandom->Uniform(0,30);
		gFSingle->SetParameter(1, pulsetime);

		for (int iBin = 0; iBin < Nbins; iBin++) {
			fAmplitude[iBin]+= fluctuate_pix*CT_num*gFSingle->Eval(fTime[iBin]);
		}
		// complete pulse generation with cross talk: After Pulse below
		for (Int_t ap_cand=0; ap_cand<CT_num; ap_cand++){
			Double_t uniform_rn = gRandom->Uniform();
			if (uniform_rn < alpha){
				Double_t APtime = gAPtime->GetRandom();
				Double_t decayconst = gFSingle->GetParameter(0);
				Double_t AP_amp = 1.-TMath::Exp(-APtime/decayconst);
				gFSingle->SetParameter(1,pulsetime+APtime);
				for(int iBin=0; iBin < Nbins; iBin++){
					fAmplitude[iBin]+= AP_amp*gFSingle->Eval(fTime[iBin]);
				}
			}
		}
	} // end of loop considering npe
}

void Waveform::MakeDarkNoise(){
	//generate Dark Noise from -500 ~ 1500

	Double_t Nexp  = fDNFreq*2000/TMath::Power(10,9);
	fNDN=gRandom->Poisson(Nexp);
	// std::cout<<"N DN: "<<fNDN;
	// Double_t ProbinWindow =
	for (Int_t iDN = 0; iDN < fNDN; iDN++) {
		Double_t DNtime=gRandom->Uniform(-500,1500);
		// std::cout<<" time : "<<DNtime<<std::endl;
		Int_t CT_num = Borel_gen(lambda);
		Double_t fluctuate_pix = gRandom->Gaus(1.,PixelNoise/TMath::Sqrt(CT_num));
		gFSingle->SetParameter(1,DNtime);
		for(int iBin=0; iBin < Nbins; iBin++){
			fAmplitude[iBin]+= fluctuate_pix*CT_num*gFSingle->Eval(fTime[iBin]);
		}
                for (Int_t ap_cand=0; ap_cand<CT_num; ap_cand++){
                        Double_t uniform_rn = gRandom->Uniform();
                        if (uniform_rn < alpha){
                                Double_t APtime = gAPtime->GetRandom();
                                Double_t decayconst = gFSingle->GetParameter(0);
                                Double_t AP_amp = 1.-TMath::Exp(-APtime/decayconst);
                                gFSingle->SetParameter(1,DNtime+APtime);
                                for(int iBin=0; iBin < Nbins; iBin++){
                                        fAmplitude[iBin]+= AP_amp*gFSingle->Eval(fTime[iBin]);
                                }
                        }
                }// end of afterpulse 

	}// end loop for dark events



}

Int_t Waveform::FindPoint(Double_t x)
{
   // Find point which is maxmum one under x
   // -1 : underflow
   //  0 : first point
   // ,,,
   // fNPoints : last bin (overflow)
   //
   // cell#     0   1   2   3   4           fNPoints-1
   //         --*---*---*---*---*------*---*---*--
   // return -1   0   1   2   3   4             fNPoints-1
   //

   Int_t point;
   // if (fFixBinSize) {        // fix bins
      if (x < fTimeMin) point = -1;
      else  if (x > fTimeMax) point = fNPoints-1;
      else point = static_cast<int>((fNPoints-1) * (x-fTimeMin) / (fTimeMax-fTimeMin));
   // } else {         //Binary search
   //    if (x < fTime[startpoint]) point = -1;
   //    else point = static_cast<int>(TMath::BinarySearch(fNPoints - startpoint, fTime + startpoint, x)) + startpoint;
   // }
   return point;
}

Double_t Waveform::GetChargeIntegration(Double_t IntStart,Double_t IntEnd){
	Double_t charge=0;
	Int_t    startpoint = FindPoint(IntStart) + 1;
   Int_t    endpoint   = FindPoint(IntEnd) + 1;
	for (int ipnt = startpoint; ipnt < endpoint; ipnt++) {
		charge+=(fAmplitude[ipnt]+fNoiseAmp[ipnt])*fBinsize;
		// std::cout<<"bin: "<<fBinsize<<" amp: "<<fAmplitude[ipnt]<<std::endl;
	}
	return charge;
}

Double_t Waveform::GetPulseHeight(Double_t IntStart,Double_t IntEnd){
	Double_t height=-1e9;
        Int_t startpoint = FindPoint(IntStart)+1;
	Int_t endpoint = FindPoint(IntEnd)+1;
	for(int ipnt=startpoint;ipnt<endpoint;ipnt++){
		if(fAmplitude[ipnt]+fNoiseAmp[ipnt]>height) height=fAmplitude[ipnt]+fNoiseAmp[ipnt];
	}
	return height;
}

Double_t Waveform::GetBaseLineVariance(Double_t BaselineStart,Double_t BaselineEnd){
	Double_t variance=0;
	Int_t    startpoint = FindPoint(BaselineStart) + 1;
   Int_t    endpoint   = FindPoint(BaselineEnd) + 1;
	for (int ipnt = startpoint; ipnt < endpoint; ipnt++) {
		variance+=TMath::Power((fAmplitude[ipnt]+fNoiseAmp[ipnt])*fBinsize,2);
	}
	return variance;
}

Double_t Waveform::GetTotalVariance(Double_t IntStart,Double_t IntEnd){
	Double_t variance=0;
	Int_t    startpoint = FindPoint(IntStart) + 1;
   Int_t    endpoint   = FindPoint(IntEnd) + 1;
	for (int ipnt = startpoint; ipnt < endpoint; ipnt++) {
		variance+=TMath::Power((fAmplitude[ipnt]+fNoiseAmp[ipnt])*fBinsize,2);
	}
	return variance;
}

Double_t Waveform::GetSignalVariance(Double_t IntStart,Double_t IntEnd){
	Double_t variance=0;
	Int_t    startpoint = FindPoint(IntStart) + 1;
   Int_t    endpoint   = FindPoint(IntEnd) + 1;
	for (int ipnt = startpoint; ipnt < endpoint; ipnt++) {
		variance+=TMath::Power(fAmplitude[ipnt]*fBinsize,2);
	}
	return variance;
}

Double_t Waveform::GetInterference(){
	Double_t interference=0;
	for (int ipnt = 0; ipnt < fNPoints; ipnt++) {
		interference += 2 * fAmplitude[ipnt] * fNoiseAmp[ipnt] * TMath::Power(fBinsize,2);
	}
	return interference;
}


Double_t Waveform::GetNoiseConst(){
	Double_t noiseconst=0;
	for (int ipnt = 0; ipnt < fNPoints; ipnt++) {
		noiseconst += fNoiseAmp[ipnt] * fNoiseAmp[ipnt] * TMath::Power(fBinsize,2);
	}
	return noiseconst;
}

Int_t    Waveform::GetRegionNpoints(Double_t Start,Double_t End){
return FindPoint(End)-FindPoint(Start);
}
// Int_t    GetSignalNpoints(Double_t IntStart,Double_t IntEnd){
// return FindPoint(IntEnd)-FindPoint(IntStart);
// }
