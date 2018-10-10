#include "Waveform.h"
#include "TXMLEngine.h"
// #include <boost/property_tree/ptree.hpp>

void Init(void);
void Event(Int_t nPhe, Waveform* wf);
// Double_t YCut(TGraph* gr);
// void YCut(TGraph* gr,Double_t *par);
void LoadSimConfig();
Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err);
Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err);
// void Differentiate(Double_t *wf,Int_t ndiff);
Double_t fitfunc(Double_t *x, Double_t *par);

void DisplayNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level);
void InitConfig(std::string &linestr,std::string &strvalue);
void ParseConfig(TXMLEngine* xml, XMLNodePointer_t node);
void WaveformGen(Waveform* wf,Int_t Npho,Double_t noiselevel,Int_t DNFreq);
// TFile* fResult=new TFile("Result.root","recreate");

void Init(void) {
   TXMLEngine* xml = new TXMLEngine;

   XMLDocPointer_t xmldoc = xml->ParseFile("./StatConfig.xml");
   if (xmldoc==0) {
      delete xml;
      return;
   }
   XMLNodePointer_t mainnode = xml->DocGetRootElement(xmldoc);
   // display recursively all nodes and subnodes
   // DisplayNode(xml, mainnode, 1);
   ParseConfig(xml, mainnode);
   // Release memory before exit
   xml->FreeDoc(xmldoc);
   delete xml;

   // LoadSimConfig();
   // std::cout<<"ScintDecay: "<<ScintDecay<<std::endl;
   gFSctime  ->SetParameter(0,ScintDecay);
   gFLEDtime ->SetParameter(0,LEDWidth);
   gFSingle  ->SetParameter(0,SPwidth);
   gFSingle  ->SetParameter(1,0);
   gFSingle  ->FixParameter(2,Gain);
   gAPtime   ->FixParameter(0,APtimeconstant);
}

void StatGainCalib(void) {
   // TCanvas *cgr1 =new TCanvas("cgr1", "cgr1",600,600);
   // TCanvas *cgrQNvar =new TCanvas("cgrQNvar", "cgrQNvar",600,600);
   Init();
   std::vector<Double_t> vecSum;
   std::vector<std::vector<Double_t>> vecNoise;
   std::vector<std::vector<Double_t>> vecInterference;


   TFile* fout = new TFile("fout.root","recreate");
   TTree* tout = new TTree("tout","tout");
   Double_t charge;
   Int_t    Npho;
   Int_t    NphoMean;
   Double_t noisevar;
   Double_t noiselevel;
   tout->Branch("charge"    ,&charge    ,"charge/D");
   tout->Branch("NphoMean"  ,&NphoMean  ,"NphoMean/I");
   tout->Branch("Npho"      ,&Npho      ,"Npho/I");
   tout->Branch("noisevar"  ,&noisevar  ,"noisevar/D");
   tout->Branch("noiselevel",&noiselevel,"noiselevel/D");

   Double_t sum;

   Int_t npheRange[2] = {RangeMin,RangeMax};

   Double_t logstep = (TMath::Log(npheRange[1]) - TMath::Log(npheRange[0]))/ (Nstep -1);

   TClonesArray* cagrQNvar = new TClonesArray("TGraphErrors",noiselist.size());
   for (int i = 0; i < noiselist.size(); i++) {
      new ((TGraphErrors*)(*cagrQNvar)[i]) TGraphErrors();
   }


   for (int istep = 0; istep < Nstep; istep++) {
      // int iphe = (int)TMath::Exp(TMath::Log(npheRange[0]) + logstep * istep);
      int iphe = (npheRange[1]-npheRange[0])*istep/Nstep+npheRange[0];
      for(int irep=0;irep<Nevent;irep++){
         // cout<<istep<<" "<<iphe<<endl;
         for (int i = 0; i < noiselist.size(); i++) {

            noiselevel=noiselist[i];
            NphoMean=iphe;
            if (LightSource=="LED") {
               Npho=gRandom->Poisson(NphoMean);
            }else if(LightSource=="Scint"){
               Npho=NphoMean;
            }
            Waveform* wf= new Waveform(Nbins,timemin,timemax);
            WaveformGen(wf,Npho,noiselevel,DNFreq);
            charge  = wf->GetChargeIntegration(IntStart,IntEnd);
            // cout<<istep<<" "<<Npho<<"charge: "<<charge<<" blvar: "<<blvar<<endl;
            wf->Differentiate(Ndiff);
            noisevar = wf->GetTotalVariance(IntStart,IntEnd);
            // Double_t blvar   = wf->GetBaseLineVariance(BaselineStart,BaselineEnd);
            // Double_t BLratio = (Double_t)wf->GetRegionNpoints(BaselineStart,BaselineEnd)/
            // (Double_t)wf->GetRegionNpoints(IntStart,IntEnd);
            // noisevar-=blvar/BLratio;
            // cout<<istep<<" "<<Npho<<"charge: "<<charge<<" blvar: "<<blvar<<" BLratio: "<<BLratio<<endl;

            Int_t index=istep*Nevent+irep;
            ((TGraphErrors*)(*cagrQNvar)[i])->SetPoint(     index,charge,noisevar);
            ((TGraphErrors*)(*cagrQNvar)[i])->SetPointError(index,0     ,noisevar*0.2);
            // if (wf->GetNDN()>0) {
            //    // std::cout<<"NDN: "<<wf->GetNDN()<<std::endl;
            //    wf->Draw();
            // }

            delete wf;
            tout->Fill();
         } // end loop for noiselist.size
      } // end loop for Nevent
   } // end loop for Nstep
   fout->cd();
   tout->Write();
   fout->Close();
   TString fname="gain1_noise5e-2.pdf";

   TF1* func= new TF1("fitfunc",fitfunc,-10,5000,3);

   for (int i = 0; i < noiselist.size(); i++) {

      Double_t par[3];
      Double_t err[3];
      TF1* pol1= new TF1("fitfunc",fitfunc,-100,5000,3);
      // pol1->SetParLimits(2,0,1);
      // pol1->SetParLimits(0,-1,20);
      ((TGraphErrors*)(*cagrQNvar)[i])->Fit("fitfunc","","");
      // Double_t par[3];
      for (int i = 0; i < 3; i++) {
         par[i]=pol1->GetParameter(i);
         err[i]=pol1->GetParError(i);
      }

      Double_t MesGain=par[1]/(1-par[1]*par[2]);
      Double_t pmul=par[1]*par[2];
      Double_t pmulerr=MultipleError(par[1],err[1],par[2],err[2]);
      Double_t MesGainerr=DivisionError(par[1],err[1],(double)1-pmul,pmulerr);
      std::cout<<"Gain: "<<MesGain<<"+-"<<MesGainerr<<std::endl;

      std::cout<<"Pmul: "<<pmul<<"+-"<<pmulerr<<std::endl;
      // ((TGraph*)(*cagrQNvar)[i])->Fit("fitfunc","","");
      ((TGraphErrors*)(*cagrQNvar)[i])->SetTitle("Q_{dint}vs Q^{2}_{drms};Q_{dint};Q^{2}_{drms}");
      ((TGraph*)(*cagrQNvar)[i])->SetMaximum(500);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMinimum(0);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMarkerStyle(20);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMarkerColor(2+i);
   }
}

Double_t fitfunc(Double_t *x, Double_t *par){
   //par[0]: offset;
   //par[1]:
   Double_t  xx        = x[0];
   // Double_t PMultiple = TMath::Sqrt(TMath::Power(par[2]*xx,2)/1+TMath::Power(par[2]*xx,2));
   // Double_t PMultiple = par[2];
   // Double_t Gain = 1;
   return par[0]+par[1]*(xx+par[2]*xx*xx);
}

Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err){
   /*calculate error of value 1/ value 2*/
   return TMath::Sqrt(TMath::Power(value1err/value2,2)+TMath::Power(value1*value2err/(value2*value2),2));
}

Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err){
   /*calculate error of value 1* value 2*/
   return TMath::Sqrt(TMath::Power(value1err* value2,2)+TMath::Power(value1*value2err,2));
}

void DisplayNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level)
{
   // this function display all accessible information about xml node and its children
   printf("%*c node: %s\n",level,' ', xml->GetNodeName(node));
   // display namespace
   XMLNsPointer_t ns = xml->GetNS(node);
   if (ns!=0)
   printf("%*c namespace: %s refer: %s\n",level+2,' ', xml->GetNSName(ns), xml->GetNSReference(ns));
   // display attributes
   XMLAttrPointer_t attr = xml->GetFirstAttr(node);
   while (attr!=0) {
      printf("%*c attr: %s value: %s\n",level+2,' ', xml->GetAttrName(attr), xml->GetAttrValue(attr));
      attr = xml->GetNextAttr(attr);
   }
   // display content (if exists)
   const char* content = xml->GetNodeContent(node);
   if (content!=0)
   printf("%*c cont: %s\n",level+2,' ', content);
   // display all child nodes
   XMLNodePointer_t child = xml->GetChild(node);
   while (child!=0) {
      DisplayNode(xml, child, level+2);
      child = xml->GetNext(child);
   }
}

void ParseConfig(TXMLEngine* xml, XMLNodePointer_t node){
   // this function display all accessible information about xml node and its children
   std::string NodeName = xml->GetNodeName(node);
   std::string ConfigName;
   std::string ConfigValue;
   if (NodeName=="config") {
      XMLNodePointer_t childNode = xml->GetChild(node);
      // std::cout<<"n: "<<NodeName<<std::endl;
      while (childNode!=0) {
         std::string ChildNodeName = xml->GetNodeName(childNode);
         // std::cout<<"cn: "<<ChildNodeName<<std::endl;
         if (ChildNodeName=="name") {
            ConfigName = xml->GetNodeContent(childNode);

         }else if(ChildNodeName=="values"){
            XMLNodePointer_t valueNode = xml->GetChild(childNode);
            while (valueNode!=0) {
               ConfigValue = xml->GetNodeContent(valueNode);
               valueNode = xml->GetNext(valueNode);
               // std::cout<<"ConfigName: "<<ConfigName<<std::endl;
               // std::cout<<"ConfigValue: "<<ConfigValue<<std::endl;
               InitConfig(ConfigName,ConfigValue);
            }

         }

         childNode = xml->GetNext(childNode);
      }


   }

   XMLNodePointer_t child = xml->GetChild(node);
   while (child!=0) {
      ParseConfig(xml, child);
      child = xml->GetNext(child);
   }
}

void InitConfig(std::string &linestr,std::string &strvalue){
   if (linestr.find(strLightSource)   !=string::npos)    LightSource = strvalue;
   if (linestr.find(strScintDecay)    !=string::npos)     ScintDecay = std::stod(strvalue);
   if (linestr.find(strLEDWidth)      !=string::npos)     LEDWidth   = std::stod(strvalue);
   if (linestr.find(strLambda)        !=string::npos)         lambda = std::stod(strvalue);
   if (linestr.find(strAlpha)         !=string::npos)          alpha = std::stod(strvalue);
   if (linestr.find(strSPwidth)       !=string::npos)        SPwidth = std::stod(strvalue);
   if (linestr.find(strAPtimeconstant)!=string::npos) APtimeconstant = std::stod(strvalue);
   if (linestr.find(strGain)          !=string::npos)           Gain = std::stod(strvalue);
   if (linestr.find(strRangeMin)      !=string::npos)       RangeMin = std::stoi(strvalue);
   if (linestr.find(strRangeMax)      !=string::npos)       RangeMax = std::stoi(strvalue);
   if (linestr.find(strNstep)         !=string::npos)          Nstep = std::stoi(strvalue);
   if (linestr.find(strNdiff)         !=string::npos)          Ndiff = std::stoi(strvalue);
   if (linestr.find(strNevent)        !=string::npos)         Nevent = std::stoi(strvalue);
   if (linestr.find(strDNFreq)        !=string::npos)         DNFreq = std::stod(strvalue);
   if (linestr.find(strNbins)         !=string::npos)         Nbins  = std::stoi(strvalue);
   if (linestr.find(strBLstart)       !=string::npos)  BaselineStart = std::stod(strvalue);
   if (linestr.find(strBLend)         !=string::npos)   BaselineEnd  = std::stod(strvalue);
   if (linestr.find(strIntStart)      !=string::npos)       IntStart = std::stod(strvalue);
   if (linestr.find(strIntEnd)        !=string::npos)        IntEnd  = std::stod(strvalue);
   if (linestr.find(strNoiseLevel)    !=string::npos)     noiselist.push_back(std::stod(strvalue));
   if (linestr.find(PixelNoise)       !=string::npos)    PixelNoise  = std::stod(strvalue);
}

void WaveformGen(Waveform* wf,Int_t Npho,Double_t noiselevel,Int_t DNFreq){
   wf->MakeEvent(Npho);
   wf->SetDarkNoiseFrequency(DNFreq);
   wf->MakeDarkNoise();
   wf->SetNoiseLevel(noiselevel);
   wf->MakeElectricNoise();
}
