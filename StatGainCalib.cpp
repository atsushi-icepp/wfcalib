#include "Waveform.h"
#include "TXMLEngine.h"
// #include <boost/property_tree/ptree.hpp>

void Init(void);
void Event(Int_t nPhe, Waveform* wf);
// Double_t YCut(TGraph* gr);
// void YCut(TGraph* gr,Double_t *par);
void LoadSimConfig();

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

void StatGainCalib(TString file_argv="fout.root") {
   Init();
   std::vector<Double_t> vecSum;
   std::vector<std::vector<Double_t>> vecNoise;
   std::vector<std::vector<Double_t>> vecInterference;


   TFile* fout = new TFile(file_argv,"recreate");
   TTree* tout = new TTree("tout","tout");
   Int_t noiseNum = noiselist.size();
   Double_t *charge = new Double_t[noiseNum];
   Double_t *noisevar = new Double_t[noiseNum];
   Double_t *noiselevel = new Double_t[noiseNum];
   Int_t *Npho = new Int_t[noiseNum];
   Double_t NphoMean;
   tout->Branch("noiselist" ,&noiseNum ,"noiselist/I");
   tout->Branch("charge"    ,charge    ,"charge[noiselist]/D");
   tout->Branch("noisevar"  ,noisevar  ,"noisevar[noiselist]/D");
   tout->Branch("noiselevel",noiselevel,"noiselevel[noiselist]/D");
   tout->Branch("NphoMean"  ,&NphoMean ,"NphoMean/D");
   tout->Branch("Npho"      ,Npho      ,"Npho[noiselist]/I");

   Double_t sum;

   Double_t npheRange[2] = {RangeMin,RangeMax};

   Double_t logstep = (TMath::Log(npheRange[1]) - TMath::Log(npheRange[0]))/ (Nstep -1);

   for (int istep = 0; istep < Nstep; istep++) {
      // int iphe = (int)TMath::Exp(TMath::Log(npheRange[0]) + logstep * istep);
      Double_t iphe = TMath::Exp(logstep*istep)*npheRange[0];
      NphoMean=iphe;
      for(int irep=0;irep<Nevent;irep++){
         // cout<<istep<<" "<<iphe<<endl;
         for (int i = 0; i < noiselist.size(); i++) {

            noiselevel[i]=noiselist[i];
            if (LightSource=="LED") {
               Npho[i]=gRandom->Poisson(NphoMean);
            }else if(LightSource=="Scint"){
               Npho[i]=(Int_t)NphoMean;
            }
            Waveform* wf= new Waveform(Nbins,timemin,timemax);
            WaveformGen(wf,Npho[i],noiselevel[i],DNFreq);
            charge[i]  = wf->GetChargeIntegration(IntStart,IntEnd);
            wf->Differentiate(Ndiff);
            noisevar[i] = wf->GetTotalVariance(IntStart,IntEnd);

            Int_t index=istep*Nevent+irep;

            delete wf;
         } // end loop for noiselist.size
         tout->Fill();
      } // end loop for Nevent
   } // end loop for Nstep
   fout->cd();
   tout->Write();
   fout->Close();
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
   if (linestr.find(strRangeMin)      !=string::npos)       RangeMin = std::stod(strvalue);
   if (linestr.find(strRangeMax)      !=string::npos)       RangeMax = std::stod(strvalue);
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
   if (linestr.find(strPixelNoise)       !=string::npos)    PixelNoise  = std::stod(strvalue);
}

void WaveformGen(Waveform* wf,Int_t Npho,Double_t noiselevel,Int_t DNFreq){
   wf->MakeEvent(Npho);
   wf->SetDarkNoiseFrequency(DNFreq);
   wf->MakeDarkNoise();
   wf->SetNoiseLevel(noiselevel);
   wf->MakeElectricNoise();
}
