#include "Waveform.h"
#include "TXMLEngine.h"

void Init(void);
void Event(Int_t nPhe, Waveform* wf);
void LoadSimConfig();
Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err);
Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err);
Double_t fitfunc(Double_t *x, Double_t *par);

void DisplayNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level);
void InitConfig(std::string &linestr,std::string &strvalue);
void ParseConfig(TXMLEngine* xml, XMLNodePointer_t node);

void Init(void) {
   TXMLEngine* xml = new TXMLEngine;

   XMLDocPointer_t xmldoc = xml->ParseFile("./SimConfig.xml");
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

   gFSctime->SetParameter(0,ScintDecay);
   gFSingle->SetParameter(0,SPwidth);
   gFSingle->SetParameter(1,0);
   gFSingle->FixParameter(2,Gain);
   gAPtime ->FixParameter(0,APtimeconstant);
}

void wf_simulation(void) {
   Init();
   std::vector<Double_t> vecSum;
   std::vector<std::vector<Double_t>> vecNoise;
   std::vector<std::vector<Double_t>> vecInterference;

   TFile* fout = new TFile("fout.root","recreate");
   TTree* tout = new TTree("tout","tout");
   Int_t noiseNum = noiselist.size();
   Double_t *charge = new Double_t[noiseNum];
   Double_t *noisevar = new Double_t[noiseNum];
   Double_t *noiselevel = new Double_t[noiseNum];
   tout->Branch("noiselist" ,&noiseNum ,"noiselist/I");
   tout->Branch("charge"    ,charge    ,"charge[noiselist]/D");
   tout->Branch("noisevar"  ,noisevar  ,"noisevar[noiselist]/D");
   tout->Branch("noiselevel",noiselevel,"noiselevel[noiselist]/D");

   Double_t sum;

   Int_t npheRange[2] = {RangeMin,RangeMax};

   Double_t logstep = (TMath::Log(npheRange[1]) - TMath::Log(npheRange[0]))/ (Nstep -1);

   for (int istep = 0; istep < Nstep; istep++) {
      int iphe = istep*(npheRange[1]-npheRange[0])/Nstep;
      for(int irep=0;irep<Nevent;irep++){
         cout<<istep<<" "<<iphe<<endl;

         for (int i = 0; i < noiselist.size(); i++) {
            noiselevel[i]=noiselist[i];
            Waveform* wf= new Waveform(Nbins,-100,1000);
            wf->SetDarkNoiseFrequency(DNFreq);
            wf->MakeEvent(iphe);
            wf->MakeDarkNoise();
            wf->SetNoiseLevel(noiselevel[i]);
            wf->MakeElectricNoise();
            charge[i]  = wf->GetCharge();
            wf->Differentiate(Ndiff);
            noisevar[i] = wf->GetTotalVariance();

            Int_t index=istep*Nevent+irep;
            delete wf;
         }
         tout->Fill();
      }
   }
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
   if (linestr.find(strLambda)        !=string::npos)         lambda = std::stod(strvalue);
   if (linestr.find(strAlpha)         !=string::npos)          alpha = std::stod(strvalue);
   if (linestr.find(strScintDecay)    !=string::npos)     ScintDecay = std::stod(strvalue);
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
   if (linestr.find(strNoiseLevel)    !=string::npos)     noiselist.push_back(std::stod(strvalue));
}
