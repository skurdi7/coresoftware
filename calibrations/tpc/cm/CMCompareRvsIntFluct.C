#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"

using namespace std;

class Shifter {
public:
  Shifter(TString sourcefilename, TString averagefilename);
  TFile *forward, *average;
  TH3F *hPosX, *hPosY, *hPosZ, *hPosR, *hPosPhi, *hPosXave, *hPosYave, *hPosZave, *hPosRave, *hPosPhiave;
  TH3F *hNegX, *hNegY, *hNegZ, *hNegR, *hNegPhi, *hNegXave, *hNegYave, *hNegZave, *hNegRave, *hNegPhiave;
  TH3F *hXBack, *hYBack, *hZBack;  
};

Shifter::Shifter(TString sourcefilename, TString averagefilename){
  //single event distortion file
  forward=TFile::Open(sourcefilename,"READ"); 

  //positive side (z = 0 to 105.5)
  hPosX=(TH3F*)forward->Get("hIntDistortionPosX");
  hPosY=(TH3F*)forward->Get("hIntDistortionPosY");
  hPosZ=(TH3F*)forward->Get("hIntDistortionPosZ");

  hPosR=(TH3F*)forward->Get("hIntDistortionPosR");
  hPosPhi=(TH3F*)forward->Get("hIntDistortionPosP");

  hPosX=(TH3F*)forward->Get("hIntDistortionPosX");
  hPosY=(TH3F*)forward->Get("hIntDistortionPosY");
  hPosZ=(TH3F*)forward->Get("hIntDistortionPosZ");

  hPosR=(TH3F*)forward->Get("hIntDistortionPosR");
  hPosPhi=(TH3F*)forward->Get("hIntDistortionPosP");

  //negative side (z = -105.5 to 0)
  hNegX=(TH3F*)forward->Get("hIntDistortionNegX");
  hNegY=(TH3F*)forward->Get("hIntDistortionNegY");
  hNegZ=(TH3F*)forward->Get("hIntDistortionNegZ");

  hNegR=(TH3F*)forward->Get("hIntDistortionNegR");
  hNegPhi=(TH3F*)forward->Get("hIntDistortionNegP");

  hNegX=(TH3F*)forward->Get("hIntDistortionNegX");
  hNegY=(TH3F*)forward->Get("hIntDistortionNegY");
  hNegZ=(TH3F*)forward->Get("hIntDistortionNegZ");

  hNegR=(TH3F*)forward->Get("hIntDistortionNegR");
  hNegPhi=(TH3F*)forward->Get("hIntDistortionNegP");

  //average distortion file
  average=TFile::Open(averagefilename,"READ"); 
  //average=TFile::Open("/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.average.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ"); 

  hPosXave=(TH3F*)average->Get("hIntDistortionPosX");
  hPosYave=(TH3F*)average->Get("hIntDistortionPosY");
  hPosZave=(TH3F*)average->Get("hIntDistortionPosZ");
  
  hPosRave=(TH3F*)average->Get("hIntDistortionPosR");
  hPosPhiave=(TH3F*)average->Get("hIntDistortionPosP");

  hNegXave=(TH3F*)average->Get("hIntDistortionNegX");
  hNegYave=(TH3F*)average->Get("hIntDistortionNegY");
  hNegZave=(TH3F*)average->Get("hIntDistortionNegZ");
  
  hNegRave=(TH3F*)average->Get("hIntDistortionNegR");
  hNegPhiave=(TH3F*)average->Get("hIntDistortionNegP");
 
  //subtract average from total distortions to study fluctuations
  hPosX->Add(hPosXave,-1);
  hPosY->Add(hPosYave,-1);
  hPosZ->Add(hPosZave,-1);
  
  hPosR->Add(hPosRave,-1);
  hPosPhi->Add(hPosPhiave,-1);

  hNegX->Add(hNegXave,-1);
  hNegY->Add(hNegYave,-1);
  hNegZ->Add(hNegZave,-1);
  
  hNegR->Add(hNegRave,-1);
  hNegPhi->Add(hNegPhiave,-1);
}

int CMCompareRvsIntFluct(int nMaxEvents = -1) {
  Shifter *shifter;
  int nEvents;

  TCanvas *integcomp=new TCanvas("integcompbyregion","CompareRvIntFluctbyRegion",1500,1500);
  TCanvas *intfluctproj=new TCanvas("intfluctproj","influctproj",1500,1500);
  intfluctproj->Divide(3,3);
  
  const char * sourceinputpattern="/sphenix/user/rcorliss/distortion_maps/2021.04/*h_Charge_*.root"; //updated
    
  //find all files that match the input string (includes wildcards)
  TFileCollection *sourcefilelist=new TFileCollection();
  sourcefilelist->Add(sourceinputpattern);
  TString sourcefilename;
  TString averagefilename = "/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.average.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";
  
  //full charge
  const char * fullchargeinputpattern="/sphenix/user/shulga/Work/IBF/DistortionMap/Files/outputFile_75Hz_G4Hits_sHijing_0-12fm_*.root";
  
  TFileCollection *fullchargefilelist=new TFileCollection();
  fullchargefilelist->Add(fullchargeinputpattern);
  TString fullchargefilename;
  TFile *fullcharge;
  
  TH3F *hFluctCharge;

  //smoothed average
  TFile *smoothedave;
  smoothedave=TFile::Open("/gpfs/mnt/gpfs02/sphenix/user/rcorliss/distortion_maps/2021.04/charge/Smoothed_Average_AA_events.root","READ");
  TH3F *hSmoothedAve;
  hSmoothedAve=(TH3F*)smoothedave->Get("h_Charge_evt_0");
  
  //how many events
  if (nMaxEvents<0){
    nEvents=sourcefilelist->GetNFiles();
  } else if(nMaxEvents<sourcefilelist->GetNFiles()){
    nEvents=nMaxEvents;
  } else {
    nEvents= sourcefilelist->GetNFiles();
  }

  for (int ifile=0;ifile < nEvents;ifile++){
    fullchargefilename=((TFileInfo*)(fullchargefilelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();
    fullcharge=TFile::Open(fullchargefilename,"READ");   

    for(int ihist=0;ihist < 2;ihist++){//should be ihist < 10 to run over all
      //for each file, find all histograms in that file.
      sourcefilename=Form("/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.file%d.h_Charge_%d.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root",ifile,ihist);
    
      //create shifter
      shifter = new Shifter(sourcefilename, averagefilename);

      hFluctCharge=(TH3F*)fullcharge->Get(Form("h_Charge_%d",ihist)); // only 0-9 available
      //nbins: x 360, y 159, z 248
      //xmin: 0.0000000, 0.20000000, -1.0550000
      //xmax: 6.2831900, 0.78000000, 1.0550000
      hFluctCharge->Add(hSmoothedAve,-1); 
    
      TFile *plots;

      plots=TFile::Open(Form("CMModelsPhiRFull_Event%d.root",(10*ifile+ihist)),"READ");

      TH3F *hCylCMModelPhiRPos[2];
      hCylCMModelPhiRPos[0]=(TH3F*)plots->Get("hCMModelR_PhiR_Pos");
      hCylCMModelPhiRPos[1]=(TH3F*)plots->Get("hCMModelPhi_PhiR_Pos");

      //for forward only

      //same range and bins for each coordinate, binned in cm
      //hardcoded numbers from average distortion file's hIntDistortionPosX and hIntDistortionNegX
      int nphi = 82;
      int nr = 54;
      int nz = 82;
    
      double minphi = -0.078539819;
      double minr = 18.884615;
      double minzPos = -1.3187500;
      double minzNeg = -106.81875;
    
      double maxphi = 6.3617253;
      double maxr = 79.115387;
      double maxzPos = 106.81875;
      double maxzNeg = 1.3187500;

      int ndiff = 300;
      int mindiff = -20;
      int maxdiff = 20;

      int nbinsint = 50;

      //compare linear model to integrated charge
      TFile *intFluct;
      TH3F *hIntFluctChargeSmallRPos, *hIntFluctChargeLargeRPos;
   
      intFluct=TFile::Open(Form("IntFluctEvent%d.root", (10*ifile + ihist)), "READ");//change 2 to 10 when looping over all ihist
    
      hIntFluctChargeSmallRPos=(TH3F*)intFluct->Get("hIntFluctChargeSmallRPos");
      hIntFluctChargeLargeRPos=(TH3F*)intFluct->Get("hIntFluctChargeLargeRPos");

      //looking at each region separately
      
      //true vs int fluct inside
      TH2F *hCompareRTruevIntFluctSmallRPosRegion[9];
      hCompareRTruevIntFluctSmallRPosRegion[0] = new TH2F("hCompareRTruevIntFluctSmallRPosRegion0", "True vs Int Fluct Region 0, Pos Side, R Inside; int fluct charge inside (ions); true shift (#mum)",nbinsint,-14e7,-4e7,nbinsint,-10,10);
      hCompareRTruevIntFluctSmallRPosRegion[1] = new TH2F("hCompareRTruevIntFluctSmallRPosRegion1", "True vs Int Fluct Region 1, Pos Side, R Inside; int fluct charge inside (ions); true shift (#mum)",nbinsint,-13e7,-2e7,nbinsint,-35,15);
      hCompareRTruevIntFluctSmallRPosRegion[2] = new TH2F("hCompareRTruevIntFluctSmallRPosRegion2", "True vs Int Fluct Region 2, Pos Side, R Inside; int fluct charge inside (ions); true shift (#mum)",nbinsint,-12e7,1e7,nbinsint,-2,2);
      hCompareRTruevIntFluctSmallRPosRegion[3] = new TH2F("hCompareRTruevIntFluctSmallRPosRegion3", "True vs Int Fluct Region 3, Pos Side, R Inside; int fluct charge inside (ions); true shift (#mum)",nbinsint,-1e8,1e7,nbinsint,-20,20);
      hCompareRTruevIntFluctSmallRPosRegion[4] = new TH2F("hCompareRTruevIntFluctSmallRPosRegion4", "True vs Int Fluct Region 4, Pos Side, R Inside; int fluct charge inside (ions); true shift (#mum)",nbinsint,-1e8,1e7,nbinsint,-25,20);
      hCompareRTruevIntFluctSmallRPosRegion[5] = new TH2F("hCompareRTruevIntFluctSmallRPosRegion5", "True vs Int Fluct Region 5, Pos Side, R Inside; int fluct charge inside (ions); true shift (#mum)",nbinsint,-1e8,1e7,nbinsint,-5,5);
      hCompareRTruevIntFluctSmallRPosRegion[6] = new TH2F("hCompareRTruevIntFluctSmallRPosRegion6", "True vs Int Fluct Region 6, Pos Side, R Inside; int fluct charge inside (ions); true shift (#mum)",nbinsint,-3e7,5e6,nbinsint,-25,30);
      hCompareRTruevIntFluctSmallRPosRegion[7] = new TH2F("hCompareRTruevIntFluctSmallRPosRegion7", "True vs Int Fluct Region 7, Pos Side, R Inside; int fluct charge inside (ions); true shift (#mum)",nbinsint,-3e7,5e6,nbinsint,-30,25);
      hCompareRTruevIntFluctSmallRPosRegion[8] = new TH2F("hCompareRTruevIntFluctSmallRPosRegion8", "True vs Int Fluct Region 8, Pos Side, R Inside; int fluct charge inside (ions); true shift (#mum)",nbinsint,-3e7,5e6,nbinsint,-5,5);

      //true vs int fluct outside
      TH2F *hCompareRTruevIntFluctLargeRPosRegion[9];
      hCompareRTruevIntFluctLargeRPosRegion[0] = new TH2F("hCompareRTruevIntFluctLargeRPosRegion0", "True vs Int Fluct Region 0, Pos Side, R Outside; int fluct charge outside (ions); true shift (#mum)",nbinsint,-13e6,3e6,nbinsint,-10,10);
      hCompareRTruevIntFluctLargeRPosRegion[1] = new TH2F("hCompareRTruevIntFluctLargeRPosRegion1", "True vs Int Fluct Region 1, Pos Side, R Outside; int fluct charge outside (ions); true shift (#mum)",nbinsint,-13e6,3e6,nbinsint,-35,15);
      hCompareRTruevIntFluctLargeRPosRegion[2] = new TH2F("hCompareRTruevIntFluctLargeRPosRegion2", "True vs Int Fluct Region 2, Pos Side, R Outside; int fluct charge outside (ions); true shift (#mum)",nbinsint,-1e7,2e6,nbinsint,-2,2);
      hCompareRTruevIntFluctLargeRPosRegion[3] = new TH2F("hCompareRTruevIntFluctLargeRPosRegion3", "True vs Int Fluct Region 3, Pos Side, R Outside; int fluct charge outside (ions); true shift (#mum)",nbinsint,-1e8,1e7,nbinsint,-25,20);
      hCompareRTruevIntFluctLargeRPosRegion[4] = new TH2F("hCompareRTruevIntFluctLargeRPosRegion4", "True vs Int Fluct Region 4, Pos Side, R Outside; int fluct charge outside (ions); true shift (#mum)",nbinsint,-1e8,1e7,nbinsint,-25,20);
      hCompareRTruevIntFluctLargeRPosRegion[5] = new TH2F("hCompareRTruevIntFluctLargeRPosRegion5", "True vs Int Fluct Region 5, Pos Side, R Outside; int fluct charge outside (ions); true shift (#mum)",nbinsint,-1e8,1e7,nbinsint,-6,6);
      hCompareRTruevIntFluctLargeRPosRegion[6] = new TH2F("hCompareRTruevIntFluctLargeRPosRegion6", "True vs Int Fluct Region 6, Pos Side, R Outside; int fluct charge outside (ions); true shift (#mum)",nbinsint,-13e7,-3e7,nbinsint,-30,30);
      hCompareRTruevIntFluctLargeRPosRegion[7] = new TH2F("hCompareRTruevIntFluctLargeRPosRegion7", "True vs Int Fluct Region 7, Pos Side, R Outside; int fluct charge outside (ions); true shift (#mum)",nbinsint,-13e7,-3e7,nbinsint,-30,30);
      hCompareRTruevIntFluctLargeRPosRegion[8] = new TH2F("hCompareRTruevIntFluctLargeRPosRegion8", "True vs Int Fluct Region 8, Pos Side, R Outside; int fluct charge outside (ions); true shift (#mum)",nbinsint,-11e7,1e7,nbinsint,-5,5);
      
      //true vs net int fluct
      TH2F *hCompareRTruevIntFluctDiffRPosRegion[9];
      hCompareRTruevIntFluctDiffRPosRegion[0] = new TH2F("hCompareRTruevIntFluctDiffRPosRegion0", "True vs Net Int Fluct Region 0, Pos Side; net int fluct charge (ions); true shift (#mum)",nbinsint,-14e7,-2e7,nbinsint,-8,10);
      hCompareRTruevIntFluctDiffRPosRegion[1] = new TH2F("hCompareRTruevIntFluctDiffRPosRegion1", "True vs Net Int Fluct Region 1, Pos Side; net int fluct charge (ions); true shift (#mum)",nbinsint,-14e7,-2e7,nbinsint,-5,15);
      hCompareRTruevIntFluctDiffRPosRegion[2] = new TH2F("hCompareRTruevIntFluctDiffRPosRegion2", "True vs Net Int Fluct Region 2, Pos Side; net int fluct charge (ions); true shift (#mum)",nbinsint,-15e7,1e7,nbinsint,-2,2);
      hCompareRTruevIntFluctDiffRPosRegion[3] = new TH2F("hCompareRTruevIntFluctDiffRPosRegion3", "True vs Net Int Fluct Region 3, Pos Side; net int fluct charge (ions); true shift (#mum)",nbinsint,-15e7,1e8,nbinsint,-25,20);
      hCompareRTruevIntFluctDiffRPosRegion[4] = new TH2F("hCompareRTruevIntFluctDiffRPosRegion4", "True vs Net Int Fluct Region 4, Pos Side; net int fluct charge (ions); true shift (#mum)",nbinsint,-15e7,1e8,nbinsint,-25,20);
      hCompareRTruevIntFluctDiffRPosRegion[5] = new TH2F("hCompareRTruevIntFluctDiffRPosRegion5", "True vs Net Int Fluct Region 5, Pos Side; net int fluct charge (ions); true shift (#mum)",nbinsint,-1e8,6e7,nbinsint,-5,5);
      hCompareRTruevIntFluctDiffRPosRegion[6] = new TH2F("hCompareRTruevIntFluctDiffRPosRegion6", "True vs Net Int Fluct Region 6, Pos Side; net int fluct charge (ions); true shift (#mum)",nbinsint,2e7,15e7,nbinsint,-25,30);
      hCompareRTruevIntFluctDiffRPosRegion[7] = new TH2F("hCompareRTruevIntFluctDiffRPosRegion7", "True vs Net Int Fluct Region 7, Pos Side; net int fluct charge (ions); true shift (#mum)",nbinsint,2e7,15e7,nbinsint,-30,25);
      hCompareRTruevIntFluctDiffRPosRegion[8] = new TH2F("hCompareRTruevIntFluctDiffRPosRegion8", "True vs Net Int Fluct Region 8, Pos Side; net int fluct charge (ions); true shift (#mum)",nbinsint,-1e7,12e7,nbinsint,-4,4);
      
      //model-true vs int fluct inside
       TH2F *hCompareRDiffvIntFluctSmallRPosRegion[9];
      hCompareRDiffvIntFluctSmallRPosRegion[0] = new TH2F("hCompareRDiffvIntFluctSmallRPosRegion0", "Model-True vs Int Fluct Region 0, Pos Side, R Inside; int fluct charge inside (ions); shift difference (#mum)",nbinsint,-13e7,-3e7,nbinsint,-10,5);
      hCompareRDiffvIntFluctSmallRPosRegion[1] = new TH2F("hCompareRDiffvIntFluctSmallRPosRegion1", "Model-True vs Int Fluct Region 1, Pos Side, R Inside; int fluct charge inside (ions); shift difference (#mum)",nbinsint,-13e7,-3e7,nbinsint,-10,10);
      hCompareRDiffvIntFluctSmallRPosRegion[2] = new TH2F("hCompareRDiffvIntFluctSmallRPosRegion2", "Model-True vs Int Fluct Region 2, Pos Side, R Inside; int fluct charge inside (ions); shift difference (#mum)",nbinsint,-12e7,1e7,nbinsint,-2,2);
      hCompareRDiffvIntFluctSmallRPosRegion[3] = new TH2F("hCompareRDiffvIntFluctSmallRPosRegion3", "Model-True vs Int Fluct Region 3, Pos Side, R Inside; int fluct charge inside (ions); shift difference (#mum)",nbinsint,-15e7,1e7,nbinsint,-10,10);
      hCompareRDiffvIntFluctSmallRPosRegion[4] = new TH2F("hCompareRDiffvIntFluctSmallRPosRegion4", "Model-True vs Int Fluct Region 4, Pos Side, R Inside; int fluct charge inside (ions); shift difference (#mum)",nbinsint,-15e7,1e7,nbinsint,-10,10);
      hCompareRDiffvIntFluctSmallRPosRegion[5] = new TH2F("hCompareRDiffvIntFluctSmallRPosRegion5", "Model-True vs Int Fluct Region 5, Pos Side, R Inside; int fluct charge inside (ions); shift difference (#mum)",nbinsint,-12e7,1e7,nbinsint,-5,5);
      hCompareRDiffvIntFluctSmallRPosRegion[6] = new TH2F("hCompareRDiffvIntFluctSmallRPosRegion6", "Model-True vs Int Fluct Region 6, Pos Side, R Inside; int fluct charge inside (ions); shift difference (#mum)",nbinsint,-4e7,5e6,nbinsint,-15,10);
      hCompareRDiffvIntFluctSmallRPosRegion[7] = new TH2F("hCompareRDiffvIntFluctSmallRPosRegion7", "Model-True vs Int Fluct Region 7, Pos Side, R Inside; int fluct charge inside (ions); shift difference (#mum)",nbinsint,-4e7,5e6,nbinsint,-10,20);
      hCompareRDiffvIntFluctSmallRPosRegion[8] = new TH2F("hCompareRDiffvIntFluctSmallRPosRegion8", "Model-True vs Int Fluct Region 8, Pos Side, R Inside; int fluct charge inside (ions); shift difference (#mum)",nbinsint,-3e7,5e6,nbinsint,-5,5);
      
      //model-true vs int fluct outside
      TH2F *hCompareRDiffvIntFluctLargeRPosRegion[9];
      hCompareRDiffvIntFluctLargeRPosRegion[0] = new TH2F("hCompareRDiffvIntFluctLargeRPosRegion0", "Model-True vs Int Fluct Region 0, Pos Side, R Outside; int fluct charge outside (ions); shift difference (#mum)",nbinsint,-12e6,2e6,nbinsint,-10,4);
      hCompareRDiffvIntFluctLargeRPosRegion[1] = new TH2F("hCompareRDiffvIntFluctLargeRPosRegion1", "Model-True vs Int Fluct Region 1, Pos Side, R Outside; int fluct charge outside (ions); shift difference (#mum)",nbinsint,-12e6,2e6,nbinsint,-12,4);
      hCompareRDiffvIntFluctLargeRPosRegion[2] = new TH2F("hCompareRDiffvIntFluctLargeRPosRegion2", "Model-True vs Int Fluct Region 2, Pos Side, R Outside; int fluct charge outside (ions); shift difference (#mum)",nbinsint,-1e7,2e6,nbinsint,-2,2);
      hCompareRDiffvIntFluctLargeRPosRegion[3] = new TH2F("hCompareRDiffvIntFluctLargeRPosRegion3", "Model-True vs Int Fluct Region 3, Pos Side, R Outside; int fluct charge outside (ions); shift difference (#mum)",nbinsint,-1e8,1e7,nbinsint,-10,10);
      hCompareRDiffvIntFluctLargeRPosRegion[4] = new TH2F("hCompareRDiffvIntFluctLargeRPosRegion4", "Model-True vs Int Fluct Region 4, Pos Side, R Outside; int fluct charge outside (ions); shift difference (#mum)",nbinsint,-1e8,1e7,nbinsint,-10,12);
      hCompareRDiffvIntFluctLargeRPosRegion[5] = new TH2F("hCompareRDiffvIntFluctLargeRPosRegion5", "Model-True vs Int Fluct Region 5, Pos Side, R Outside; int fluct charge outside (ions); shift difference (#mum)",nbinsint,-8e7,1e7,nbinsint,-5,5);
      hCompareRDiffvIntFluctLargeRPosRegion[6] = new TH2F("hCompareRDiffvIntFluctLargeRPosRegion6", "Model-True vs Int Fluct Region 6, Pos Side, R Outside; int fluct charge outside (ions); shift difference (#mum)",nbinsint,-13e7,-4e7,nbinsint,-15,10);
      hCompareRDiffvIntFluctLargeRPosRegion[7] = new TH2F("hCompareRDiffvIntFluctLargeRPosRegion7", "Model-True vs Int Fluct Region 7, Pos Side, R Outside; int fluct charge outside (ions); shift difference (#mum)",nbinsint,-13e7,-3e7,nbinsint,-10,25);
      hCompareRDiffvIntFluctLargeRPosRegion[8] = new TH2F("hCompareRDiffvIntFluctLargeRPosRegion8", "Model-True vs Int Fluct Region 8, Pos Side, R Outside; int fluct charge outside (ions); shift difference (#mum)",nbinsint,-12e7,1e7,nbinsint,-6,6);
      
      //model-true vs net int fluct
      TH2F *hCompareRDiffvIntFluctDiffRPosRegion[9];
      hCompareRDiffvIntFluctDiffRPosRegion[0] = new TH2F("hCompareRDiffvIntFluctDiffRPosRegion0", "Model-True vs Net Int Fluct Region 0, Pos Side; net int fluct charge (ions); shift difference (#mum)",nbinsint,-14e7,-3e7,nbinsint,-10,5);
      hCompareRDiffvIntFluctDiffRPosRegion[1] = new TH2F("hCompareRDiffvIntFluctDiffRPosRegion1", "Model-True vs Net Int Fluct Region 1, Pos Side; net int fluct charge (ions); shift difference (#mum)",nbinsint,-14e7,-3e7,nbinsint,-10,10);
      hCompareRDiffvIntFluctDiffRPosRegion[2] = new TH2F("hCompareRDiffvIntFluctDiffRPosRegion2", "Model-True vs Net Int Fluct Region 2, Pos Side; net int fluct charge (ions); shift difference (#mum)",nbinsint,-15e7,1e7,nbinsint,-2,2);
      hCompareRDiffvIntFluctDiffRPosRegion[3] = new TH2F("hCompareRDiffvIntFluctDiffRPosRegion3", "Model-True vs Net Int Fluct Region 3, Pos Side; net int fluct charge (ions); shift difference (#mum)",nbinsint,-15e7,1e8,nbinsint,-8,8);
      hCompareRDiffvIntFluctDiffRPosRegion[4] = new TH2F("hCompareRDiffvIntFluctDiffRPosRegion4", "Model-True vs Net Int Fluct Region 4, Pos Side; net int fluct charge (ions); shift difference (#mum)",nbinsint,-15e7,1e8,nbinsint,-10,12);
      hCompareRDiffvIntFluctDiffRPosRegion[5] = new TH2F("hCompareRDiffvIntFluctDiffRPosRegion5", "Model-True vs Net Int Fluct Region 5, Pos Side; net int fluct charge (ions); shift difference (#mum)",nbinsint,-1e8,8e7,nbinsint,-5,5);
      hCompareRDiffvIntFluctDiffRPosRegion[6] = new TH2F("hCompareRDiffvIntFluctDiffRPosRegion6", "Model-True vs Net Int Fluct Region 6, Pos Side; net int fluct charge (ions); shift difference (#mum)",nbinsint,1e7,15e7,nbinsint,-15,10);
      hCompareRDiffvIntFluctDiffRPosRegion[7] = new TH2F("hCompareRDiffvIntFluctDiffRPosRegion7", "Model-True vs Net Int Fluct Region 7, Pos Side; net int fluct charge (ions); shift difference (#mum)",nbinsint,1e7,15e7,nbinsint,-10,20);
      hCompareRDiffvIntFluctDiffRPosRegion[8] = new TH2F("hCompareRDiffvIntFluctDiffRPosRegion8", "Model-True vs Net Int Fluct Region 8, Pos Side; net int fluct charge (ions); shift difference (#mum)",nbinsint,-1e7,12e7,nbinsint,-6,6);

      //int fluct by region
      // small r
      TH3F *hIntFluctSmallRPosRegion[9];
      hIntFluctSmallRPosRegion[0] = new TH3F("hIntFluctSmallRPosRegion0","Int Fluct Charge (ions), Pos Side, R Inside, Region 0; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctSmallRPosRegion[1] = new TH3F("hIntFluctSmallRPosRegion1","Int Fluct Charge (ions), Pos Side, R Inside, Region 1; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctSmallRPosRegion[2] = new TH3F("hIntFluctSmallRPosRegion2","Int Fluct Charge (ions), Pos Side, R Inside, Region 2; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctSmallRPosRegion[3] = new TH3F("hIntFluctSmallRPosRegion3","Int Fluct Charge (ions), Pos Side, R Inside, Region 3; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctSmallRPosRegion[4] = new TH3F("hIntFluctSmallRPosRegion4","Int Fluct Charge (ions), Pos Side, R Inside, Region 4; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctSmallRPosRegion[5] = new TH3F("hIntFluctSmallRPosRegion5","Int Fluct Charge (ions), Pos Side, R Inside, Region 5; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctSmallRPosRegion[6] = new TH3F("hIntFluctSmallRPosRegion6","Int Fluct Charge (ions), Pos Side, R Inside, Region 6; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctSmallRPosRegion[7] = new TH3F("hIntFluctSmallRPosRegion7","Int Fluct Charge (ions), Pos Side, R Inside, Region 7; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);   
      hIntFluctSmallRPosRegion[8] = new TH3F("hIntFluctSmallRPosRegion8","Int Fluct Charge (ions), Pos Side, R Inside, Region 8; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);

       // large r
      TH3F *hIntFluctLargeRPosRegion[9];
      hIntFluctLargeRPosRegion[0] = new TH3F("hIntFluctLargeRPosRegion0","Int Fluct Charge (ions), Pos Side, R Outside, Region 0; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctLargeRPosRegion[1] = new TH3F("hIntFluctLargeRPosRegion1","Int Fluct Charge (ions), Pos Side, R Outside, Region 1; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctLargeRPosRegion[2] = new TH3F("hIntFluctLargeRPosRegion2","Int Fluct Charge (ions), Pos Side, R Outside, Region 2; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctLargeRPosRegion[3] = new TH3F("hIntFluctLargeRPosRegion3","Int Fluct Charge (ions), Pos Side, R Outside, Region 3; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctLargeRPosRegion[4] = new TH3F("hIntFluctLargeRPosRegion4","Int Fluct Charge (ions), Pos Side, R Outside, Region 4; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctLargeRPosRegion[5] = new TH3F("hIntFluctLargeRPosRegion5","Int Fluct Charge (ions), Pos Side, R Outside, Region 5; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctLargeRPosRegion[6] = new TH3F("hIntFluctLargeRPosRegion6","Int Fluct Charge (ions), Pos Side, R Outside, Region 6; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctLargeRPosRegion[7] = new TH3F("hIntFluctLargeRPosRegion7","Int Fluct Charge (ions), Pos Side, R Outside, Region 7; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);   
      hIntFluctLargeRPosRegion[8] = new TH3F("hIntFluctLargeRPosRegion8","Int Fluct Charge (ions), Pos Side, R Outside, Region 8; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);

        // net
      TH3F *hIntFluctDiffRPosRegion[9];
      hIntFluctDiffRPosRegion[0] = new TH3F("hIntFluctDiffRPosRegion0","Net Int Fluct Charge (ions), Pos Side, Region 0; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctDiffRPosRegion[1] = new TH3F("hIntFluctDiffRPosRegion1","Net Int Fluct Charge (ions), Pos Side, Region 1; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctDiffRPosRegion[2] = new TH3F("hIntFluctDiffRPosRegion2","Net Int Fluct Charge (ions), Pos Side, Region 2; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctDiffRPosRegion[3] = new TH3F("hIntFluctDiffRPosRegion3","Net Int Fluct Charge (ions), Pos Side, Region 3; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctDiffRPosRegion[4] = new TH3F("hIntFluctDiffRPosRegion4","Net Int Fluct Charge (ions), Pos Side, Region 4; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctDiffRPosRegion[5] = new TH3F("hIntFluctDiffRPosRegion5","Net Int Fluct Charge (ions), Pos Side, Region 5; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctDiffRPosRegion[6] = new TH3F("hIntFluctDiffRPosRegion6","Net Int Fluct Charge (ions), Pos Side, Region 6; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
      hIntFluctDiffRPosRegion[7] = new TH3F("hIntFluctDiffRPosRegion7","Net Int Fluct Charge (ions), Pos Side, Region 7; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);   
      hIntFluctDiffRPosRegion[8] = new TH3F("hIntFluctDiffRPosRegion8","Net Int Fluct Charge (ions), Pos Side, Region 8; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos); 

      // true vs int fluct localized plots
      // inner R, region 4, r,z plot
      TH2F *hLocalRTruevIntFluctSmallRPosRegion4_RZ = new TH2F("hLocalRTruevIntFluctSmallRPosRegion4_RZ","R,Z Close-up of True vs Int Fluct, Pos Side, R Inside, Region 4; r (cm); z (cm)",nr,minr,maxr, nz,minzPos,maxzPos);
      // inner R, region 4, phi,z plot
      TH2F *hLocalRTruevIntFluctSmallRPosRegion4_PhiZ = new TH2F("hLocalRTruevIntFluctSmallRPosRegion4_PhiZ","Phi,Z Close-up of True vs Int Fluct, Pos Side, R Inside, Region 4; r (cm); z (cm)",nphi,minphi,maxphi, nz,minzPos,maxzPos);
      // inner R, region 6, lower peak, r,z
      TH2F *hLocalRTruevIntFluctSmallRPosRegion6_Lower = new TH2F("hLocalRTruevIntFluctSmallRPosRegion6_Lower","R,Z Lower Close-up of True vs Int Fluct, Pos Side, R Inside, Region 6; r (cm); z (cm)",nr,minr,maxr, nz,minzPos,maxzPos);
      // inner R, region 6, upper peak, r,z
      TH2F *hLocalRTruevIntFluctSmallRPosRegion6_Upper = new TH2F("hLocalRTruevIntFluctSmallRPosRegion6_Upper","R,Z Upper Close-up of True vs Int Fluct, Pos Side, R Inside, Region 6; r (cm); z (cm)",nr,minr,maxr, nz,minzPos,maxzPos);
      // inner R, region 8, r,z
      TH2F *hLocalRTruevIntFluctSmallRPosRegion8 = new TH2F("hLocalRTruevIntFluctSmallRPosRegion8","R,Z Close-up of True vs Int Fluct, Pos Side, R Inside, Region 8; r (cm); z (cm)",nr,minr,maxr, nz,minzPos,maxzPos);
      // outer R, region 4, r,z
      TH2F *hLocalRTruevIntFluctLargeRPosRegion4 = new TH2F("hLocalRTruevIntFluctLargeRPosRegion4","R,Z Close-up of True vs Int Fluct, Pos Side, R Outside, Region 4; r (cm); z (cm)",nr,minr,maxr, nz,minzPos,maxzPos);
      // net, region 4, r,z
      TH2F *hLocalRTruevIntFluctDiffRPosRegion4 = new TH2F("hLocalRTruevIntFluctDiffRPosRegion4","R,Z Close-up of True vs Net Int Fluct, Pos Side, Region 4; r (cm); z (cm)",nr,minr,maxr, nz,minzPos,maxzPos);
      
      for(int i = 1; i < nphi - 1; i++){
	double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
	for(int j = 1; j < nr - 1; j++){
	  double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin
	  for(int k = 1; k < nz - 1; k++){
	    double zPos = minzPos + ((maxzPos - minzPos)/(1.0*nz))*(k+0.5); //center of bin
	    double zNeg = minzNeg + ((maxzNeg - minzNeg)/(1.0*nz))*(k+0.5); //center of bin

	    //positive
	   
	    double shifttrueCylPos[2];
	    double shiftrecoCylPhiRPos[2];
	    double differenceCylPhiRPos[2];

	    int binPhiRPos = hCylCMModelPhiRPos[0]->FindBin(phi,r,zPos);

	    double intfluctchargeSmallRPos, intfluctchargeLargeRPos, intfluctchargeDiffRPos;

	    if(r >= 70.0){
	      if(zPos <= 15.0){ //0
		//int fluct
		intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
		intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
		intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;

		hIntFluctSmallRPosRegion[0]->Fill(phi,r,zPos,intfluctchargeSmallRPos);
		hIntFluctLargeRPosRegion[0]->Fill(phi,r,zPos,intfluctchargeLargeRPos);
		hIntFluctDiffRPosRegion[0]->Fill(phi,r,zPos,intfluctchargeDiffRPos);
		  
		//true
		shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
		shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);

		//model-true
		shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
		differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 

		shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
		differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    
		//fill region 0 plots
		//true
		hCompareRTruevIntFluctSmallRPosRegion[0]->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctLargeRPosRegion[0]->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctDiffRPosRegion[0]->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);
		//model-true
		hCompareRDiffvIntFluctSmallRPosRegion[0]->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctLargeRPosRegion[0]->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctDiffRPosRegion[0]->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);
			
	      }else if((zPos > 15.0) && (zPos < 90.0)){ //1
		//int fluct
		intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
		intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
		intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;

		hIntFluctSmallRPosRegion[1]->Fill(phi,r,zPos,intfluctchargeSmallRPos);
		hIntFluctLargeRPosRegion[1]->Fill(phi,r,zPos,intfluctchargeLargeRPos);
		hIntFluctDiffRPosRegion[1]->Fill(phi,r,zPos,intfluctchargeDiffRPos);
		
		//true
		shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
		shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);

		//model-true
		shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
		differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 

		shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
		differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    
		//fill region 1 plots
		//true
		hCompareRTruevIntFluctSmallRPosRegion[1]->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctLargeRPosRegion[1]->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctDiffRPosRegion[1]->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);
		//model-true
		hCompareRDiffvIntFluctSmallRPosRegion[1]->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctLargeRPosRegion[1]->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctDiffRPosRegion[1]->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);
		
	      }else if(zPos >= 90.0){ //2
		//int fluct
		intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
		intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
		intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;
		
		hIntFluctSmallRPosRegion[2]->Fill(phi,r,zPos,intfluctchargeSmallRPos);
		hIntFluctLargeRPosRegion[2]->Fill(phi,r,zPos,intfluctchargeLargeRPos);
		hIntFluctDiffRPosRegion[2]->Fill(phi,r,zPos,intfluctchargeDiffRPos);

		//true
		shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
		shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);

		//model-true
		shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
		differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 

		shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
		differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    
		//fill region 2 plots
		//true
		hCompareRTruevIntFluctSmallRPosRegion[2]->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctLargeRPosRegion[2]->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctDiffRPosRegion[2]->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);
		//model-true
		hCompareRDiffvIntFluctSmallRPosRegion[2]->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctLargeRPosRegion[2]->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctDiffRPosRegion[2]->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);
		
	      }
	    } else if((r > 30.0) && (r < 70.0)) {
	      if(zPos <= 15.0){ //3
		//int fluct
		intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
		intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
		intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;

		hIntFluctSmallRPosRegion[3]->Fill(phi,r,zPos,intfluctchargeSmallRPos);
		hIntFluctLargeRPosRegion[3]->Fill(phi,r,zPos,intfluctchargeLargeRPos);
		hIntFluctDiffRPosRegion[3]->Fill(phi,r,zPos,intfluctchargeDiffRPos);
		
		//true
		shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
		shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);

		//model-true
		shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
		differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 

		shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
		differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    
		//fill region 3 plots
		//true
		hCompareRTruevIntFluctSmallRPosRegion[3]->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctLargeRPosRegion[3]->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctDiffRPosRegion[3]->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);
		//model-true
		hCompareRDiffvIntFluctSmallRPosRegion[3]->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctLargeRPosRegion[3]->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctDiffRPosRegion[3]->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);
		
	      }else if((zPos > 15.0) && (zPos < 90.0)){//4
		//int fluct
		intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
		intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
		intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;

		hIntFluctSmallRPosRegion[4]->Fill(phi,r,zPos,intfluctchargeSmallRPos);
		hIntFluctLargeRPosRegion[4]->Fill(phi,r,zPos,intfluctchargeLargeRPos);
		hIntFluctDiffRPosRegion[4]->Fill(phi,r,zPos,intfluctchargeDiffRPos);
		
		//true
		shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
		shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);

		//model-true
		shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
		differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 

		shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
		differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    
		//fill region 4 plots
		//true
		hCompareRTruevIntFluctSmallRPosRegion[4]->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctLargeRPosRegion[4]->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctDiffRPosRegion[4]->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);
		//model-true
		hCompareRDiffvIntFluctSmallRPosRegion[4]->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctLargeRPosRegion[4]->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctDiffRPosRegion[4]->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);

		//localized plots
		if((-1.0 < shifttrueCylPos[0]) && (shifttrueCylPos[0] < 1.0) && (-9e7 < intfluctchargeSmallRPos) && (intfluctchargeSmallRPos < -6e7)){
		  hLocalRTruevIntFluctSmallRPosRegion4_RZ->Fill(r,zPos,1);
		  hLocalRTruevIntFluctSmallRPosRegion4_PhiZ->Fill(phi,zPos,1);
		}else if((-1.0 < shifttrueCylPos[0]) && (shifttrueCylPos[0] < 1.0) && (-3e7 < intfluctchargeLargeRPos) && (intfluctchargeLargeRPos < 0.0)){
		  hLocalRTruevIntFluctLargeRPosRegion4->Fill(r,zPos,1);
		  //	}else if((-1.0 < shifttrueCylPos[0]) && (shifttrueCylPos[0] < 1.0) && (-9e7 < intfluctchargeDiffRPos) && (intfluctchargeDiffRPos < -3e7)){
		}else if((-1.0 < shifttrueCylPos[0]) && (shifttrueCylPos[0] < 1.0) && (-9e7 < intfluctchargeDiffRPos) && (intfluctchargeDiffRPos < -2e7)){
		  hLocalRTruevIntFluctDiffRPosRegion4->Fill(r,zPos,1);
		}
		
	      }else if(zPos >= 90.0){//5
		//int fluct
		intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
		intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
		intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;

		hIntFluctSmallRPosRegion[5]->Fill(phi,r,zPos,intfluctchargeSmallRPos);
		hIntFluctLargeRPosRegion[5]->Fill(phi,r,zPos,intfluctchargeLargeRPos);
		hIntFluctDiffRPosRegion[5]->Fill(phi,r,zPos,intfluctchargeDiffRPos);
		
		//true
		shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
		shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);

		//model-true
		shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
		differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 

		shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
		differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    
		//fill region 5 plots
		//true
		hCompareRTruevIntFluctSmallRPosRegion[5]->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctLargeRPosRegion[5]->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctDiffRPosRegion[5]->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);
		//model-true
		hCompareRDiffvIntFluctSmallRPosRegion[5]->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctLargeRPosRegion[5]->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctDiffRPosRegion[5]->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);
		
	      }
	    }else if(r <= 30.0) {
	      if(zPos <= 15.0){ //6
		//int fluct
		intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
		intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
		intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;

		hIntFluctSmallRPosRegion[6]->Fill(phi,r,zPos,intfluctchargeSmallRPos);
		hIntFluctLargeRPosRegion[6]->Fill(phi,r,zPos,intfluctchargeLargeRPos);
		hIntFluctDiffRPosRegion[6]->Fill(phi,r,zPos,intfluctchargeDiffRPos);
		
		//true
		shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
		shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);

		//model-true
		shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
		differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 

		shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
		differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    
		//fill region 6 plots
		//true
		hCompareRTruevIntFluctSmallRPosRegion[6]->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctLargeRPosRegion[6]->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctDiffRPosRegion[6]->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);
		//model-true
		hCompareRDiffvIntFluctSmallRPosRegion[6]->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctLargeRPosRegion[6]->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctDiffRPosRegion[6]->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);

		//localized plots
		if((6.0 < shifttrueCylPos[0]) && (shifttrueCylPos[0] < 8.0) && (-1e6 < intfluctchargeSmallRPos) && (intfluctchargeSmallRPos < -1e6)){
		  hLocalRTruevIntFluctSmallRPosRegion6_Lower->Fill(r,zPos,1);
		}else if((-8.0 < shifttrueCylPos[0]) && (shifttrueCylPos[0] < -5.0) && (-1e6 < intfluctchargeSmallRPos) && (intfluctchargeSmallRPos < -1e6)){
		  hLocalRTruevIntFluctSmallRPosRegion6_Upper->Fill(r,zPos,1);
		}
		
	      }else if((zPos > 15.0) && (zPos < 90.0)){ //7
		//int fluct
		intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
		intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
		intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;

		hIntFluctSmallRPosRegion[7]->Fill(phi,r,zPos,intfluctchargeSmallRPos);
		hIntFluctLargeRPosRegion[7]->Fill(phi,r,zPos,intfluctchargeLargeRPos);
		hIntFluctDiffRPosRegion[7]->Fill(phi,r,zPos,intfluctchargeDiffRPos);
		
		//true
		shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
		shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);

		//model-true
		shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
		differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 

		shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
		differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    
		//fill region 7 plots
		//true
		hCompareRTruevIntFluctSmallRPosRegion[7]->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctLargeRPosRegion[7]->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctDiffRPosRegion[7]->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);
		//model-true
		hCompareRDiffvIntFluctSmallRPosRegion[7]->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctLargeRPosRegion[7]->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctDiffRPosRegion[7]->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);
		
	      }else if(zPos >= 90.0){ //8
		//int fluct
		intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
		intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
		intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;

		hIntFluctSmallRPosRegion[8]->Fill(phi,r,zPos,intfluctchargeSmallRPos);
		hIntFluctLargeRPosRegion[8]->Fill(phi,r,zPos,intfluctchargeLargeRPos);
		hIntFluctDiffRPosRegion[8]->Fill(phi,r,zPos,intfluctchargeDiffRPos);
		
		//true
		shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
		shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);

		//model-true
		shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
		differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 

		shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
		differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    
		//fill region 8 plots
		//true
		hCompareRTruevIntFluctSmallRPosRegion[8]->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctLargeRPosRegion[8]->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
		hCompareRTruevIntFluctDiffRPosRegion[8]->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);
		//model-true
		hCompareRDiffvIntFluctSmallRPosRegion[8]->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctLargeRPosRegion[8]->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
		hCompareRDiffvIntFluctDiffRPosRegion[8]->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);

			//localized plot
		if((1.0 < shifttrueCylPos[0]) && (shifttrueCylPos[0] < 2.0) && (-1e6 < intfluctchargeSmallRPos) && (intfluctchargeSmallRPos < -1e6)){
		  hLocalRTruevIntFluctSmallRPosRegion8->Fill(r,zPos,1);
		}
		
	      }
	    }
	  }
	}
      }
    
      gStyle->SetOptStat("oue");

      //set up canvas
      //int fluct inner R xz proj for each region
      hIntFluctSmallRPosRegion[0]->Project3D("yz")->SetStats(0);
      hIntFluctSmallRPosRegion[1]->Project3D("yz")->SetStats(0);
      hIntFluctSmallRPosRegion[2]->Project3D("yz")->SetStats(0);
      hIntFluctSmallRPosRegion[3]->Project3D("yz")->SetStats(0);
      hIntFluctSmallRPosRegion[4]->Project3D("yz")->SetStats(0);
      hIntFluctSmallRPosRegion[5]->Project3D("yz")->SetStats(0);
      hIntFluctSmallRPosRegion[6]->Project3D("yz")->SetStats(0);
      hIntFluctSmallRPosRegion[7]->Project3D("yz")->SetStats(0);
      hIntFluctSmallRPosRegion[8]->Project3D("yz")->SetStats(0);

      //int fluct outer R xz proj for each region
      hIntFluctLargeRPosRegion[0]->Project3D("yz")->SetStats(0);
      hIntFluctLargeRPosRegion[1]->Project3D("yz")->SetStats(0);
      hIntFluctLargeRPosRegion[2]->Project3D("yz")->SetStats(0);
      hIntFluctLargeRPosRegion[3]->Project3D("yz")->SetStats(0);
      hIntFluctLargeRPosRegion[4]->Project3D("yz")->SetStats(0);
      hIntFluctLargeRPosRegion[5]->Project3D("yz")->SetStats(0);
      hIntFluctLargeRPosRegion[6]->Project3D("yz")->SetStats(0);
      hIntFluctLargeRPosRegion[7]->Project3D("yz")->SetStats(0);
      hIntFluctLargeRPosRegion[8]->Project3D("yz")->SetStats(0);

      //net int fluct xz proj for each region
      hIntFluctDiffRPosRegion[0]->Project3D("yz")->SetStats(0);
      hIntFluctDiffRPosRegion[1]->Project3D("yz")->SetStats(0);
      hIntFluctDiffRPosRegion[2]->Project3D("yz")->SetStats(0);
      hIntFluctDiffRPosRegion[3]->Project3D("yz")->SetStats(0);
      hIntFluctDiffRPosRegion[4]->Project3D("yz")->SetStats(0);
      hIntFluctDiffRPosRegion[5]->Project3D("yz")->SetStats(0);
      hIntFluctDiffRPosRegion[6]->Project3D("yz")->SetStats(0);
      hIntFluctDiffRPosRegion[7]->Project3D("yz")->SetStats(0);
      hIntFluctDiffRPosRegion[8]->Project3D("yz")->SetStats(0);
     
      // intfluctproj->cd();
      //inner r
      intfluctproj->cd(1);
      hIntFluctSmallRPosRegion[0]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(2);
      hIntFluctSmallRPosRegion[1]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(3);
      hIntFluctSmallRPosRegion[2]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(4);
      hIntFluctSmallRPosRegion[3]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(5);
      hIntFluctSmallRPosRegion[4]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(6);
      hIntFluctSmallRPosRegion[5]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(7);
      hIntFluctSmallRPosRegion[6]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(8);
      hIntFluctSmallRPosRegion[7]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(9);
      hIntFluctSmallRPosRegion[8]->Project3D("xz")->Draw("colz");

      intfluctproj->Print(Form("intfluctprojEvent%d.pdf(",(10*ifile + ihist)),"pdf");

      //outer r
      intfluctproj->cd(1);
      hIntFluctLargeRPosRegion[0]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(2);
      hIntFluctLargeRPosRegion[1]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(3);
      hIntFluctLargeRPosRegion[2]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(4);
      hIntFluctLargeRPosRegion[3]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(5);
      hIntFluctLargeRPosRegion[4]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(6);
      hIntFluctLargeRPosRegion[5]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(7);
      hIntFluctLargeRPosRegion[6]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(8);
      hIntFluctLargeRPosRegion[7]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(9);
      hIntFluctLargeRPosRegion[8]->Project3D("xz")->Draw("colz");

      intfluctproj->Print(Form("intfluctprojEvent%d.pdf",(10*ifile + ihist)),"pdf");

      //net
      intfluctproj->cd(1);
      hIntFluctDiffRPosRegion[0]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(2);
      hIntFluctDiffRPosRegion[1]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(3);
      hIntFluctDiffRPosRegion[2]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(4);
      hIntFluctDiffRPosRegion[3]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(5);
      hIntFluctDiffRPosRegion[4]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(6);
      hIntFluctDiffRPosRegion[5]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(7);
      hIntFluctDiffRPosRegion[6]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(8);
      hIntFluctDiffRPosRegion[7]->Project3D("xz")->Draw("colz");
      intfluctproj->cd(9);
      hIntFluctDiffRPosRegion[8]->Project3D("xz")->Draw("colz");

      intfluctproj->Print(Form("intfluctprojEvent%d.pdf)",(10*ifile + ihist)),"pdf");
      
      TPad *integcomptitlepad = new TPad("integcomptitlepad","",0.0,0.96,1.0,1.0);
      TPad *integcompplots = new TPad("integcompplotspad","",0.0,0.0,1.0,0.96);

      TLatex *integcomptitle = new TLatex(0.0,0.0,"");

      integcomptitle->SetNDC();
      integcomptitle->SetTextSize(0.4);

      integcomp->cd();
      integcompplots->Draw();
      integcomptitlepad->Draw();

      integcompplots->Divide(3,3);

      //true, inner
      integcompplots->cd(1);
      hCompareRTruevIntFluctSmallRPosRegion[0]->Draw("colz");
      integcompplots->cd(2);
      hCompareRTruevIntFluctSmallRPosRegion[1]->Draw("colz");
      integcompplots->cd(3);
      hCompareRTruevIntFluctSmallRPosRegion[2]->Draw("colz");
      integcompplots->cd(4);
      hCompareRTruevIntFluctSmallRPosRegion[3]->Draw("colz");
      integcompplots->cd(5);
      hCompareRTruevIntFluctSmallRPosRegion[4]->Draw("colz");
      integcompplots->cd(6);
      hCompareRTruevIntFluctSmallRPosRegion[5]->Draw("colz");
      integcompplots->cd(7);
      hCompareRTruevIntFluctSmallRPosRegion[6]->Draw("colz");
      integcompplots->cd(8);
      hCompareRTruevIntFluctSmallRPosRegion[7]->Draw("colz");
      integcompplots->cd(9);
      hCompareRTruevIntFluctSmallRPosRegion[8]->Draw("colz");

      integcomptitlepad->cd();
      integcomptitlepad->Clear();
      integcomptitle->DrawLatex(0.01,0.4,Form("Event %d: True vs Int Fluct Charge, Inner R", (10*ifile + ihist))); 
      integcomptitle->Draw();
      integcomp->Print(Form("integcompEvent%d.pdf(",(10*ifile + ihist)),"pdf");
      
      //true, outer
      integcompplots->cd(1);
      hCompareRTruevIntFluctLargeRPosRegion[0]->Draw("colz");
      integcompplots->cd(2);
      hCompareRTruevIntFluctLargeRPosRegion[1]->Draw("colz");
      integcompplots->cd(3);
      hCompareRTruevIntFluctLargeRPosRegion[2]->Draw("colz");
      integcompplots->cd(4);
      hCompareRTruevIntFluctLargeRPosRegion[3]->Draw("colz");
      integcompplots->cd(5);
      hCompareRTruevIntFluctLargeRPosRegion[4]->Draw("colz");
      integcompplots->cd(6);
      hCompareRTruevIntFluctLargeRPosRegion[5]->Draw("colz");
      integcompplots->cd(7);
      hCompareRTruevIntFluctLargeRPosRegion[6]->Draw("colz");
      integcompplots->cd(8);
      hCompareRTruevIntFluctLargeRPosRegion[7]->Draw("colz");
      integcompplots->cd(9);
      hCompareRTruevIntFluctLargeRPosRegion[8]->Draw("colz");

      integcomptitlepad->cd();
      integcomptitlepad->Clear();
      integcomptitle->DrawLatex(0.01,0.4,Form("Event %d: True vs Int Fluct Charge, Outer R", (10*ifile + ihist))); 
      integcomptitle->Draw();
      integcomp->Print(Form("integcompEvent%d.pdf",(10*ifile + ihist)),"pdf");

      //true, net
      integcompplots->cd(1);
      hCompareRTruevIntFluctDiffRPosRegion[0]->Draw("colz");
      integcompplots->cd(2);
      hCompareRTruevIntFluctDiffRPosRegion[1]->Draw("colz");
      integcompplots->cd(3);
      hCompareRTruevIntFluctDiffRPosRegion[2]->Draw("colz");
      integcompplots->cd(4);
      hCompareRTruevIntFluctDiffRPosRegion[3]->Draw("colz");
      integcompplots->cd(5);
      hCompareRTruevIntFluctDiffRPosRegion[4]->Draw("colz");
      integcompplots->cd(6);
      hCompareRTruevIntFluctDiffRPosRegion[5]->Draw("colz");
      integcompplots->cd(7);
      hCompareRTruevIntFluctDiffRPosRegion[6]->Draw("colz");
      integcompplots->cd(8);
      hCompareRTruevIntFluctDiffRPosRegion[7]->Draw("colz");
      integcompplots->cd(9);
      hCompareRTruevIntFluctDiffRPosRegion[8]->Draw("colz");

      integcomptitlepad->cd();
      integcomptitlepad->Clear();
      integcomptitle->DrawLatex(0.01,0.4,Form("Event %d: True vs Int Fluct Charge, Net", (10*ifile + ihist))); 
      integcomptitle->Draw();
      integcomp->Print(Form("integcompEvent%d.pdf",(10*ifile + ihist)),"pdf");

      //model-true, inner
      integcompplots->cd(1);
      hCompareRDiffvIntFluctSmallRPosRegion[0]->Draw("colz");
      integcompplots->cd(2);
      hCompareRDiffvIntFluctSmallRPosRegion[1]->Draw("colz");
      integcompplots->cd(3);
      hCompareRDiffvIntFluctSmallRPosRegion[2]->Draw("colz");
      integcompplots->cd(4);
      hCompareRDiffvIntFluctSmallRPosRegion[3]->Draw("colz");
      integcompplots->cd(5);
      hCompareRDiffvIntFluctSmallRPosRegion[4]->Draw("colz");
      integcompplots->cd(6);
      hCompareRDiffvIntFluctSmallRPosRegion[5]->Draw("colz");
      integcompplots->cd(7);
      hCompareRDiffvIntFluctSmallRPosRegion[6]->Draw("colz");
      integcompplots->cd(8);
      hCompareRDiffvIntFluctSmallRPosRegion[7]->Draw("colz");
      integcompplots->cd(9);
      hCompareRDiffvIntFluctSmallRPosRegion[8]->Draw("colz");

      integcomptitlepad->cd();
      integcomptitlepad->Clear();
      integcomptitle->DrawLatex(0.01,0.4,Form("Event %d: Model-True vs Int Fluct Charge, Inner R", (10*ifile + ihist))); 
      integcomptitle->Draw();
      integcomp->Print(Form("integcompEvent%d.pdf",(10*ifile + ihist)),"pdf");

      //model-true, outer
      integcompplots->cd(1);
      hCompareRDiffvIntFluctLargeRPosRegion[0]->Draw("colz");
      integcompplots->cd(2);
      hCompareRDiffvIntFluctLargeRPosRegion[1]->Draw("colz");
      integcompplots->cd(3);
      hCompareRDiffvIntFluctLargeRPosRegion[2]->Draw("colz");
      integcompplots->cd(4);
      hCompareRDiffvIntFluctLargeRPosRegion[3]->Draw("colz");
      integcompplots->cd(5);
      hCompareRDiffvIntFluctLargeRPosRegion[4]->Draw("colz");
      integcompplots->cd(6);
      hCompareRDiffvIntFluctLargeRPosRegion[5]->Draw("colz");
      integcompplots->cd(7);
      hCompareRDiffvIntFluctLargeRPosRegion[6]->Draw("colz");
      integcompplots->cd(8);
      hCompareRDiffvIntFluctLargeRPosRegion[7]->Draw("colz");
      integcompplots->cd(9);
      hCompareRDiffvIntFluctLargeRPosRegion[8]->Draw("colz");

      integcomptitlepad->cd();
      integcomptitlepad->Clear();
      integcomptitle->DrawLatex(0.01,0.4,Form("Event %d: Model-True vs Int Fluct Charge, Outer R", (10*ifile + ihist))); 
      integcomptitle->Draw();
      integcomp->Print(Form("integcompEvent%d.pdf",(10*ifile + ihist)),"pdf");

      //model-true, net
      integcompplots->cd(1);
      hCompareRDiffvIntFluctDiffRPosRegion[0]->Draw("colz");
      integcompplots->cd(2);
      hCompareRDiffvIntFluctDiffRPosRegion[1]->Draw("colz");
      integcompplots->cd(3);
      hCompareRDiffvIntFluctDiffRPosRegion[2]->Draw("colz");
      integcompplots->cd(4);
      hCompareRDiffvIntFluctDiffRPosRegion[3]->Draw("colz");
      integcompplots->cd(5);
      hCompareRDiffvIntFluctDiffRPosRegion[4]->Draw("colz");
      integcompplots->cd(6);
      hCompareRDiffvIntFluctDiffRPosRegion[5]->Draw("colz");
      integcompplots->cd(7);
      hCompareRDiffvIntFluctDiffRPosRegion[6]->Draw("colz");
      integcompplots->cd(8);
      hCompareRDiffvIntFluctDiffRPosRegion[7]->Draw("colz");
      integcompplots->cd(9);
      hCompareRDiffvIntFluctDiffRPosRegion[8]->Draw("colz");

      integcomptitlepad->cd();
      integcomptitlepad->Clear();
      integcomptitle->DrawLatex(0.01,0.4,Form("Event %d: Model-True vs Int Fluct Charge, Net", (10*ifile + ihist))); 
      integcomptitle->Draw();
      integcomp->Print(Form("integcompEvent%d.pdf",(10*ifile + ihist)),"pdf");

      // localizations
      integcompplots->cd(1);
      hLocalRTruevIntFluctSmallRPosRegion4_RZ->Draw("colz");
      integcompplots->cd(2);
      hLocalRTruevIntFluctSmallRPosRegion4_PhiZ->Draw("colz");
      integcompplots->cd(3);
      hLocalRTruevIntFluctSmallRPosRegion6_Lower->Draw("colz");
      integcompplots->cd(4);
      hLocalRTruevIntFluctSmallRPosRegion6_Upper->Draw("colz");
      integcompplots->cd(5);
      hLocalRTruevIntFluctSmallRPosRegion8->Draw("colz");
      integcompplots->cd(6);
      hLocalRTruevIntFluctLargeRPosRegion4->Draw("colz");
      integcompplots->cd(7);
      hLocalRTruevIntFluctDiffRPosRegion4->Draw("colz");
      integcompplots->cd(8)->Clear();
      integcompplots->cd(9)->Clear();

      integcomptitlepad->cd();
      integcomptitlepad->Clear();
      integcomptitle->DrawLatex(0.01,0.4,Form("Event %d: True vs Int Fluct Charge Localized on Peaks", (10*ifile + ihist))); 
      integcomptitle->Draw();
      //integcomp->Print(Form("integcompEvent%d.pdf",(10*ifile + ihist)),"pdf");

      //
      
      integcomp->Print(Form("integcompEvent%d.pdf)",(10*ifile + ihist)),"pdf");
    }
  }
  
  return 0;
}
