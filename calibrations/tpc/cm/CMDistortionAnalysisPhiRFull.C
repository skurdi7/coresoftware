//step 3 with phi,r coords
#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"

using namespace std;

//void CMModels{int nphi, double minphi, double maxphi, int nr,double minr,double maxr, int nz, double minzPos,double maxzPos, double minzNeg,double maxzNeg, TH2F *hCartesianAveShiftPhiRPos[3], TH2F *hCylindricalAveShiftPhiRPos[2], TH2F *hCartesianAveShiftPhiRNeg[3], TH2F *hCylindricalAveShiftPhiRNeg[2]};

void WriteIntFluctFile(int ifile, int ihist, int nphi, double minphi, double maxphi, int nr,double minr,double maxr, int nz, double minzPos,double maxzPos, double minzNeg,double maxzNeg, TH3F *hFluctCharge);

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

int CMDistortionAnalysisPhiRFull(int nMaxEvents = -1) {
  Shifter *shifter;
  int nbins = 35; 
  double low = -80.0;
  double high = 80.0;
  double deltaXPos, deltaYPos, deltaZPos, deltaRPos, deltaPhiPos;
  double deltaXNeg, deltaYNeg, deltaZNeg, deltaRNeg, deltaPhiNeg;
  int nEvents; 
  
  TCanvas *canvas=new TCanvas("canvas","CMDistortionAnalysisPhiRNeg",2000,3000);
  TCanvas *integ=new TCanvas("integ","IntegratedFluctAnalysis",1000,1000);
  TCanvas *integcomp=new TCanvas("integcomp","CompareRvIntFluct",1500,1000);
  
  int nsumbins = 20;
  int minsum = -10;
  int maxsum = 10;
  
  //set up summary plots
  //positive 
  TH1F *hDifferenceMeanRPos = new TH1F("hDifferenceMeanR_Pos", "Average Difference between R Model and True of All Events, Positive Side (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
  TH1F *hDifferenceStdDevRPos = new TH1F("hDifferenceStdDevR_Pos", "Std Dev of Difference between R Model and True of All Events, Positive Side (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
  TH1F *hTrueMeanRPos = new TH1F("hTrueMeanR_Pos", "Mean True R Distortion Model of All Events, Positive Side (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
  TH1F *hTrueStdDevRPos = new TH1F("hTrueStdDevR_Pos", "Std Dev of True R Distortion Model of All Events, Positive Side (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
  TH1F *hDifferenceMeanPhiPos = new TH1F("hDifferenceMeanPhi_Pos", "Average Difference between Phi Model and True of All Events, Positive Side (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
  TH1F *hDifferenceStdDevPhiPos = new TH1F("hDifferenceStdDevPhi_Pos", "Std Dev of Difference between Phi Model and True of All Events, Positive Side (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    
  TH1F *hTrueMeanPhiPos = new TH1F("hTrueMeanPhi_Pos", "Mean True Phi Distortion Model of All Events, Positive Side (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
  TH1F *hTrueStdDevPhiPos = new TH1F("hTrueStdDevPhi_Pos", "Std Dev of True Phi Distortion Model of All Events, Positive Side (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);

  //negative
  TH1F *hDifferenceMeanRNeg = new TH1F("hDifferenceMeanR_Neg", "Average Difference between R Model and True of All Events, Negative Side (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
  TH1F *hDifferenceStdDevRNeg = new TH1F("hDifferenceStdDevR_Neg", "Std Dev of Difference between R Model and True of All Events, Negative Side (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
  TH1F *hTrueMeanRNeg = new TH1F("hTrueMeanR_Neg", "Mean True R Distortion Model of All Events, Negative Side (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
  TH1F *hTrueStdDevRNeg = new TH1F("hTrueStdDevR_Neg", "Std Dev of True R Distortion Model of All Events, Negative Side (R > 30); #Delta R (#mum)", nsumbins, minsum, maxsum);
    
  TH1F *hDifferenceMeanPhiNeg = new TH1F("hDifferenceMeanPhi_Neg", "Average Difference between Phi Model and True of All Events, Negative Side (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
  TH1F *hDifferenceStdDevPhiNeg = new TH1F("hDifferenceStdDevPhi_Neg", "Std Dev of Difference between Phi Model and True of All Events, Negative Side (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
    
  TH1F *hTrueMeanPhiNeg = new TH1F("hTrueMeanPhi_Neg", "Mean True Phi Distortion Model of All Events, Negative Side (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);
  TH1F *hTrueStdDevPhiNeg = new TH1F("hTrueStdDevPhi_Neg", "Std Dev of True Phi Distortion Model of All Events, Negative Side (R > 30); #Delta Phi (#mum)", nsumbins, minsum, maxsum);

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
    //sourcefilename=((TFileInfo*)(sourcefilelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();
    sourcefilename=Form("/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.file%d.h_Charge_%d.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root",ifile,ihist);
    
    //create shifter
    shifter = new Shifter(sourcefilename, averagefilename);

    hFluctCharge=(TH3F*)fullcharge->Get(Form("h_Charge_%d",ihist)); // only 0-9 available
    //nbins: x 360, y 159, z 248
    //xmin: 0.0000000, 0.20000000, -1.0550000
    //xmax: 6.2831900, 0.78000000, 1.0550000
    hFluctCharge->Add(hSmoothedAve,-1); 
    
    TFile *plots;

    plots=TFile::Open(Form("CMModelsPhiRFull_Event%d.root",ifile),"READ");

    TH3F *hCartCMModelPhiRPos[3];
    hCartCMModelPhiRPos[0]=(TH3F*)plots->Get("hCMModelX_PhiR_Pos");
    hCartCMModelPhiRPos[1]=(TH3F*)plots->Get("hCMModelY_PhiR_Pos");
    hCartCMModelPhiRPos[2]=(TH3F*)plots->Get("hCMModelZ_PhiR_Pos");

    TH3F *hCylCMModelPhiRPos[2];
    hCylCMModelPhiRPos[0]=(TH3F*)plots->Get("hCMModelR_PhiR_Pos");
    hCylCMModelPhiRPos[1]=(TH3F*)plots->Get("hCMModelPhi_PhiR_Pos");

    TH3F *hCartCMModelPhiRNeg[3];
    hCartCMModelPhiRNeg[0]=(TH3F*)plots->Get("hCMModelX_PhiR_Neg");
    hCartCMModelPhiRNeg[1]=(TH3F*)plots->Get("hCMModelY_PhiR_Neg");
    hCartCMModelPhiRNeg[2]=(TH3F*)plots->Get("hCMModelZ_PhiR_Neg");

    TH3F *hCylCMModelPhiRNeg[2];
    hCylCMModelPhiRNeg[0]=(TH3F*)plots->Get("hCMModelR_PhiR_Neg");
    hCylCMModelPhiRNeg[1]=(TH3F*)plots->Get("hCMModelPhi_PhiR_Neg");
    
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

    // WriteIntFluctFile(ifile, ihist, nphi,   minphi,   maxphi,   nr,  minr,  maxr,   nz,   minzPos,  maxzPos,   minzNeg,  maxzNeg, hFluctCharge);

    //return 0;
    
    //positive
    TH1F *hCartesianShiftDifferencePhiRPos[3];
    hCartesianShiftDifferencePhiRPos[0] = new TH1F("hShiftDifferenceX_PhiR_Pos", "Difference between CM Model X and True, Phi,R binning, Positive Side (R > 30); #Delta X (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifferencePhiRPos[1] = new TH1F("hShiftDifferenceY_PhiR_Pos", "Difference between CM Model Y and True, Phi,R binning, Positive Side (R > 30); #Delta Y (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifferencePhiRPos[2] = new TH1F("hShiftDifferenceZ_PhiR_Pos", "Difference between CM Model Z and True, Phi,R binning, Positive Side (R > 30); #Delta Z (#mum)", ndiff, mindiff, maxdiff);
    
    TH1F *hCylindricalShiftDifferencePhiRPos[2];
    hCylindricalShiftDifferencePhiRPos[0] = new TH1F("hShiftDifferenceR_PhiR_Pos", "Difference between CM Model R and True, Phi,R binning, Positive Side (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    hCylindricalShiftDifferencePhiRPos[1] = new TH1F("hShiftDifferencePhi_PhiR_Pos", "Difference between CM Model Phi and True, Phi,R binning, Positive Side (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);

    TH1F *hRShiftTruePos = new TH1F("hRShiftTruePos", "True R Distortion Model, Positive Side (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    TH1F *hPhiShiftTruePos = new TH1F("hPhiShiftTruePos", "True Phi Distortion Model, Positive Side (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);
  
     TH2F *hCartesianDiffPhiRPos[6];
    hCartesianDiffPhiRPos[0] = new TH2F("hDiffXYX_PhiR_Pos", "Difference in PhiR for CM Model X, Phi,R binning, Positive Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianDiffPhiRPos[1] = new TH2F("hDiffRZX_PhiR_Pos", "Difference in RZ for CM Model X, Phi,R binning, Positive Side; z (cm); r (cm)", nz,minzPos,maxzPos,nr,minr,maxr);
    hCartesianDiffPhiRPos[2] = new TH2F("hDiffXYY_PhiR_Pos", "Difference in PhiR for CM Model Y, Phi,R binning, Positive Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianDiffPhiRPos[3] = new TH2F("hDiffRZY_PhiR_Pos", "Difference in RZ for CM Model Y, Phi,R binning, Positive Side; z (cm); r (cm)", nz,minzPos,maxzPos,nr,minr,maxr);
    hCartesianDiffPhiRPos[4] = new TH2F("hDiffXYZ_PhiR_Pos", "Difference in PhiR for CM Model Z, Phi,R binning, Positive Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianDiffPhiRPos[5] = new TH2F("hDiffRZZ_PhiR_Pos", "Difference in RZ for CM Model Z, Phi,R binning, Positive Side; z (cm); r (cm)", nz,minzPos,maxzPos,nr,minr,maxr);
    
    TH2F *hCylindricalDiffPhiRPos[4];
    hCylindricalDiffPhiRPos[0] = new TH2F("hDiffXYR_PhiR_Pos", "Difference in PhiR for CM Model R, Phi,R binning, Positive Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalDiffPhiRPos[1] = new TH2F("hDiffRZR_PhiR_Pos", "Difference in RZ for CM Model R, Phi,R binning, Positive Side; z (cm); r (cm)",nz,minzPos,maxzPos,nr,minr,maxr);
    hCylindricalDiffPhiRPos[2] = new TH2F("hDiffXYPhi_PhiR_Pos", "Difference in PhiR for CM Model Phi, Phi,R binning, Positive Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalDiffPhiRPos[3] = new TH2F("hDiffRZPhi_PhiR_Pos", "Difference in RZ for CM Model Phi, Phi,R binning, Positive Side; z (cm); r (cm)",nz,minzPos,maxzPos,nr,minr,maxr);
  
    TH2F *hCartesianAveDiffPhiRPos[6];
    hCartesianAveDiffPhiRPos[0] = new TH2F("hAveDiffXYX_PhiR_Pos", "X Model - Truth Averaged Over z, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianAveDiffPhiRPos[1] = new TH2F("hAveDiffRZX_PhiR_Pos", "X Model - Truth Averaged Over phi, Phi,R binning, Positive Side (#mum); z (cm); r (cm)", nz,minzPos,maxzPos,nr,minr,maxr);
    hCartesianAveDiffPhiRPos[2] = new TH2F("hAveDiffXYY_PhiR_Pos", "Y Model - Truth Averaged Over z, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianAveDiffPhiRPos[3] = new TH2F("hAveDiffRZY_PhiR_Pos", "Y Model - Truth Averaged Over phi, Phi,R binning, Positive Side (#mum); z (cm); r (cm)", nz,minzPos,maxzPos,nr,minr,maxr);
    hCartesianAveDiffPhiRPos[4] = new TH2F("hAveDiffXYZ_PhiR_Pos", "Z Model - Truth Averaged Over z, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianAveDiffPhiRPos[5] = new TH2F("hAveDiffRZZ_PhiR_Pos", "Z Model - Truth Averaged Over phi, Phi,R binning, Positive Side (#mum); z (cm); r (cm)", nz,minzPos,maxzPos,nr,minr,maxr);
    
     TH2F *hCylindricalAveDiffPhiRPos[4];
    hCylindricalAveDiffPhiRPos[0] = new TH2F("hAveDiffXYR_PhiR_Pos", "R Model - Truth Averaged Over z, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalAveDiffPhiRPos[1] = new TH2F("hAveDiffRZR_PhiR_Pos", "R Model - Truth Averaged Over phi, Phi,R binning, Positive Side (#mum); z (cm); r (cm)", nz,minzPos,maxzPos,nr,minr,maxr);
    hCylindricalAveDiffPhiRPos[2] = new TH2F("hAveDiffXYPhi_PhiR_Pos", "Phi Model - Truth Averaged Over z, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalAveDiffPhiRPos[3] = new TH2F("hAveDiffRZPhi_PhiR_Pos", "Phi Model - Truth Averaged Over phi, Phi,R binning, Positive Side (#mum); z (cm); r (cm)", nz,minzPos,maxzPos,nr,minr,maxr);

    TH2F *hSamplePerBinRZPos = new TH2F("hSamplePerBinRZPos", "Filling each rz bin, Positive Side; z (cm); r (cm)", nz,minzPos,maxzPos,nr,minr,maxr);

    TH2F *hSamplePerBinPhiR = new TH2F("hSamplePerBinPhiR", "Filling each PhiR bin; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);

    TH2F *hCompareRTrue_PhiRPos = new TH2F("hCompareRTrue_PhiR_Pos", "Compare Difference from R Model and True, Phi,R binning, Positive Side (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);
    TH2F *hComparePhiTrue_PhiRPos = new TH2F("hComparePhiTrue_PhiR_Pos", "Compare Difference from Phi Model and True, Phi,R binning, Positive Side (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);

    TH2F *hRDiffvR_PhiRPos = new TH2F("hRDiffvR_PhiR_Pos", "Difference between R Model and True vs. r, Phi,R binning, Positive Side (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvZ_PhiRPos = new TH2F("hRDiffvZ_PhiR_Pos", "Difference between R Model and True vs. z, Phi,R binning, Positive Side (R > 30); z (cm); shift difference (#mum)",nz,minzPos,maxzPos,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvPhi_PhiRPos = new TH2F("hRDiffvPhi_PhiR_Pos", "Difference between R Model and True vs. phi, Phi,R binning, Positive Side (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    TH2F *hPhiDiffvR_PhiRPos = new TH2F("hPhiDiffvR_PhiR_Pos", "Difference between Phi Model and True vs. r, Phi,R binning, Positive Side (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvZ_PhiRPos = new TH2F("hPhiDiffvZ_PhiR_Pos", "Difference between Phi Model and True vs. z, Phi,R binning, Positive Side (R > 30); z (cm); shift difference (#mum)",nz,minzPos,maxzPos,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvPhi_PhiRPos = new TH2F("hPhiDiffvPhi_PhiR_Pos", "Difference between Phi Model and True vs. phi, Phi,R binning, Positive Side (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    //Negative
    TH1F *hCartesianShiftDifferencePhiRNeg[3];
    hCartesianShiftDifferencePhiRNeg[0] = new TH1F("hShiftDifferenceX_PhiR_Neg", "Difference between CM Model X and True, Phi,R binning, Negative Side (R > 30); #Delta X (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifferencePhiRNeg[1] = new TH1F("hShiftDifferenceY_PhiR_Neg", "Difference between CM Model Y and True, Phi,R binning, Negative Side (R > 30); #Delta Y (#mum)", ndiff, mindiff, maxdiff);
    hCartesianShiftDifferencePhiRNeg[2] = new TH1F("hShiftDifferenceZ_PhiR_Neg", "Difference between CM Model Z and True, Phi,R binning, Negative Side (R > 30); #Delta Z (#mum)", ndiff, mindiff, maxdiff);
    
    TH1F *hCylindricalShiftDifferencePhiRNeg[2];
    hCylindricalShiftDifferencePhiRNeg[0] = new TH1F("hShiftDifferenceR_PhiR_Neg", "Difference between CM Model R and True, Phi,R binning, Negative Side (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    hCylindricalShiftDifferencePhiRNeg[1] = new TH1F("hShiftDifferencePhi_PhiR_Neg", "Difference between CM Model Phi and True, Phi,R binning, Negative Side (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);

    TH1F *hRShiftTrueNeg = new TH1F("hRShiftTrueNeg", "True R Distortion Model, Negative Side (R > 30); #Delta R (#mum)", ndiff, mindiff, maxdiff);
    TH1F *hPhiShiftTrueNeg = new TH1F("hPhiShiftTrueNeg", "True Phi Distortion Model, Negative Side (R > 30); #Delta Phi (#mum)", ndiff, mindiff, maxdiff);
  
     TH2F *hCartesianDiffPhiRNeg[6];
    hCartesianDiffPhiRNeg[0] = new TH2F("hDiffXYX_PhiR_Neg", "Difference in PhiR for CM Model X, Phi,R binning, Negative Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianDiffPhiRNeg[1] = new TH2F("hDiffRZX_PhiR_Neg", "Difference in RZ for CM Model X, Phi,R binning, Negative Side; z (cm); r (cm)", nz,minzNeg,maxzNeg,nr,minr,maxr);
    hCartesianDiffPhiRNeg[2] = new TH2F("hDiffXYY_PhiR_Neg", "Difference in PhiR for CM Model Y, Phi,R binning, Negative Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianDiffPhiRNeg[3] = new TH2F("hDiffRZY_PhiR_Neg", "Difference in RZ for CM Model Y, Phi,R binning, Negative Side; z (cm); r (cm)", nz,minzNeg,maxzNeg,nr,minr,maxr);
    hCartesianDiffPhiRNeg[4] = new TH2F("hDiffXYZ_PhiR_Neg", "Difference in PhiR for CM Model Z, Phi,R binning, Negative Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianDiffPhiRNeg[5] = new TH2F("hDiffRZZ_PhiR_Neg", "Difference in RZ for CM Model Z, Phi,R binning, Negative Side; z (cm); r (cm)", nz,minzNeg,maxzNeg,nr,minr,maxr);
    
    TH2F *hCylindricalDiffPhiRNeg[4];
    hCylindricalDiffPhiRNeg[0] = new TH2F("hDiffXYR_PhiR_Neg", "Difference in PhiR for CM Model R, Phi,R binning, Negative Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalDiffPhiRNeg[1] = new TH2F("hDiffRZR_PhiR_Neg", "Difference in RZ for CM Model R, Phi,R binning, Negative Side; z (cm); r (cm)",nz,minzNeg,maxzNeg,nr,minr,maxr);
    hCylindricalDiffPhiRNeg[2] = new TH2F("hDiffXYPhi_PhiR_Neg", "Difference in PhiR for CM Model Phi, Phi,R binning, Negative Side; phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalDiffPhiRNeg[3] = new TH2F("hDiffRZPhi_PhiR_Neg", "Difference in RZ for CM Model Phi, Phi,R binning, Negative Side; z (cm); r (cm)",nz,minzNeg,maxzNeg,nr,minr,maxr);
  
    TH2F *hCartesianAveDiffPhiRNeg[6];
    hCartesianAveDiffPhiRNeg[0] = new TH2F("hAveDiffXYX_PhiR_Neg", "X Model - Truth Averaged Over z, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianAveDiffPhiRNeg[1] = new TH2F("hAveDiffRZX_PhiR_Neg", "X Model - Truth Averaged Over phi, Phi,R binning, Negative Side (#mum); z (cm); r (cm)", nz,minzNeg,maxzNeg,nr,minr,maxr);
    hCartesianAveDiffPhiRNeg[2] = new TH2F("hAveDiffXYY_PhiR_Neg", "Y Model - Truth Averaged Over z, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianAveDiffPhiRNeg[3] = new TH2F("hAveDiffRZY_PhiR_Neg", "Y Model - Truth Averaged Over phi, Phi,R binning, Negative Side (#mum); z (cm); r (cm)", nz,minzNeg,maxzNeg,nr,minr,maxr);
    hCartesianAveDiffPhiRNeg[4] = new TH2F("hAveDiffXYZ_PhiR_Neg", "Z Model - Truth Averaged Over z, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCartesianAveDiffPhiRNeg[5] = new TH2F("hAveDiffRZZ_PhiR_Neg", "Z Model - Truth Averaged Over phi, Phi,R binning, Negative Side (#mum); z (cm); r (cm)", nz,minzNeg,maxzNeg,nr,minr,maxr);
    
     TH2F *hCylindricalAveDiffPhiRNeg[4];
    hCylindricalAveDiffPhiRNeg[0] = new TH2F("hAveDiffXYR_PhiR_Neg", "R Model - Truth Averaged Over z, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalAveDiffPhiRNeg[1] = new TH2F("hAveDiffRZR_PhiR_Neg", "R Model - Truth Averaged Over phi, Phi,R binning, Negative Side (#mum); z (cm); r (cm)", nz,minzNeg,maxzNeg,nr,minr,maxr);
    hCylindricalAveDiffPhiRNeg[2] = new TH2F("hAveDiffXYPhi_PhiR_Neg", "Phi Model - Truth Averaged Over z, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nphi,minphi,maxphi,nr,minr,maxr);
    hCylindricalAveDiffPhiRNeg[3] = new TH2F("hAveDiffRZPhi_PhiR_Neg", "Phi Model - Truth Averaged Over phi, Phi,R binning, Negative Side (#mum); z (cm); r (cm)", nz,minzNeg,maxzNeg,nr,minr,maxr);

    TH2F *hSamplePerBinRZNeg = new TH2F("hSamplePerBinRZNeg", "Filling each rz bin, Negative Side; z (cm); r (cm)", nz,minzNeg,maxzNeg,nr,minr,maxr);

    TH2F *hCompareRTrue_PhiRNeg = new TH2F("hCompareRTrue_PhiR_Neg", "Compare Difference from R Model and True, Phi,R binning, Negative Side (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);
    TH2F *hComparePhiTrue_PhiRNeg = new TH2F("hComparePhiTrue_PhiR_Neg", "Compare Difference from Phi Model and True, Phi,R binning, Negative Side (R > 30, 10 < z < 90); reco shift (#mum); true shift (#mum)",nbins,-550,550,nbins,-550,550);

    TH2F *hRDiffvR_PhiRNeg = new TH2F("hRDiffvR_PhiR_Neg", "Difference between R Model and True vs. r, Phi,R binning, Negative Side (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvZ_PhiRNeg = new TH2F("hRDiffvZ_PhiR_Neg", "Difference between R Model and True vs. z, Phi,R binning, Negative Side (R > 30); z (cm); shift difference (#mum)",nz,minzNeg,maxzNeg,ndiff,mindiff,maxdiff);
    TH2F *hRDiffvPhi_PhiRNeg = new TH2F("hRDiffvPhi_PhiR_Neg", "Difference between R Model and True vs. phi, Phi,R binning, Negative Side (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);

    TH2F *hPhiDiffvR_PhiRNeg = new TH2F("hPhiDiffvR_PhiR_Neg", "Difference between Phi Model and True vs. r, Phi,R binning, Negative Side (R > 30, 10 < z < 90); r (cm); shift difference (#mum)",nr,minr,maxr,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvZ_PhiRNeg = new TH2F("hPhiDiffvZ_PhiR_Neg", "Difference between Phi Model and True vs. z, Phi,R binning, Negative Side (R > 30); z (cm); shift difference (#mum)",nz,minzNeg,maxzNeg,ndiff,mindiff,maxdiff);
    TH2F *hPhiDiffvPhi_PhiRNeg = new TH2F("hPhiDiffvPhi_PhiR_Neg", "Difference between Phi Model and True vs. phi, Phi,R binning, Negative Side (R > 30, 10 < z < 90); phi (rad); shift difference (#mum)",nphi,minphi,maxphi,ndiff,mindiff,maxdiff);
    
    //compare linear model to space charge
    TH2F *hCompareRTruevFluctNeg = new TH2F("hCompareRTruevFluct", "Compare True R Distortion Fluctuation and True Charge Fluctuation, Phi,R binning, Negative Side (R > 30); fluct charge (#mum); true shift (#mum)",nbins,-1e4,1e4,nbins,-30,30);
    TH2F *hCompareRDiffvFluctNeg = new TH2F("hCompareRDiffvFluct", "Compare Difference between R Model and True R vs True Charge Fluctuation, Phi,R binning, Negative Side (R > 30); fluct charge (#mum); shift difference (#mum)",nbins,-1e4,1e4,nbins,-30,30);

    //compare linear model to integrated charge
    TH2F *hCompareRTruevIntFluctSmallRPos = new TH2F("hCompareRTruevIntFluctSmallRPos", "Compare True R Distortion Fluctuation and True Integrated Charge Fluctuation, Positive Side, R Inside; int fluct charge (#mum); true shift (#mum)",nbins,-1e6,1e6,nbins,-3e3,3e3);
    TH2F *hCompareRTruevIntFluctLargeRPos = new TH2F("hCompareRTruevIntFluctLargeRPos", "Compare True R Distortion Fluctuation and True Integrated Charge Fluctuation, Positive Side, R Outside; int fluct charge (#mum); true shift (#mum)",nbins,-1e6,1e6,nbins,-3e3,3e3);
    TH2F *hCompareRTruevIntFluctDiffRPos = new TH2F("hCompareRTruevIntFluctDiffRPos", "Compare True R Distortion Fluctuation and True Integrated Charge Fluctuation, Positive Side, Difference of R In and Out; int fluct charge (#mum); true shift (#mum)",nbins,-1e6,1e6,nbins,-3e3,3e3);

    TH2F *hCompareRDiffvIntFluctSmallRPos = new TH2F("hCompareRDiffvIntFluctSmallRPos", "Compare Difference between R Model and True R vs True Integrated Charge Fluctuation, Positive Side, R Inside; int fluct charge (#mum); shift difference (#mum)",nbins,-1e6,1e6,nbins,-3e3,3e3);
    TH2F *hCompareRDiffvIntFluctLargeRPos = new TH2F("hCompareRDiffvIntFluctLargeRPos", "Compare Difference between R Model and True R vs True Integrated Charge Fluctuation, Positive Side, R Outside; int fluct charge (#mum); shift difference (#mum)",nbins,-1e6,1e6,nbins,-3e3,3e3);
    TH2F *hCompareRDiffvIntFluctDiffRPos = new TH2F("hCompareRDiffvIntFluctDiffRPos", "Compare Difference between R Model and True R vs True Integrated Charge Fluctuation, Positive Side, Difference of R In and Out; int fluct charge (#mum); shift difference (#mum)",nbins,-1e6,1e6,nbins,-3e3,3e3);

  
    //TH1F *hFluc = new TH1F("hFluc", "Fluctuation Charge", 1000, 1, 1e7); 

    TFile *intFluct;
    TH3F *hIntFluctChargeSmallRPos, *hIntFluctChargeLargeRPos;
   
    intFluct=TFile::Open(Form("IntFluctEvent%d.root", (2*ifile + ihist)), "READ");//change 2 to 10 when looping over all ihist

    hIntFluctChargeSmallRPos=(TH3F*)intFluct->Get("hIntFluctChargeSmallRPos");
    hIntFluctChargeLargeRPos=(TH3F*)intFluct->Get("hIntFluctChargeLargeRPos");
    				  
    for(int i = 1; i < nphi - 1; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
      for(int j = 1; j < nr - 1; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin
	for(int k = 1; k < nz - 1; k++){
	  double zPos = minzPos + ((maxzPos - minzPos)/(1.0*nz))*(k+0.5); //center of bin
	  double zNeg = minzNeg + ((maxzNeg - minzNeg)/(1.0*nz))*(k+0.5); //center of bin

	  //positive
	  double shifttrueCartPos[3];
	  double shifttrueCylPos[2];

	  double shiftrecoCartPhiRPos[3];
	  double differenceCartPhiRPos[3];

	  double shiftrecoCylPhiRPos[2];
	  double differenceCylPhiRPos[2];

	  double differenceR_PhiRPos, differencePhi_PhiRPos;	  

	  int binPhiRPos = hCartCMModelPhiRPos[0]->FindBin(phi,r,zPos);

	  
	  if((r > 30.0) && (r < 76.0)){
	    //x y and z
	    shifttrueCartPos[0] = (shifter->hPosX->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
	    shifttrueCartPos[1] = (shifter->hPosY->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron 
	    shifttrueCartPos[2] = (shifter->hPosZ->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
	    //r and phi
	    shifttrueCylPos[0] = (shifter->hPosR->Interpolate(phi,r,zPos))*(1e4); //convert from cm to micron
	    shifttrueCylPos[1] = (shifter->hPosPhi->Interpolate(phi,r,zPos))*(1e4);
	    hRShiftTruePos->Fill(shifttrueCylPos[0]);
	    hPhiShiftTruePos->Fill(shifttrueCylPos[1]);
	    
	    for(int l = 0; l < 3; l ++){
	      shiftrecoCartPhiRPos[l] =  (hCartCMModelPhiRPos[l]->GetBinContent(binPhiRPos))*(1e4);
	  
	      differenceCartPhiRPos[l] = shiftrecoCartPhiRPos[l] - shifttrueCartPos[l]; 

	      hCartesianShiftDifferencePhiRPos[l]->Fill(differenceCartPhiRPos[l]);
	    }

	    //r
	    shiftrecoCylPhiRPos[0] =  (hCylCMModelPhiRPos[0]->GetBinContent(binPhiRPos))*(1e4);
	    differenceCylPhiRPos[0] = shiftrecoCylPhiRPos[0] - shifttrueCylPos[0]; 
	    hCylindricalShiftDifferencePhiRPos[0]->Fill(differenceCylPhiRPos[0]);
	      
	    //phi
	    shiftrecoCylPhiRPos[1] = r*(1e4)*(hCylCMModelPhiRPos[1]->GetBinContent(binPhiRPos));
	    differenceCylPhiRPos[1] = (shiftrecoCylPhiRPos[1] - shifttrueCylPos[1]); 
	    hCylindricalShiftDifferencePhiRPos[1]->Fill(differenceCylPhiRPos[1]);

	    //x
	    hCartesianDiffPhiRPos[0]->Fill(phi,r, differenceCartPhiRPos[0]);
	    hCartesianDiffPhiRPos[1]->Fill(zPos,r, differenceCartPhiRPos[0]);
	    //y
	    hCartesianDiffPhiRPos[2]->Fill(phi,r, differenceCartPhiRPos[1]);	  
	    hCartesianDiffPhiRPos[3]->Fill(zPos,r, differenceCartPhiRPos[1]);
	    //z
	    hCartesianDiffPhiRPos[4]->Fill(phi,r, differenceCartPhiRPos[2]);
	    hCartesianDiffPhiRPos[5]->Fill(zPos,r, differenceCartPhiRPos[2]);

	    //r
	    hCylindricalDiffPhiRPos[0]->Fill(phi,r, differenceCylPhiRPos[0]);
	    hCylindricalDiffPhiRPos[1]->Fill(zPos,r, differenceCylPhiRPos[0]);

	    hCompareRTrue_PhiRPos->Fill(shiftrecoCylPhiRPos[0],shifttrueCylPos[0]);

	    hRDiffvR_PhiRPos->Fill(r,differenceCylPhiRPos[0],1);
	    hRDiffvPhi_PhiRPos->Fill(phi,differenceCylPhiRPos[0],1);
	    hRDiffvZ_PhiRPos->Fill(zPos,differenceCylPhiRPos[0],1);

	    //phi 
	    hCylindricalDiffPhiRPos[2]->Fill(phi,r, differenceCylPhiRPos[1]);
	    hCylindricalDiffPhiRPos[3]->Fill(zPos,r, differenceCylPhiRPos[1]);
	    	    
	    hComparePhiTrue_PhiRPos->Fill(shiftrecoCylPhiRPos[1],shifttrueCylPos[1]);
	    
	    hPhiDiffvR_PhiRPos->Fill(r,differenceCylPhiRPos[1],1);
	    hPhiDiffvPhi_PhiRPos->Fill(phi,differenceCylPhiRPos[1],1);
	    hPhiDiffvZ_PhiRPos->Fill(zPos,differenceCylPhiRPos[1],1);

	    hSamplePerBinRZPos->Fill(zPos,r,1);

	    //fluct
	    double intfluctchargeSmallRPos = hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
	    double intfluctchargeLargeRPos = hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
	    double intfluctchargeDiffRPos = intfluctchargeSmallRPos - intfluctchargeLargeRPos;

	    hCompareRTruevIntFluctSmallRPos->Fill(intfluctchargeSmallRPos,shifttrueCylPos[0]);
	    hCompareRTruevIntFluctLargeRPos->Fill(intfluctchargeLargeRPos,shifttrueCylPos[0]);
	    hCompareRTruevIntFluctDiffRPos->Fill(intfluctchargeDiffRPos,shifttrueCylPos[0]);

	    hCompareRDiffvIntFluctSmallRPos->Fill(intfluctchargeSmallRPos,differenceCylPhiRPos[0]);
	    hCompareRDiffvIntFluctLargeRPos->Fill(intfluctchargeLargeRPos,differenceCylPhiRPos[0]);
	    hCompareRDiffvIntFluctDiffRPos->Fill(intfluctchargeDiffRPos,differenceCylPhiRPos[0]);
	  }
	    
	  //negative
	  double shifttrueCartNeg[3];
	  double shifttrueCylNeg[2];
	  
	  double shiftrecoCartPhiRNeg[3];
	  double differenceCartPhiRNeg[3];

	  double shiftrecoCylPhiRNeg[2];
	  double differenceCylPhiRNeg[2];

	  double differenceR_PhiRNeg, differencePhi_PhiRNeg;

	  double fluctchargeNeg;

	  int binPhiRNeg = hCartCMModelPhiRNeg[0]->FindBin(phi,r,zNeg);

	  if((r > 30.0) && (r < 76.0)){
	    //x y and z
	    shifttrueCartNeg[0] = (shifter->hNegX->Interpolate(phi,r,zNeg))*(1e4); //convert from cm to micron
	    shifttrueCartNeg[1] = (shifter->hNegY->Interpolate(phi,r,zNeg))*(1e4); //convert from cm to micron 
	    shifttrueCartNeg[2] = (shifter->hNegZ->Interpolate(phi,r,zNeg))*(1e4); //convert from cm to micron
	    //r and phi
	    shifttrueCylNeg[0] = (shifter->hNegR->Interpolate(phi,r,zNeg))*(1e4); //convert from cm to micron
	    shifttrueCylNeg[1] = (shifter->hNegPhi->Interpolate(phi,r,zNeg))*(1e4);
	    hRShiftTrueNeg->Fill(shifttrueCylNeg[0]);
	    hPhiShiftTrueNeg->Fill(shifttrueCylNeg[1]);

	    
	    fluctchargeNeg = (hFluctCharge->Interpolate(phi,r/100.,zNeg/100.));//convert from cm to micron
	   	    
	    for(int l = 0; l < 3; l ++){
	      shiftrecoCartPhiRNeg[l] =  (hCartCMModelPhiRNeg[l]->GetBinContent(binPhiRNeg))*(1e4);
	  
	      differenceCartPhiRNeg[l] = shiftrecoCartPhiRNeg[l] - shifttrueCartNeg[l]; 

	      hCartesianShiftDifferencePhiRNeg[l]->Fill(differenceCartPhiRNeg[l]);
	    }

	    //r
	    shiftrecoCylPhiRNeg[0] =  (hCylCMModelPhiRNeg[0]->GetBinContent(binPhiRNeg))*(1e4);
	    differenceCylPhiRNeg[0] = shiftrecoCylPhiRNeg[0] - shifttrueCylNeg[0]; 
	    hCylindricalShiftDifferencePhiRNeg[0]->Fill(differenceCylPhiRNeg[0]);
	      
	    //phi
	    shiftrecoCylPhiRNeg[1] = r*(1e4)*(hCylCMModelPhiRNeg[1]->GetBinContent(binPhiRNeg));
	    differenceCylPhiRNeg[1] = (shiftrecoCylPhiRNeg[1] - shifttrueCylNeg[1]); 
	    hCylindricalShiftDifferencePhiRNeg[1]->Fill(differenceCylPhiRNeg[1]);

	    //x
	    hCartesianDiffPhiRNeg[0]->Fill(phi,r, differenceCartPhiRNeg[0]);
	    hCartesianDiffPhiRNeg[1]->Fill(zNeg,r, differenceCartPhiRNeg[0]);
	    //y
	    hCartesianDiffPhiRNeg[2]->Fill(phi,r, differenceCartPhiRNeg[1]);	  
	    hCartesianDiffPhiRNeg[3]->Fill(zNeg,r, differenceCartPhiRNeg[1]);
	    //z
	    hCartesianDiffPhiRNeg[4]->Fill(phi,r, differenceCartPhiRNeg[2]);
	    hCartesianDiffPhiRNeg[5]->Fill(zNeg,r, differenceCartPhiRNeg[2]);

	    //r
	    hCylindricalDiffPhiRNeg[0]->Fill(phi,r, differenceCylPhiRNeg[0]);
	    hCylindricalDiffPhiRNeg[1]->Fill(zNeg,r, differenceCylPhiRNeg[0]);

	    hCompareRTrue_PhiRNeg->Fill(shiftrecoCylPhiRNeg[0],shifttrueCylNeg[0]);

	    hRDiffvR_PhiRNeg->Fill(r,differenceCylPhiRNeg[0],1);
	    hRDiffvPhi_PhiRNeg->Fill(phi,differenceCylPhiRNeg[0],1);
	    hRDiffvZ_PhiRNeg->Fill(zNeg,differenceCylPhiRNeg[0],1);

	    //phi 
	    hCylindricalDiffPhiRNeg[2]->Fill(phi,r, differenceCylPhiRNeg[1]);
	    hCylindricalDiffPhiRNeg[3]->Fill(zNeg,r, differenceCylPhiRNeg[1]);
	    	    
	    hComparePhiTrue_PhiRNeg->Fill(shiftrecoCylPhiRNeg[1],shifttrueCylNeg[1]);
	    
	    hPhiDiffvR_PhiRNeg->Fill(r,differenceCylPhiRNeg[1],1);
	    hPhiDiffvPhi_PhiRNeg->Fill(phi,differenceCylPhiRNeg[1],1);
	    hPhiDiffvZ_PhiRNeg->Fill(zNeg,differenceCylPhiRNeg[1],1);

	    hSamplePerBinRZNeg->Fill(zNeg,r,1);
	    
	    hSamplePerBinPhiR->Fill(phi,r,1);

	    //fluct
	    hCompareRTruevFluctNeg->Fill(fluctchargeNeg,shifttrueCylNeg[0]);
	    hCompareRDiffvFluctNeg->Fill(fluctchargeNeg,differenceCylPhiRNeg[0]);
	    
	  }
	}
      }
    }
  
    //average over z
    for (int m = 0; m < 6; m = m+2){
      hCartesianAveDiffPhiRPos[m]->Divide(hCartesianDiffPhiRPos[m],hSamplePerBinPhiR);
      hCartesianAveDiffPhiRNeg[m]->Divide(hCartesianDiffPhiRNeg[m],hSamplePerBinPhiR);
    }
    for (int m = 0; m < 4; m = m+2){
      hCylindricalAveDiffPhiRPos[m]->Divide(hCylindricalDiffPhiRPos[m],hSamplePerBinPhiR);
      hCylindricalAveDiffPhiRNeg[m]->Divide(hCylindricalDiffPhiRNeg[m],hSamplePerBinPhiR);
    }
    
    //average over phi
    for (int m = 1; m < 6; m = m+2){
      hCartesianAveDiffPhiRPos[m]->Divide(hCartesianDiffPhiRPos[m],hSamplePerBinRZPos);
      hCartesianAveDiffPhiRNeg[m]->Divide(hCartesianDiffPhiRNeg[m],hSamplePerBinRZNeg);
    }
    for (int m = 1; m < 4; m = m+2){
      hCylindricalAveDiffPhiRPos[m]->Divide(hCylindricalDiffPhiRPos[m],hSamplePerBinRZPos);
      hCylindricalAveDiffPhiRNeg[m]->Divide(hCylindricalDiffPhiRNeg[m],hSamplePerBinRZNeg);
    }

      
    //summary plots
    //positive
    hDifferenceMeanRPos->Fill(hCylindricalShiftDifferencePhiRPos[0]->GetMean(1));
    hDifferenceStdDevRPos->Fill(hCylindricalShiftDifferencePhiRPos[0]->GetStdDev(1));

    hTrueMeanRPos->Fill(hRShiftTruePos->GetMean(1));
    hTrueStdDevRPos->Fill(hRShiftTruePos->GetStdDev(1));
    
    hDifferenceMeanPhiPos->Fill(hCylindricalShiftDifferencePhiRPos[1]->GetMean(1));
    hDifferenceStdDevPhiPos->Fill(hCylindricalShiftDifferencePhiRPos[1]->GetStdDev(1));

    hTrueMeanPhiPos->Fill(hPhiShiftTruePos->GetMean(1));
    hTrueStdDevPhiPos->Fill(hPhiShiftTruePos->GetStdDev(1));

    //negative
    hDifferenceMeanRNeg->Fill(hCylindricalShiftDifferencePhiRNeg[0]->GetMean(1));
    hDifferenceStdDevRNeg->Fill(hCylindricalShiftDifferencePhiRNeg[0]->GetStdDev(1));

    hTrueMeanRNeg->Fill(hRShiftTrueNeg->GetMean(1));
    hTrueStdDevRNeg->Fill(hRShiftTrueNeg->GetStdDev(1));
    
    hDifferenceMeanPhiNeg->Fill(hCylindricalShiftDifferencePhiRNeg[1]->GetMean(1));
    hDifferenceStdDevPhiNeg->Fill(hCylindricalShiftDifferencePhiRNeg[1]->GetStdDev(1));

    hTrueMeanPhiNeg->Fill(hPhiShiftTrueNeg->GetMean(1));
    hTrueStdDevPhiNeg->Fill(hPhiShiftTrueNeg->GetStdDev(1));

    //remove stat boxes
    for (int m = 0; m < 6; m++){
      hCartesianAveDiffPhiRPos[m]->SetStats(0);
      hCartesianAveDiffPhiRNeg[m]->SetStats(0);
    }
    for (int m = 0; m < 4; m++){
      hCylindricalAveDiffPhiRPos[m]->SetStats(0);
      hCylindricalAveDiffPhiRNeg[m]->SetStats(0);
    }

    //positive
    hCompareRTrue_PhiRPos->SetStats(0);
    hComparePhiTrue_PhiRPos->SetStats(0);

    hRDiffvR_PhiRPos->SetStats(0);
    hRDiffvZ_PhiRPos->SetStats(0);
    hRDiffvPhi_PhiRPos->SetStats(0);
  
    hPhiDiffvR_PhiRPos->SetStats(0);
    hPhiDiffvZ_PhiRPos->SetStats(0);
    hPhiDiffvPhi_PhiRPos->SetStats(0);

    //negative
    hCompareRTrue_PhiRNeg->SetStats(0);
    hComparePhiTrue_PhiRNeg->SetStats(0);

    hRDiffvR_PhiRNeg->SetStats(0);
    hRDiffvZ_PhiRNeg->SetStats(0);
    hRDiffvPhi_PhiRNeg->SetStats(0);
  
    hPhiDiffvR_PhiRNeg->SetStats(0);
    hPhiDiffvZ_PhiRNeg->SetStats(0);
    hPhiDiffvPhi_PhiRNeg->SetStats(0);

    hCompareRTruevFluctNeg->SetStats(0);
    hCompareRDiffvFluctNeg->SetStats(0);

    hIntFluctChargeSmallRPos->Project3D("yz")->SetStats(0);
    hIntFluctChargeLargeRPos->Project3D("yz")->SetStats(0);
    hFluctCharge->Project3D("yz")->SetStats(0);

    gStyle->SetOptStat("oue");
    
    /* hCompareRTruevIntFluctSmallRPos->SetStats(0);
    hCompareRTruevIntFluctLargeRPos->SetStats(0);
    hCompareRTruevIntFluctDiffRPos->SetStats(0);

    hCompareRDiffvIntFluctSmallRPos->SetStats(0);
    hCompareRDiffvIntFluctLargeRPos->SetStats(0);
    hCompareRDiffvIntFluctDiffRPos->SetStats(0);
    */
    
    
    //integrated fluct plots
    /* integcomp->Divide(3,2);
    integcomp->cd(1);
    hCompareRTruevIntFluctSmallRPos->Draw("colz");
    integcomp->cd(2);
    hCompareRTruevIntFluctLargeRPos->Draw("colz");
    integcomp->cd(3);
    hCompareRTruevIntFluctDiffRPos->Draw("colz");
    integcomp->cd(4);
    hCompareRDiffvIntFluctSmallRPos->Draw("colz");
    integcomp->cd(5);
    hCompareRDiffvIntFluctLargeRPos->Draw("colz");
    integcomp->cd(6);
    hCompareRDiffvIntFluctDiffRPos->Draw("colz");

    if((2*ifile + ihist) == 0){ 
      integcomp->Print("integcomp.pdf(","pdf");
    } else if ((2*ifile + ihist) == nEvents - 1){
      integcomp->Print("integcomp.pdf)","pdf");
    } else {
      integcomp->Print("integcomp.pdf","pdf");
      }*/
    
    TPad *integtitlepad = new TPad("integtitlepad","",0.0,0.96,1.0,1.0);
    TPad *integplots = new TPad("integplotspad","",0.0,0.0,1.0,0.96);

    TLatex *integtitle = new TLatex(0.0,0.0,"");

    integtitle->SetNDC();
    integtitle->SetTextSize(0.4);

    integ->cd();
    
    integplots->Draw();
    integtitlepad->Draw();
    
    integ->Divide(2,2);
    integ->cd(1);
    hIntFluctChargeSmallRPos->Project3D("yz")->Draw("colz");
    integ->cd(2);
    hIntFluctChargeLargeRPos->Project3D("yz")->Draw("colz");
    
    //hCompareRTruevIntFluctSmallRPos->Draw("colz");
    
    //hCompareRTruevIntFluctLargeRPos->Draw("colz");
    integ->cd(3);
    hFluctCharge->Project3D("yz")->Draw("colz");
    integ->cd(4)->Clear();

    integtitlepad->cd();
    integtitlepad->Clear();
    integtitle->DrawLatex(0.4,0.4,Form("Event %d", (2*ifile + ihist))); 

    if(ifile == 0){
      integ->Print("IntegratedFluctAnalysis.pdf(","pdf");
    } else if((ifile == 1) || (ifile == nEvents - 1)){
      integ->Print("IntegratedFluctAnalysis.pdf","pdf");
    } else{
      integ->Print("IntegratedFluctAnalysis.pdf)","pdf");
    }
    
    TPad *c1=new TPad("c1","",0.0,0.8,1.0,0.93); //.13 height each
    TPad *c2=new TPad("c2","",0.0,0.64,1.0,0.77);
    TPad *c3=new TPad("c3","",0.0,0.48,1.0,0.61);
    TPad *c4=new TPad("c4","",0.0,0.32,1.0,0.45);
    TPad *c5=new TPad("c5","",0.0,0.16,1.0,0.29);
    TPad *c6=new TPad("c6","",0.0,0.0,1.0,0.13);
    
    TPad *titlepad=new TPad("titlepad","",0.0,0.96,1.0,1.0); //0.04 height

    TPad *stitlepad1=new TPad("stitlepad1","",0.0,0.93,1.0,0.96); //0.03 height
    TPad *stitlepad2=new TPad("stitlepad2","",0.0,0.77,1.0,0.8);
    TPad *stitlepad3=new TPad("stitlepad3","",0.0,0.61,1.0,0.64);
    TPad *stitlepad4=new TPad("stitlepad4","",0.0,0.45,1.0,0.48);
    TPad *stitlepad5=new TPad("stitlepad5","",0.0,0.29,1.0,0.32);
    TPad *stitlepad6=new TPad("stitlepad6","",0.0,0.13,1.0,0.16);
    
    TLatex * title = new TLatex(0.0,0.0,"");
    title->SetNDC();
    title->SetTextSize(0.32);
    
    TLatex *stitle[6];
    for(int i = 0; i < 6; i++){
      stitle[i] = new TLatex(0.0,0.0,"");
      stitle[i]->SetNDC();
      stitle[i]->SetTextSize(0.35);
    }
    
    canvas->cd();
    c1->Draw();
    stitlepad1->Draw();
    c2->Draw();
    stitlepad2->Draw();
    c3->Draw();
    stitlepad3->Draw();
    c4->Draw();
    stitlepad4->Draw();
    c5->Draw();
    stitlepad5->Draw();
    c6->Draw();
    stitlepad6->Draw();
    titlepad->Draw();

    //x plots
    c1->Divide(4,1);
    c1->cd(1);
    hCartesianAveDiffPhiRNeg[0]->Draw("colz");
    c1->cd(2);
    // hCartesianAveDiffPhiRNeg[1]->Draw("colz");
    hCompareRTruevIntFluctSmallRPos->Draw("colz");
    c1->cd(3);
    //hCartesianShiftDifferencePhiRNeg[0]->Draw();
    hCompareRTruevIntFluctLargeRPos->Draw("colz");
    //c1->cd(4)->Clear();  
    c1->cd(4);
    //hCMmodelSliceRvTrue->Draw("colz");
    //hSamplePerBinRZNeg->Draw("colz");
    //hCompareRTruevFluctNeg->Draw("colz");
    //hFluc->Draw();
    hCompareRTruevIntFluctDiffRPos->Draw("colz");
    
    //y plots
    c2->Divide(4,1);
    c2->cd(1);
    hCartesianAveDiffPhiRNeg[2]->Draw("colz");
    c2->cd(2);
    //hCartesianAveDiffPhiRNeg[3]->Draw("colz");
    hCompareRDiffvIntFluctSmallRPos->Draw("colz");
    c2->cd(3);
    //hCartesianShiftDifferencePhiRNeg[1]->Draw();
    //c2->cd(4)->Clear();
    hCompareRDiffvIntFluctLargeRPos->Draw("colz");
    c2->cd(4);
    //hStripesPerBin->Draw("colz");
    //hSamplePerBinPhiR->Draw("colz");
    //hCompareRDiffvFluctNeg->Draw("colz");
    hCompareRDiffvIntFluctDiffRPos->Draw("colz");
    
    //r cart
    c3->Divide(4,1);
    c3->cd(1);
    hCylindricalAveDiffPhiRNeg[0]->Draw("colz");
    c3->cd(2);
    hCylindricalAveDiffPhiRNeg[1]->Draw("colz");
    c3->cd(3);
    hCylindricalShiftDifferencePhiRNeg[0]->Draw();
    c3->cd(4);
    hRShiftTrueNeg->Draw();
    
    //phi cart
    c4->Divide(4,1);
    c4->cd(1);
    hCylindricalAveDiffPhiRNeg[2]->Draw("colz");
    c4->cd(2);
    hCylindricalAveDiffPhiRNeg[3]->Draw("colz");
    c4->cd(3);
    hCylindricalShiftDifferencePhiRNeg[1]->Draw();
    c4->cd(4);
    hPhiShiftTrueNeg->Draw();

    //r to true comparison
    c5->Divide(4,1);
    c5->cd(1);
    hCompareRTrue_PhiRNeg->Draw("colz");
    c5->cd(2);
    hRDiffvR_PhiRNeg->Draw("colz");
    c5->cd(3);
    hRDiffvZ_PhiRNeg->Draw("colz");
    c5->cd(4);
    hRDiffvPhi_PhiRNeg->Draw("colz");

    //phi to true comparison
    c6->Divide(4,1);
    c6->cd(1);
    hComparePhiTrue_PhiRNeg->Draw("colz");
    c6->cd(2);
    hPhiDiffvR_PhiRNeg->Draw("colz");
    c6->cd(3);
    hPhiDiffvZ_PhiRNeg->Draw("colz");
    c6->cd(4);
    hPhiDiffvPhi_PhiRNeg->Draw("colz");

    titlepad->cd();
    titlepad->Clear();
    title->DrawLatex(0.01,0.4,Form("Event %d; %s", (2*ifile + ihist), sourcefilename.Data())); 
    title->Draw();
    
    stitlepad1->cd();
    stitlepad1->Clear();
    stitle[0]->DrawLatex(0.45,0.2,"X Model"); 
    stitle[0]->Draw();
     
    stitlepad2->cd();
    stitlepad2->Clear();
    stitle[1]->DrawLatex(0.45,0.2,"Y Model"); 
    stitle[1]->Draw();

    stitlepad3->cd();
    stitlepad3->Clear();
    stitle[2]->DrawLatex(0.45,0.2,"R Model"); 
    stitle[2]->Draw();

    stitlepad4->cd();
    stitlepad4->Clear();
    stitle[3]->DrawLatex(0.45,0.2,"Phi Model"); 
    stitle[3]->Draw();

    stitlepad5->cd();
    stitlepad5->Clear();
    stitle[4]->DrawLatex(0.4,0.2,"Comparing R Model to True"); 
    stitle[4]->Draw();

    stitlepad6->cd();
    stitlepad6->Clear();
    stitle[5]->DrawLatex(0.4,0.2,"Comparing Phi Model to True"); 
    stitle[5]->Draw();

    if(ifile == 0){ 
      //if(ifile == 1){
      canvas->Print("CMDistortionAnalysisPhiRNeg.pdf(","pdf");
    } else if((ifile == 1) || (ifile == nEvents - 1)){
      canvas->Print("CMDistortionAnalysisPhiRNeg.pdf","pdf");
    }
    }
  }

  TCanvas *summary = new TCanvas("summary","ShiftPlotsSummary",2000,3000);

  TPad *sumtitlepad = new TPad("sumtitlepad","",0.0,0.96,1.0,1.0);
  TPad *sumplots = new TPad("sumplotspad","",0.0,0.0,1.0,0.96);

  TLatex *sumtitle = new TLatex(0.0,0.0,"");

  sumtitle->SetNDC();
  sumtitle->SetTextSize(0.4);

  summary->cd();
  sumplots->Draw();
  sumtitlepad->Draw();

  sumplots->Divide(4,6);
  sumplots->cd(1);
  hDifferenceMeanRNeg->Draw();
  sumplots->cd(2);
  hDifferenceStdDevRNeg->Draw();
  sumplots->cd(3);
  hTrueMeanRNeg->Draw();
  sumplots->cd(4);
  hTrueStdDevRNeg->Draw();
  sumplots->cd(5);
  hDifferenceMeanPhiNeg->Draw();
  sumplots->cd(6);
  hDifferenceStdDevPhiNeg->Draw();
  sumplots->cd(7);
  hTrueMeanPhiNeg->Draw();
  sumplots->cd(8);
  hTrueStdDevPhiNeg->Draw();
  sumplots->cd(9);
  sumplots->cd(10)->Clear();
  sumplots->cd(11)->Clear();
  sumplots->cd(12)->Clear();
  sumplots->cd(13)->Clear();
  sumplots->cd(14)->Clear();
  sumplots->cd(15)->Clear();
  sumplots->cd(16)->Clear();
  sumplots->cd(17)->Clear();
  sumplots->cd(18)->Clear();
  sumplots->cd(19)->Clear();
  sumplots->cd(20)->Clear();
  sumplots->cd(21)->Clear();
  sumplots->cd(22)->Clear();
  sumplots->cd(23)->Clear();
  sumplots->cd(24)->Clear();

  sumtitlepad->cd();
  sumtitlepad->Clear();
  sumtitle->DrawLatex(0.4,0.4,"Summary of Events"); 
  summary->Print("CMDistortionAnalysisPhiRNeg.pdf)","pdf");

  return 0;
}


void WriteIntFluctFile(int ifile, int ihist, int nphi, double minphi, double maxphi, int nr,double minr,double maxr, int nz, double minzPos,double maxzPos, double minzNeg,double maxzNeg, TH3F *hFluctCharge){
    //2 hist: int of everything w small r n equal or larger z on same phi ; int of everthing w larger r ...
    //compare linear model to integrated fluctuation charge
    
    TH3F *hIntFluctChargeSmallRPos =new TH3F("hIntFluctChargeSmallRPos", "Integrated Fluctuation Charge, Positive Side, R Inside; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos); 
    TH3F *hIntFluctChargeSmallRNeg =new TH3F("hIntFluctChargeSmallRNeg", "Integrated Fluctuation Charge, Negative Side, R Inside; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzNeg,maxzNeg);

    TH3F *hIntFluctChargeLargeRPos =new TH3F("hIntFluctChargeLargeRPos", "Integrated Fluctuation Charge, Positive Side, R Outside; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos); 
    TH3F *hIntFluctChargeLargeRNeg =new TH3F("hIntFluctChargeLargeRNeg", "Integrated Fluctuation Charge, Negative Side, R Outside; phi (rad); r (cm); z (cm)", nphi,minphi,maxphi, nr,minr,maxr, nz,minzNeg,maxzNeg);

    
    int minbinR = hFluctCharge->GetYaxis()->FindBin(minr/100.); 
    int maxbinZPos = hFluctCharge->GetZaxis()->FindBin(maxzPos/100.);
    int minbinZNeg = hFluctCharge->GetZaxis()->FindBin(minzNeg/100.); 

    int maxbinR = hFluctCharge->GetYaxis()->FindBin(maxr/100.); 

    int minbinZPos = hFluctCharge->GetZaxis()->FindBin(minzPos/100.);
    int maxbinZNeg = hFluctCharge->GetZaxis()->FindBin(maxzNeg/100.); 
    
    /*  int nphi = 82;
    int nr = 54;
    int nz = 82;
    
    double minphi = -0.078539819;
    double minzPos = -1.3187500;
    double minzNeg = -106.81875;
    
    double maxphi = 6.3617253;
    double maxr = 79.115387;
    double maxzPos = 106.81875;*/
    
    for(int i = 1; i < nphi - 1; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin
      for(int j = 1; j < nr - 1; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); 
	for(int k = 1; k < nz - 1; k++){
	  double zPos = minzPos + ((maxzPos - minzPos)/(1.0*nz))*(k+0.5); 
	  double zNeg = minzNeg + ((maxzNeg - minzNeg)/(1.0*nz))*(k+0.5); 

	  double intfluctchargePos, intfluctchargeNeg;

	  int binPhi = hFluctCharge->GetXaxis()->FindBin(phi); 
	  int binR = hFluctCharge->GetYaxis()->FindBin(r/100.); 
	  int binZPos = hFluctCharge->GetZaxis()->FindBin(zPos/100.);
	  int binZNeg = hFluctCharge->GetZaxis()->FindBin(zNeg/100.); 
	  
	  
	  hIntFluctChargeSmallRPos->Fill(phi,r,zPos,hFluctCharge->Integral(binPhi, binPhi, minbinR, binR, binZPos, maxbinZPos));
	  
	    //intfluctchargePos =  hIntFluctChargeSmallRPos->Interpolate(phi,r,zPos);
	    //hCompareRTruevIntFluctSmallRPos->Fill(intfluctchargePos,);

	    /* hIntFluctChargeSmallRNeg->Fill(phi,r,zNeg,hFluctCharge->Integral());
	    intfluctchargeNeg =  hIntFluctChargeSmallRNeg->Interpolate(phi,r,zNeg);
	    hCompareRTruevIntFluctSmallRNeg->Fill(intfluctchargeNeg,);*/



	  hIntFluctChargeLargeRPos->Fill(phi,r,zPos,hFluctCharge->Integral(binPhi, binPhi, binR, maxbinR, binZPos, maxbinZPos));
	    //intfluctchargePos =  hIntFluctChargeLargeRPos->Interpolate(phi,r,zPos);
	    //hCompareRTruevIntFluctLargeRPos->Fill(intfluctchargePos,);
	    
	    /*  hIntFluctChargeLargeRNeg->Fill(phi,r,zNeg,hFluctCharge->Integral());
	    intfluctchargeNeg =  hIntFluctChargeLargeRNeg->Interpolate(phi,r,zNeg);
	    hCompareRTruevIntFluctLargeRNeg->Fill(intfluctchargeNeg,);*/


	}
      }
      cout << "finished with phi: " << phi << endl;
    }

    TFile *integ;

    integ=TFile::Open(Form("IntFluctEvent%d.root", (2*ifile+ihist)), "RECREATE");

    hIntFluctChargeSmallRPos->Write();
    hIntFluctChargeLargeRPos->Write();

    integ->Close();

}
