// step 2 with phi,r coords
#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"
#include <TTime.h>

using namespace std;

int CMDistortionRecoPhiRFull(int nMaxEvents = -1) {
  int nbins = 35; 
  double low = -80.0;
  double high = 80.0;
  double deltaXPos, deltaYPos, deltaZPos, deltaRPos, deltaPhiPos;
  double deltaXNeg, deltaYNeg, deltaZNeg, deltaRNeg, deltaPhiNeg;
  int nEvents;
    
  //take in events
  const char * inputpattern="/sphenix/u/skurdi/CMCalibration/cmDistHitsFullTree_Event*.root"; 
  
  //find all files that match the input string (includes wildcards)
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  TString sourcefilename;
  
  //how many events
  if (nMaxEvents<0){
    nEvents=filelist->GetNFiles();
  } else if(nMaxEvents<filelist->GetNFiles()){
    nEvents=nMaxEvents;
  } else {
    nEvents= filelist->GetNFiles();
  }
  
  TCanvas *canvas1=new TCanvas("canvas1","CMDistortionReco1",1200,800);
  canvas1->Divide(3,2);

  //canvas for time plot
  TCanvas *canvas=new TCanvas("canvas","CMDistortionReco2",400,400);
  
  TVector3 *position, *newpositionPos, *newpositionNeg;
  position = new TVector3(1.,1.,1.);
  newpositionPos = new TVector3(1.,1.,1.);
  newpositionNeg = new TVector3(1.,1.,1.);

  //histogram to compare times
  TH1F *hTimePerEvent = new TH1F("hTimePerEvent","Time Per Event; time (ms)",20,0,10000);
    
  for (int ifile=0;ifile < nEvents;ifile++){
    //call to TTime before opening ttree
    TTime now;
    now=gSystem->Now();
    unsigned long before = now;
    
    //get data from ttree
    sourcefilename=((TFileInfo*)(filelist->GetList()->At(ifile)))->GetCurrentUrl()->GetFile();
    
    char const *treename="cmDistHitsTree";
    TFile *input=TFile::Open(sourcefilename, "READ");
    TTree *inTree=(TTree*)input->Get("tree");
    
    inTree->SetBranchAddress("position",&position);
    inTree->SetBranchAddress("newpositionPos",&newpositionPos);
    inTree->SetBranchAddress("newpositionNeg",&newpositionNeg);

    //hardcoded numbers from average distortion file's hIntDistortionPosX
    int nbinsphi = 30; //when using 35, blank spots at around r = 22 cm, phi just above n below pi
    double lowphi = -0.078539819;
    double highphi = 6.3617253;
    int nbinsr = 30; // when using 35, blank stripe around r = 58 cm
    double lowr = 0.0;
    double highr = 90.0;
    
    //for forward only

    TH2F *hStripesPerBinPhiR = new TH2F("hStripesPerBinPhiR","CM Stripes Per Bin (z in stripes); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);

    //positive
    TH2F *hCartesianForwardPhiRPos[3];
    hCartesianForwardPhiRPos[0] = new TH2F("hForwardX_PhiR_Pos","X Shift Forward of Stripe Centers, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    hCartesianForwardPhiRPos[1] = new TH2F("hForwardY_PhiR_Pos","Y Shift Forward of Stripe Centers, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    hCartesianForwardPhiRPos[2] = new TH2F("hForwardZ_PhiR_Pos","Z Shift Forward of Stripe Centers, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    
     TH2F *hCylindricalForwardPhiRPos[2];
    hCylindricalForwardPhiRPos[0] = new TH2F("hForwardR_PhiR_Pos","Radial Shift Forward of Stripe Centers, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    hCylindricalForwardPhiRPos[1] = new TH2F("hForwardPhi_PhiR_Pos","Phi Shift Forward of Stripe Centers, Phi,R binning, Positive Side (rad); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);

    //negative
    TH2F *hCartesianForwardPhiRNeg[3];
    hCartesianForwardPhiRNeg[0] = new TH2F("hForwardX_PhiR_Neg","X Shift Forward of Stripe Centers, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    hCartesianForwardPhiRNeg[1] = new TH2F("hForwardY_PhiR_Neg","Y Shift Forward of Stripe Centers, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    hCartesianForwardPhiRNeg[2] = new TH2F("hForwardZ_PhiR_Neg","Z Shift Forward of Stripe Centers, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    
     TH2F *hCylindricalForwardPhiRNeg[2];
    hCylindricalForwardPhiRNeg[0] = new TH2F("hForwardR_PhiR_Neg","Radial Shift Forward of Stripe Centers, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    hCylindricalForwardPhiRNeg[1] = new TH2F("hForwardPhi_PhiR_Neg","Phi Shift Forward of Stripe Centers, Phi,R binning, Negative Side (rad); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    
    for (int i=0;i<inTree->GetEntries();i++){
      inTree->GetEntry(i);

      double r = position->Perp();
      double phi = position->Phi();

      if(position->Phi() < 0.0){
	phi = position->Phi() + 2.0*TMath::Pi(); 
      }
      
      hStripesPerBinPhiR->Fill(phi,r,1);

      //positive
      deltaXPos = (newpositionPos->X() - position->X())*(1e4); //convert from cm to micron 
      deltaYPos = (newpositionPos->Y() - position->Y())*(1e4);
      deltaZPos = (newpositionPos->Z() - position->Z())*(1e4);

      deltaRPos = (newpositionPos->Perp() - position->Perp())*(1e4);
      deltaPhiPos = newpositionPos->DeltaPhi(*position);

      hCartesianForwardPhiRPos[0]->Fill(phi,r,deltaXPos);
      hCartesianForwardPhiRPos[1]->Fill(phi,r,deltaYPos);
      hCartesianForwardPhiRPos[2]->Fill(phi,r,deltaZPos);

      hCylindricalForwardPhiRPos[0]->Fill(phi,r,deltaRPos);
      hCylindricalForwardPhiRPos[1]->Fill(phi,r,deltaPhiPos);

      //negative
      deltaXNeg = (newpositionNeg->X() - position->X())*(1e4); //convert from cm to micron 
      deltaYNeg = (newpositionNeg->Y() - position->Y())*(1e4);
      deltaZNeg = (newpositionNeg->Z() - position->Z())*(1e4);

      deltaRNeg = (newpositionNeg->Perp() - position->Perp())*(1e4);
      deltaPhiNeg = newpositionNeg->DeltaPhi(*position);

      hCartesianForwardPhiRNeg[0]->Fill(phi,r,deltaXNeg);
      hCartesianForwardPhiRNeg[1]->Fill(phi,r,deltaYNeg);
      hCartesianForwardPhiRNeg[2]->Fill(phi,r,deltaZNeg);

      hCylindricalForwardPhiRNeg[0]->Fill(phi,r,deltaRNeg);
      hCylindricalForwardPhiRNeg[1]->Fill(phi,r,deltaPhiNeg);
    }

    //positive
    TH2F *hCartesianAveShiftPhiRPos[3];
    hCartesianAveShiftPhiRPos[0] = new TH2F("AveShiftX_PhiR_Pos","Average of CM Model X over Stripes per Bin, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 
    hCartesianAveShiftPhiRPos[1] = new TH2F("AveShiftY_PhiR_Pos","Average of CM Model Y over Stripes per Bin, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 
    hCartesianAveShiftPhiRPos[2] = new TH2F("AveShiftZ_PhiR_Pos","Average of CM Model Z over Stripes per Bin, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 

    TH2F *hCylindricalAveShiftPhiRPos[2];
     hCylindricalAveShiftPhiRPos[0] = new TH2F("AveShiftR_PhiR_Pos","Average of CM Model R over Stripes per Bin, Phi,R binning, Positive Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 
    hCylindricalAveShiftPhiRPos[1] = new TH2F("AveShiftPhi_PhiR_Pos","Average of CM Model Phi over Stripes per Bin, Phi,R binning, Positive Side (rad); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);

    //negative
    TH2F *hCartesianAveShiftPhiRNeg[3];
    hCartesianAveShiftPhiRNeg[0] = new TH2F("AveShiftX_PhiR_Neg","Average of CM Model X over Stripes per Bin, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 
    hCartesianAveShiftPhiRNeg[1] = new TH2F("AveShiftY_PhiR_Neg","Average of CM Model Y over Stripes per Bin, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 
    hCartesianAveShiftPhiRNeg[2] = new TH2F("AveShiftZ_PhiR_Neg","Average of CM Model Z over Stripes per Bin, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 

    TH2F *hCylindricalAveShiftPhiRNeg[2];
     hCylindricalAveShiftPhiRNeg[0] = new TH2F("AveShiftR_PhiR_Neg","Average of CM Model R over Stripes per Bin, Phi,R binning, Negative Side (#mum); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr); 
    hCylindricalAveShiftPhiRNeg[1] = new TH2F("AveShiftPhi_PhiR_Neg","Average of CM Model Phi over Stripes per Bin, Phi,R binning, Negative Side (rad); phi (rad); r (cm)",nbinsphi,lowphi,highphi,nbinsr,lowr,highr);
    
    for (int i = 0; i < 3; i ++){
      hCartesianAveShiftPhiRPos[i]->Divide(hCartesianForwardPhiRPos[i],hStripesPerBinPhiR);
      hCartesianAveShiftPhiRNeg[i]->Divide(hCartesianForwardPhiRNeg[i],hStripesPerBinPhiR);
    }
    
    //pos
    hCylindricalAveShiftPhiRPos[0]->Divide(hCylindricalForwardPhiRPos[0],hStripesPerBinPhiR);
    hCylindricalAveShiftPhiRPos[1]->Divide(hCylindricalForwardPhiRPos[1],hStripesPerBinPhiR);

    //neg
    hCylindricalAveShiftPhiRNeg[0]->Divide(hCylindricalForwardPhiRNeg[0],hStripesPerBinPhiR);
    hCylindricalAveShiftPhiRNeg[1]->Divide(hCylindricalForwardPhiRNeg[1],hStripesPerBinPhiR);
  
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

    //positive
    TH3F *hCartesianCMModelPhiRPos[3];
    hCartesianCMModelPhiRPos[0]=new TH3F("hCMModelX_PhiR_Pos", "CM Model: X Shift Forward of Stripe Centers, Phi,R binning, Positive Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos); //rad, cm, cm
    hCartesianCMModelPhiRPos[1]=new TH3F("hCMModelY_PhiR_Pos", "CM Model: Y Shift Forward of Stripe Centers, Phi,R binning, Positive Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
    hCartesianCMModelPhiRPos[2]=new TH3F("hCMModelZ_PhiR_Pos", "CM Model: Z Shift Forward of Stripe Centers, Phi,R binning, Positive Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);

    TH3F *hCylindricalCMModelPhiRPos[2];
    hCylindricalCMModelPhiRPos[0]=new TH3F("hCMModelR_PhiR_Pos", "CM Model: Radial Shift Forward of Stripe Centers, Phi,R binning, Positive Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);
    hCylindricalCMModelPhiRPos[1]=new TH3F("hCMModelPhi_PhiR_Pos", "CM Model: Phi Shift Forward of Stripe Centers, Phi,R binning, Positive Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzPos,maxzPos);

    //negative
    TH3F *hCartesianCMModelPhiRNeg[3];
    hCartesianCMModelPhiRNeg[0]=new TH3F("hCMModelX_PhiR_Neg", "CM Model: X Shift Forward of Stripe Centers, Phi,R binning, Negative Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzNeg,maxzNeg); //rad, cm, cm
    hCartesianCMModelPhiRNeg[1]=new TH3F("hCMModelY_PhiR_Neg", "CM Model: Y Shift Forward of Stripe Centers, Phi,R binning, Negative Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzNeg,maxzNeg);
    hCartesianCMModelPhiRNeg[2]=new TH3F("hCMModelZ_PhiR_Neg", "CM Model: Z Shift Forward of Stripe Centers, Phi,R binning, Negative Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzNeg,maxzNeg);

    TH3F *hCylindricalCMModelPhiRNeg[2];
    hCylindricalCMModelPhiRNeg[0]=new TH3F("hCMModelR_PhiR_Neg", "CM Model: Radial Shift Forward of Stripe Centers, Phi,R binning, Negative Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzNeg,maxzNeg);
    hCylindricalCMModelPhiRNeg[1]=new TH3F("hCMModelPhi_PhiR_Neg", "CM Model: Phi Shift Forward of Stripe Centers, Phi,R binning, Negative Side", nphi,minphi,maxphi, nr,minr,maxr, nz,minzNeg,maxzNeg); 
      
    double xshiftPhiRPos, yshiftPhiRPos, zshiftPhiRPos, rshiftPhiRPos, phishiftPhiRPos;
    double xshiftPhiRNeg, yshiftPhiRNeg, zshiftPhiRNeg, rshiftPhiRNeg, phishiftPhiRNeg;
  
    for(int i = 0; i < nphi; i++){
      double phi = minphi + ((maxphi - minphi)/(1.0*nphi))*(i+0.5); //center of bin

      for(int j = 0; j < nr; j++){
	double r = minr + ((maxr - minr)/(1.0*nr))*(j+0.5); //center of bin
	
	for(int k = 0; k < nz; k++){
	  double zPos = minzPos + ((maxzPos - minzPos)/(1.0*nz))*(k+0.5); //center of bin
	  double zNeg = minzNeg + ((maxzNeg - minzNeg)/(1.0*nz))*(k+0.5); //center of bin

	  //positive
	  xshiftPhiRPos=(hCartesianAveShiftPhiRPos[0]->Interpolate(phi,r))*(1e-4);//coordinate of your stripe
	  yshiftPhiRPos=(hCartesianAveShiftPhiRPos[1]->Interpolate(phi,r))*(1e-4);//convert micron to cm
	  zshiftPhiRPos=(hCartesianAveShiftPhiRPos[2]->Interpolate(phi,r))*(1e-4);

	  rshiftPhiRPos=(hCylindricalAveShiftPhiRPos[0]->Interpolate(phi,r))*(1e-4);
	  phishiftPhiRPos=hCylindricalAveShiftPhiRPos[1]->Interpolate(phi,r);
	  
	  hCartesianCMModelPhiRPos[0]->Fill(phi,r,zPos,xshiftPhiRPos*(1-zPos/105.5));
	  hCartesianCMModelPhiRPos[1]->Fill(phi,r,zPos,yshiftPhiRPos*(1-zPos/105.5));
	  hCartesianCMModelPhiRPos[2]->Fill(phi,r,zPos,zshiftPhiRPos*(1-zPos/105.5));

	  hCylindricalCMModelPhiRPos[0]->Fill(phi,r,zPos,rshiftPhiRPos*(1-zPos/105.5));
	  hCylindricalCMModelPhiRPos[1]->Fill(phi,r,zPos,phishiftPhiRPos*(1-zPos/105.5));

	  //negative
	  xshiftPhiRNeg=(hCartesianAveShiftPhiRNeg[0]->Interpolate(phi,r))*(1e-4);//coordinate of your stripe
	  yshiftPhiRNeg=(hCartesianAveShiftPhiRNeg[1]->Interpolate(phi,r))*(1e-4);//convert micron to cm
	  zshiftPhiRNeg=(hCartesianAveShiftPhiRNeg[2]->Interpolate(phi,r))*(1e-4);

	  rshiftPhiRNeg=(hCylindricalAveShiftPhiRNeg[0]->Interpolate(phi,r))*(1e-4);
	  phishiftPhiRNeg=hCylindricalAveShiftPhiRNeg[1]->Interpolate(phi,r);
	  
	  hCartesianCMModelPhiRNeg[0]->Fill(phi,r,zNeg,xshiftPhiRNeg*(1+zNeg/105.5));
	  hCartesianCMModelPhiRNeg[1]->Fill(phi,r,zNeg,yshiftPhiRNeg*(1+zNeg/105.5));
	  hCartesianCMModelPhiRNeg[2]->Fill(phi,r,zNeg,zshiftPhiRNeg*(1+zNeg/105.5));

	  hCylindricalCMModelPhiRNeg[0]->Fill(phi,r,zNeg,rshiftPhiRNeg*(1+zNeg/105.5));
	  hCylindricalCMModelPhiRNeg[1]->Fill(phi,r,zNeg,phishiftPhiRNeg*(1+zNeg/105.5));
	}
      }
    }
    
    TFile *plots;

    plots=TFile::Open(Form("CMModelsPhiRFull_Event%d.root",ifile),"RECREATE");

    for(int i = 0; i < 3; i++){
      hCartesianForwardPhiRPos[i]->Write();
      hCartesianAveShiftPhiRPos[i]->Write();
      hCartesianCMModelPhiRPos[i]->Write();

      hCartesianForwardPhiRNeg[i]->Write();
      hCartesianAveShiftPhiRNeg[i]->Write();
      hCartesianCMModelPhiRNeg[i]->Write();
    }

    for(int i = 0; i < 2; i++){
      hCylindricalForwardPhiRPos[i]->Write();
      hCylindricalAveShiftPhiRPos[i]->Write();
      hCylindricalCMModelPhiRPos[i]->Write();

      hCylindricalForwardPhiRNeg[i]->Write();
      hCylindricalAveShiftPhiRNeg[i]->Write();
      hCylindricalCMModelPhiRNeg[i]->Write();
    }
    
    plots->Close();

    
    //call to TTime after outputting TH3Fs
    now=gSystem->Now();
    unsigned long after = now;
    
    hTimePerEvent->Fill(after-before);
    
    //to check histograms
    for (int i = 0; i < 3; i++){
      hCartesianForwardPhiR[i]->SetStats(0);
    }

    canvas1->cd(1);
    hCartesianForwardPhiRPos[0]->Draw("colz");
    canvas1->cd(2);
    hCartesianForwardPhiRPos[1]->Draw("colz");
    canvas1->cd(3);
    hCartesianForwardPhiRPos[2]->Draw("colz");
    canvas1->cd(4);
    hCartesianForwardPhiRNeg[0]->Draw("colz");
    canvas1->cd(5);
    hCartesianForwardPhiRNeg[1]->Draw("colz");
    canvas1->cd(6);
    hCartesianForwardPhiRNeg[2]->Draw("colz");
  
    if(ifile == 0){ 
      canvas1->Print("CMDistortionReco1.pdf(","pdf");
    } else if (ifile == nEvents - 1){
      canvas1->Print("CMDistortionReco1.pdf)","pdf");
    } else {
      canvas1->Print("CMDistortionReco1.pdf","pdf");
    }
    
  }

  canvas->cd();
  hTimePerEvent->Draw();
  canvas->Print("CMDistortionReco2.pdf","pdf");

  return 0;
}
