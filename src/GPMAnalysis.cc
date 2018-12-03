#include "GPMAnalysis.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4String.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <math.h>

#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>

GPMAnalysis* GPMAnalysis::singleton = 0;

GPMAnalysis* GPMAnalysis::GetInstance() 
{
   if (singleton == 0) {
      static GPMAnalysis analysis;
      singleton = &analysis;
   }
   return singleton;
}

GPMAnalysis::~GPMAnalysis() 
{}

GPMAnalysis::GPMAnalysis() 
{
   m_ROOT_file = 0;
}

void GPMAnalysis::PrepareNewRun(const G4Run*) 
{   
   eventCounter = 0;
   n_gamma = 0;
   n_electron = 0;
   n_positron = 0;
   G4cout << "Enter GPM Analysis" << G4endl;

  // m_ROOT_file = new TFile("./Results/WL_5keV.root", "RECREATE");
   m_ROOT_file = new TFile("./50keV_GPM.root", "RECREATE");
   if (m_ROOT_file) {
      G4cout << "ROOT file (WL_50keV.root) made" << G4endl;
   } else {
      G4cout << "ROOT file (WL_50keV.root) not made" << G4endl;
   }

   //Create Histograms 
   m_ROOT_histo1 = new TH1D("Sil_Edep","Energy Deposition in Scint",100,0,20);
   //Create tree and Branch 
   t = new TTree("Events", "Simulated Event Data");

   t->Branch("nEvent", &nEvent, "nEvnet/I");
   t->Branch("nScint_Op", &Op_Scint, "Op_Scint/I");
   t->Branch("nCer_Op", &Op_Cer, "Op_Cer/I");
   t->Branch("ScintEdep", &scintEdep, "sintEdep/D");
   t->Branch("nGammaScint", &nGammaScint, "nGammaScint/I");
   t->Branch("nElecGasGap1", &nElec_GasGap1, "nElec_GasGap1/I");
   t->Branch("edep", &edep);
   t->Branch("edepI", &edepI);
   t->Branch("edepTime", &edepTime);
   t->Branch("primaryEne", &primaryEne, "primaryEne/D");
   t->Branch("zInteraction", &zInteraction, "zInteraction/D");
   
   t->Branch("EleGap", &eleGap, "EleGap/I");
   t->Branch("PosGap", &posGap, "PosGap/I");
   t->Branch("ChargeGap", &chargeGap, "ChargeGap/I");

   t->Branch("gapTrackPart", &gapTrackPart);
   t->Branch("gapTrackCharge", &gapTrackCharge);
   t->Branch("gapTrackGeneration", &gapTrackGeneration);
   t->Branch("gapTrackEne", &gapTrackEne);
   // stuff
   t->Branch("nOpPhoton",&nOpPhoton, "nOpPhoton/I");
   t->Branch("pdgCode",&pdgCode);
   t->Branch("kineticEnergy",&kineticEnergy);
   t->Branch("positionX",&positionX);
   t->Branch("positionY",&positionY);
   t->Branch("positionZ",&positionZ);
   t->Branch("momentumX",&momentumX);
   t->Branch("momentumY",&momentumY);
   t->Branch("momentumZ",&momentumZ);
   t->Branch("Wavelength",&wavelength);
   t->Branch("dR", &dR);
   // stuff
  
}

void GPMAnalysis::PrepareNewEvent(const G4Event*) 
{
   isNewEvent = true;
   //Clearing stuff branch
   nEvent = 1;
   scintEdep = 0;
   nElec_GasGap1 = 0;
   nGammaScint = 0;
   
   Op_Cer = 0;
   Op_Scint = 0;
   edep.clear();
   edepI.clear();
   edepTime.clear();

   newTrack = true;
   EindexId.clear();
   GindexId.clear();
   
   eleGap = 0;
   posGap = 0;
   chargeGap = 0;
   
   nOpPhoton = 0;
   pdgCode.clear();
   kineticEnergy.clear();
   dR.clear();
   positionX.clear();
   positionY.clear();
   positionZ.clear();
   momentumX.clear();
   momentumY.clear();
   momentumZ.clear();
   wavelength.clear();
   
   gapTrackPart.clear();
   gapTrackCharge.clear();
   gapTrackGeneration.clear();
   gapTrackGenZ.clear();
   gapTrackEne.clear();
}

void GPMAnalysis::EndOfEvent(const G4Event*) 
{
   //fill histograms and branches 
   
   t->Fill();
}

void GPMAnalysis::EndOfRun(const G4Run* aRun) 
{
   // Write / close root files
   G4int numEvents = aRun->GetNumberOfEvent();

   m_ROOT_file->cd();
   t->Write();
   G4cout << "Writing ROOT files ..." << G4endl;
   m_ROOT_file->Write();
   G4cout << "Closing ROOT files ..." << G4endl;
   m_ROOT_file->Close();
   delete m_ROOT_file;
   
}
void GPMAnalysis::AddSecondary(const G4ParticleDefinition* part) 
{
   if (part == G4Gamma::Gamma())            { ++n_gamma; }
   else if (part == G4Electron::Electron()) { ++n_electron; }
   else if (part == G4Positron::Positron()) { ++n_positron; }
}

void GPMAnalysis::AddGapSecondary(const G4ParticleDefinition* part, G4int gapNum) 
{
   gapNum-- ; 
   if (part == G4Gamma::Gamma())            { ++n_gapGamma[gapNum]; }
   else if (part == G4Electron::Electron()) { ++n_gapElectron[gapNum]; }
   else if (part == G4Positron::Positron()) { ++n_gapPositron[gapNum]; }
}

void GPMAnalysis::AddScintEDep(G4double Sedep) // In SteppingAction
{
  scintEdep += Sedep;
  m_ROOT_histo1->Fill(Sedep);
}

void GPMAnalysis::SetEnergyDeposition(std::string someVolume, G4double someEdep, G4double someEdepI, G4double someTime) // in GasGapSensitiveDetector
{
   edep.push_back(someEdep);
   edepI.push_back(someEdepI);
   edepTime.push_back(someTime);
}

void GPMAnalysis::ScintPCProd(G4int PDG, G4String V_volumeName, G4String C_volumeName, G4int trackIndex) // In SteppingAction
{
   if ((PDG == 11 || PDG == -11) && (C_volumeName == "GasGap1") && (V_volumeName == "PhotoCathodeLog")) {
      Eit = std::find(EindexId.begin(), EindexId.end(), trackIndex);
      if (Eit == EindexId.end()) {
         EindexId.push_back(trackIndex);
         nElec_GasGap1 += 1;
      }
   } else if ((PDG == 22) && (C_volumeName == "PhotoCathode") && (V_volumeName == "ScintillatorLog")) {
      Git = std::find(GindexId.begin(), GindexId.end(), trackIndex);
      if (Git == GindexId.end()) {
         GindexId.push_back(trackIndex);
         nGammaScint += 1;
      }
   }
}

void GPMAnalysis::OpPhoton(G4String p_Name) //in StackingAction
{
   if (p_Name == "Scintillation") Op_Scint++;
   if (p_Name == "Cerenkov") Op_Cer++;
}


void GPMAnalysis::SaveGarfieldQuantities(
      G4int aPdgCode,
      G4double aKineticEnergy,
      G4double aPositionX,
      G4double aPositionY,
      G4double aPositionZ,
      G4double aMomentumX,
      G4double aMomentumY,
      G4double aMomentumZ,
      G4int opPhoton) 
{
   aPdgCode = 11;
   nOpPhoton += opPhoton;
   pdgCode.push_back(aPdgCode);
   kineticEnergy.push_back(aKineticEnergy);
   positionX.push_back(aPositionX);
   positionY.push_back(aPositionY);
   positionZ.push_back(aPositionZ);
   momentumX.push_back(aMomentumX);
   momentumY.push_back(aMomentumY);
   momentumZ.push_back(aMomentumZ);
   
   dR.push_back(sqrt(pow(aPositionX,2) + pow(aPositionY,2)));
   double Joules = 1.60218*pow(10,-19);
   double Plank = 6.62607*pow(10,-34);
   double Clight = 3*pow(10,8);

   double WL = 1/((aKineticEnergy*Joules)/(Plank*Clight));
   wavelength.push_back(WL);   
}

void GPMAnalysis::SavePrimary(G4double primaryene, G4double zinteraction) 
{
   primaryEne = primaryene;
   zInteraction = zinteraction; 
}

void GPMAnalysis::SaveGapTrack(             // in GasGasSensitiveDetector
      G4int gapPart,
      G4int gapCharge,
      G4int generation,
      std::string genprocess,
      std::string genvolume,
      G4double genz,
      std::string volname,
      G4double kinene) 
{
   if (genprocess == "primary") return;
   if (gapCharge != 0) chargeGap = 1;
   if (gapPart == 11) eleGap = 1;
   if (gapPart == -11) posGap = 1;

   gapTrackPart.push_back(gapPart);
   gapTrackCharge.push_back(gapCharge);
   gapTrackGeneration.push_back(generation);
   gapTrackGenZ.push_back(genz);
   gapTrackEne.push_back(kinene);
}
/*
void GPMAnalysis::WorldPCFrac(G4int IntoWorld, G4int IntoPC) 
{
   WoldFrac = IntoWorld;
   PCFrac = IntoPC;
}
*/
