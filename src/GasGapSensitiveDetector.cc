#include "GasGapSensitiveDetector.hh"
#include "GPMAnalysis.hh"
#include "GasGapHit.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4HCofThisEvent.hh"
#include "G4HCtable.hh"
#include "G4VProcess.hh"

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"

GasGapSensitiveDetector::GasGapSensitiveDetector(G4String SDname) :
   G4VSensitiveDetector(SDname),
   charge(0),
   primaryene(0.), 
   zinteraction(0), 
   contaPrimary(0),
   contaInteraction(0),
   contaSec(0),
   contaSec_B(0),
   contaTrack(0),
   contaGar(0)
{
   G4cout << "**************************************" << G4endl;
   G4cout << "** Creating SD With name " << SDname << "**" << G4endl;
   G4cout << "**************************************" << G4endl;

   G4String myCollectionName = "GasGapHitCollection";
   collectionName.insert(myCollectionName);

   ttTrack.clear();
   ttTrack_B.clear();
   ttTrack_Gar.clear();
   postTrack.clear();
   container.clear();
   zinteraction = -5;
}

GasGapSensitiveDetector::~GasGapSensitiveDetector() {
}

G4bool GasGapSensitiveDetector::ProcessHits(G4Step *step, G4TouchableHistory *) {
   G4double edep = step->GetTotalEnergyDeposit();
   //if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
   //if (edep == 0.) return false;

   G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
   G4StepPoint* thePrePoint = step->GetPreStepPoint();
   G4VPhysicalVolume* thePrePV = touchable->GetVolume();
   G4double pre_energy = thePrePoint->GetKineticEnergy();
   
   G4TouchableHandle postTouchable = step->GetPostStepPoint()->GetTouchableHandle();
   G4StepPoint* thePostPoint = step->GetPostStepPoint();
   G4VPhysicalVolume* thePostPV = postTouchable->GetVolume();
   G4double post_energy = thePostPoint->GetKineticEnergy();

   if (thePrePV == thePostPV || thePostPV->GetName() != "GasGap1" || thePrePV->GetName() != "PC") return false;

   //if (post_energy != 0 ) return false; 

   G4int copyNo = touchable->GetVolume(0)->GetCopyNo();
   G4int layerIndex = copyNo;
   G4String volName = touchable->GetVolume(0)->GetName();

   G4Track* track = step->GetTrack();
   G4int trackIndex = track->GetTrackID();

  // Eit = std::find(EindexId.begin(), EindexId.end(), trackIndex);
  // if (Eit == EindexId.end()) {
  //   EindexId.push_back(trackIndex);
  //   Trackhold += 1;
  // } else {
  //   return false; 
  // }
   //G4cout<< "Number of Tracks in Scint :" << Trackhold << G4endl;
   G4double energy = step->GetPreStepPoint()->GetKineticEnergy();

   double t = track->GetGlobalTime();
   double x = track->GetPosition().getX();
   double y = track->GetPosition().getY();
   double z = track->GetPosition().getZ();

   double px = track->GetMomentum().getX();
   double py = track->GetMomentum().getY();
   double pz = track->GetMomentum().getZ();

   G4int pdg = track->GetParticleDefinition()->GetPDGEncoding();
   G4int charge = track->GetParticleDefinition()->GetPDGCharge();
   G4StepPoint* point = step->GetPostStepPoint();
   const G4VProcess* proc = point->GetProcessDefinedStep();
   const G4String procname = proc->GetProcessName();

   G4double edepI = edep-step->GetNonIonizingEnergyDeposit();
   G4int parentIndex = track->GetParentID();
   double genz = track->GetVertexPosition().getZ();

  // G4cout << "================================================== " << G4endl;
  // G4cout << " Particle            " << track->GetDefinition() << G4endl;
  // G4cout << " Track Number        " << trackIndex << G4endl;
  // G4cout << " Parent Track Number " << parentIndex << G4endl;
  // G4cout << " Pre PV              " << thePrePV->GetName() << G4endl;
  // G4cout << " Post PV             " << thePostPV->GetName() << G4endl;
  // G4cout << " post Energy         " << post_energy << G4endl;
  // G4cout << " pre Energy          " << energy << G4endl;
  // G4cout << "================================================== " << G4endl;

   G4String genprocess;
   if (track->GetCreatorProcess()!=0) {
      const G4VProcess* genproc = track->GetCreatorProcess();
      genprocess = genproc->GetProcessName();
   } else {
      genprocess = "primary";
   }

   const G4LogicalVolume * genLogVolume = track->GetLogicalVolumeAtVertex();
   G4String genvolume = genLogVolume->GetName();

   if ((*step->GetSecondary()).size()>0  && trackIndex == 1 && contaInteraction == 0) {
      zinteraction = z;
    
      contaInteraction = 1; 
   }

   container[trackIndex] = parentIndex; 

   G4int generation = GetGeneration(trackIndex);

   if (trackIndex == 1 && (volName == "FakeBottom" || volName == "FakeTop") && contaPrimary == 0) {
      primaryene = energy;
      contaPrimary = 1;
   }
   
   GPMAnalysis::GetInstance()->SetEnergyDeposition(volName, edep, edepI, t);

 //  if (volName == "GasGap1" || volName == "GasGap2") {
     if (volName == "PC") {
      for (G4int T = 0; T < ttTrack.size(); T++) {
         if (ttTrack[T] == trackIndex) { 
            contaSec = 9999;
            break;
         }
      }

      if (contaSec != 9999) {
         GPMAnalysis::GetInstance()->SaveGapTrack(pdg, charge, generation, genprocess, genvolume, genz, volName, energy); 
         contaSec = trackIndex;
         ttTrack.push_back(trackIndex);
      }

     // if (volName == "GasGap1") {
      if (volName == "PC") {
         for (G4int T=0; T < ttTrack_Gar.size(); T++) {
            if (ttTrack_Gar[T] == trackIndex) { 
               contaSec = 9999;
               break;
            }
         }
         if (contaSec != 9999) {
            GPMAnalysis::GetInstance()->SaveGarfieldQuantities(pdg, energy, x, y, z, px, py, pz,1);
            contaGar = trackIndex;
            ttTrack_Gar.push_back(trackIndex);
         }
      }
   }

   G4String isPri = step->GetTrack()->GetTrackID() == 1 ? "Yes" : "No";

   hitMap_t::iterator it = hitMap.find(layerIndex);
   GasGapHit* aHit = 0;
   if (it != hitMap.end()) {
      aHit = it->second;
   } else {
      aHit = new GasGapHit(layerIndex);
      hitMap.insert(std::make_pair(layerIndex, aHit));
      hitCollection->insert(aHit);
   }
   aHit->AddEdep(edep);
   return true;
}

void GasGapSensitiveDetector::Initialize(G4HCofThisEvent* HCE) {
   hitCollection = new GasGapHitCollection(GetName(), collectionName[0]);
   static G4int HCID = -1;
   if (HCID<0) HCID = GetCollectionID(0); 
   HCE->AddHitsCollection(HCID, hitCollection);

   hitMap.clear();
}

void GasGapSensitiveDetector::EndOfEvent(G4HCofThisEvent*) {
   G4int factor = 5;

   GPMAnalysis::GetInstance()->SavePrimary(primaryene, zinteraction);

   ttTrack.clear();
   ttTrack_B.clear();
   ttTrack_Gar.clear();
   postTrack.clear();
   container.clear();

   contaPrimary = 0;
   contaInteraction = 0;
   contaSec = 0;
   contaSec_B = 0;
   contaTrack = 0;
   contaGar = 0;
   primaryene = 0.;
   zinteraction =-5.;
}

G4int GasGapSensitiveDetector::GetGeneration(G4int index) {
   G4int i = 1;
   while (container[index] != 0) {
      index = container[index];
      i++;
   }
   return i;
}
