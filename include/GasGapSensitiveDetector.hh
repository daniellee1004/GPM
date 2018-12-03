#ifndef GasGapSensitiveDetector_h
#define GasGapSensitiveDetector_h 1

#include "GasGapHit.hh"
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4TouchableHistory;

class GasGapSensitiveDetector : public G4VSensitiveDetector {
   public:
      GasGapSensitiveDetector(G4String SDname);
      ~GasGapSensitiveDetector();

   public:
      G4bool ProcessHits(G4Step *step, G4TouchableHistory * ROhist);

      void Initialize(G4HCofThisEvent* HCE);
      void EndOfEvent(G4HCofThisEvent* HCE);

      G4int GetGeneration(G4int index);

   private:
      G4int charge;
      G4double primaryene;
      G4double zinteraction;
      G4int contaPrimary;
      G4int contaInteraction;
      G4int contaSec;
      G4int contaSec_B;
      G4int contaTrack;
      G4int contaGar;
      
      G4int Trackhold;
      std::vector<G4int> EindexId;
      std::vector<G4int>::iterator Eit;

      std::vector<G4int> ttTrack;
      std::vector<G4int> ttTrack_B;
      std::vector<G4int> ttTrack_Gar;
      std::vector<G4int> postTrack;

      typedef std::map<G4int, GasGapHit*> hitMap_t;
      hitMap_t hitMap;
      GasGapHitCollection* hitCollection;
      std::map<G4int, G4int> container;
};

#endif

