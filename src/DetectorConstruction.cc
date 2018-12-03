//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file /src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "GasGapSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"

#include "G4GeometryManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
  fExpHall_x = fExpHall_y = fExpHall_z = 1.0*m;
  //fTank_x    = fTank_y    = fTank_z    =  5.0*m;
  //fBubble_x  = fBubble_y  = fBubble_z  =  0.5*m;
  fGem_x  = fGem_y  =  0.1*m;
  fScint_z = 50.0*mm;
  fPC_z = 20.0*nm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

// ------------- Materials -------------

  G4double density(0.), fractionMass(0.);
  G4int nelements, natoms, numel;

  G4NistManager* manager = G4NistManager::Instance();
  G4Element* elH = manager->FindOrBuildElement(1);
  G4Element* elC = manager->FindOrBuildElement(6);
  G4Element* elTl = manager->FindOrBuildElement(81);
  G4Element* elSi = manager->FindOrBuildElement(14);
  G4Element* elO = manager->FindOrBuildElement(8);
  G4Element* elK = manager->FindOrBuildElement(19);
  G4Element* elCs = manager->FindOrBuildElement(55);
  G4Element* elSb = manager->FindOrBuildElement(51);

  G4Material* Vaccum = manager->FindOrBuildMaterial("G4_Galactic");
  G4Material* Argon = manager->FindOrBuildMaterial("G4_Ar");
  G4Material* CarbonDioxide = manager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4Material* Copper = manager->FindOrBuildMaterial("G4_Cu");
  G4Material* Kapton = manager->FindOrBuildMaterial("G4_KAPTON");
  G4Material* NaI = manager->FindOrBuildMaterial("G4_SODIUM_IODIDE");

// KCsSb 
  G4double mixtureDensity1 = (0.862 * 33.3/100.0 + 1.873* 33.3/100.0 + 6.691*33.3/100.0);
  G4Material* KCsSb = new G4Material("PC", mixtureDensity1, numel=3);
  KCsSb->AddElement(elK, 1); 
  KCsSb->AddElement(elCs, 1); 
  KCsSb->AddElement(elSb, 1); 

//SiO2
  density = 2.200*g/cm3;
  G4Material* SiO2 =  new G4Material("quartz",density, numel=2);
  SiO2->AddElement(elSi, natoms=1);
  SiO2->AddElement(elO , natoms=2);

//Epoxy
  density = 1.2*g/cm3;
  G4Material* Epoxy = new G4Material("Epoxy" , density, numel=2);
  Epoxy->AddElement(elH, natoms=2);
  Epoxy->AddElement(elC, natoms=2);

//FR4 (Glass + Epoxy)
  density = 1.86*g/cm3;
  G4Material* FR4 = new G4Material("FR4"  , density, numel=2);
  FR4->AddMaterial(Epoxy, fractionMass=0.472);
  FR4->AddMaterial(SiO2, fractionMass=0.528);

//Gas
  G4double mixtureDensity = (Argon->GetDensity() * 70/100.0 + CarbonDioxide->GetDensity()* 30/100.0);
  G4Material *ArCO2 = new G4Material("Ar/CO2", mixtureDensity, 2);
  ArCO2->AddMaterial(Argon, 0.7);
  ArCO2->AddMaterial(CarbonDioxide, 0.3);

//Scint
  G4Material *NaITl = new G4Material("NaITl", 3.67*g/cm3, 2);
  NaITl->AddMaterial(NaI,99.6*perCent);
  NaITl->AddElement(elTl,0.4*perCent);
//
// ------------ Generate & Add Material Properties Table ------------
//

// Scint
  const G4int nEntries = 46;
  G4double PhotonEnergy[nEntries] =
             {3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV,
              
              3.2627*eV, 3.1791*eV, 3.0996*eV, 3.024*eV,  2.952*eV,
              2.8833*eV, 2.8178*eV, 2.7552*eV, 2.6953*eV, 2.6379*eV,
              2.583*eV,  2.5303*eV, 2.48*eV,   2.431*eV,  2.3843*eV,
              2.3393*eV, 2.296*eV,  2.254*eV,  2.214*eV,  2.1751*eV,
              2.1376*eV, 2.1014*eV, 2.0664*eV, 2.0325*eV, 1.9997*eV,
              1.9678*eV, 1.9372*eV, 1.9074*eV, 1.8785*eV, 1.8505*eV,
              1.8233*eV, 1.7969*eV, 1.7712*eV, 1.7462*eV, 1.722*eV,
              1.6984*eV, 1.6754*eV, 1.6531*eV};
  G4double RefractiveIndex1[nEntries] =
             {1.8523,   1.8527,   1.8535,   1.8540,
              1.8545,   1.8550,   1.8555,   1.8563,
           
               1.84208, 1.83556,  1.82965, 1.82427,  1.81936,
              1.81486, 1.81072,  1.80691, 1.80339,  1.80013,
              1.79711, 1.7943,  1.79168, 1.78923,  1.78695,
              1.78481, 1.7828,  1.78091, 1.77914,  1.77747,
              1.77589, 1.7744,  1.773, 1.77166,  1.7704,
              1.7692, 1.76807,  1.76699, 1.76596,  1.76498,
              1.76405, 1.76316, 1.76231, 1.7615, 1.76072,
              1.75998, 1.75927, 1.75859};
  G4double Absorption1[nEntries] =                  //http://hypernews.slac.stanford.edu:5090/HyperNews/geant4/get/opticalphotons/54.html
             {5.0*m, 5.0*m, 5.0*m, 5.0*m, 
              5.0*m, 5.0*m, 5.0*m, 5.0*m,

               5.0*m, 5.0*m, 5.0*m, 5.0*m, 5.0*m,
              5.0*m, 5.0*m, 5.0*m, 5.0*m, 5.0*m,
              5.0*m, 5.0*m, 5.0*m, 5.0*m, 5.0*m,
              5.0*m, 5.0*m, 5.0*m, 5.0*m, 5.0*m,
              5.0*m, 5.0*m, 5.0*m, 5.0*m, 5.0*m,
              5.0*m, 5.0*m, 5.0*m, 5.0*m, 5.0*m,
              5.0*m, 5.0*m, 5.0*m, 5.0*m, 5.0*m,
              5.0*m, 5.0*m, 5.0*m};
  G4double ScintilFast[nEntries] =
           // { 0.019443816, 0.004432, 0.053991, 0.241971, 0.398942, 
           //   0.000134, 0.058783631 ,0.02713094,
            { 0.100481, 0.0921945, 0.085731856, 0.0821045, 0.077436129, 
              0.058783631, 0.045218178 ,0.02713094,
             //{0.081731856, 0.100158264, 0.108523626, 0.113045444, 0.11304544,
              0.103731856, 0.106158264, 0.108523626, 0.113045444, 0.11304544,
              0.104001809, 0.08907981, 0.077436129, 0.058783631, 0.045218178,
              0.03572236, 0.027130907, 0.019443816, 0.013000226, 0.009156681,
              0.004521818, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0};

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
  myMPT1->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,     nEntries);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD",38000./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);

  NaITl->SetMaterialPropertiesTable(myMPT1);
  NaITl->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4double refractiveIndex2[] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00};

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PhotonEnergy, refractiveIndex2, nEntries);

  Vaccum->SetMaterialPropertiesTable(myMPT2);

  const G4int nEntries2 = 13; 
  G4double photonEnergy_PC[nEntries2] = 
        { 3.265*eV,  3.0114*eV, 2.9540*eV, 2.8198*eV, 2.6398*eV, 2.4814*eV, 2.3409*eV, 
           2.2155*eV, 2.1029*eV, 2.001*eV, 1.9539*eV, 1.9087*eV, 1.8246*eV};
  G4double refractiveIndex_PC[nEntries2] = 
        { 1.92, 2.38, 2.61, 2.7, 3., 3., 3.23, 3.12, 3.01, 2.96, 2.95, 2.95, 2.96};
  G4double IrefractiveIndex_PC[nEntries2] = 
        { 1.69, 1.71, 1.53, 1.5, 1.34, 1.06, 0.86, 0.53, 0.42, 0.37, 0.35, 0.34, 0.33};
  G4double eff_SPC[nEntries2] = {1. , 1. ,1. ,1. ,1.
                    ,1. ,1.,1. ,1. ,1. ,1. ,1. ,1.};
  //G4double eff_SPC[nEntries2] = {0.302855808 ,0.272603508 ,0.256462766 ,0.239815267 ,0.196891892 
  //                  ,0.157011213 ,0.088752156 ,0.043294997 ,0.002349051 ,0.002200259 ,0.0 ,0.0 ,0.0};

  G4MaterialPropertiesTable* myPC = new G4MaterialPropertiesTable();
  myPC->AddProperty("IRINDEX", photonEnergy_PC, IrefractiveIndex_PC, nEntries2);
  myPC->AddProperty("RINDEX", photonEnergy_PC, refractiveIndex_PC, nEntries2);
  myPC->AddProperty("EFFICIENCY", photonEnergy_PC, eff_SPC, nEntries2);
  KCsSb->SetMaterialPropertiesTable(myPC);



//
// ------------ Visual ---------------
//
  G4VisAttributes *gemAttributes = new G4VisAttributes(G4Color::Green()) ;
  gemAttributes->SetForceWireframe(true);
  G4VisAttributes *GasGap = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.5)) ;
  GasGap->SetForceSolid(true) ;
  G4VisAttributes *PCathode = new G4VisAttributes(G4Color::Yellow()) ;
  //PCathode->SetForceWireframe(true);
  PCathode->SetForceSolid(true) ;

      
// ----------- SD -------------------
  G4SDManager* sdman = G4SDManager::GetSDMpointer();  
  GasGapSensitiveDetector* gg_sensitive = new GasGapSensitiveDetector("/PhotoCathode");
  sdman->AddNewDetector(gg_sensitive);

// ------------- Volumes --------------

// The experimental Hall ---------------------
//
  G4Box* expHall_box = new G4Box("World",0.5*fExpHall_x,0.5*fExpHall_y,0.5*fExpHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,Vaccum,"World",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

  G4VisAttributes *worldAttributes = new G4VisAttributes;
  worldAttributes->SetVisibility(true);
  expHall_log->SetVisAttributes(worldAttributes); 

// The Scintillator --------------------------
  G4cout << "Scint ZTRANS " << ZTrans << G4endl ;
  G4Box* Scint_box = new G4Box("Scint",fGem_x,fGem_y,fScint_z);
  ZTrans += Scint_box->GetZHalfLength();

  G4LogicalVolume* Scint_log
    = new G4LogicalVolume(Scint_box,NaITl,"Scint",0,0,0);
   // = new G4LogicalVolume(Scint_box,water,"Scint",0,0,0);

  G4VPhysicalVolume* Scint_phys
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,ZTrans),Scint_log,"Scint",
                        expHall_log,false,0);
  Scint_log->SetVisAttributes(new G4VisAttributes(*gemAttributes));
 // Scint_log->SetSensitiveDetector(gg_sensitive);
  ZTrans += Scint_box->GetZHalfLength();

// PC ----------------------------------------
//
  G4Box* PC_box = new G4Box("PC",fGem_x,fGem_y,fPC_z);
  ZTrans += PC_box->GetZHalfLength();

  G4LogicalVolume* PC_log
    = new G4LogicalVolume(PC_box,KCsSb,"PC",0,0,0);

  G4VPhysicalVolume* PC_phys
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,ZTrans),PC_log,"PC",
                        expHall_log,false,0);
  PC_log->SetVisAttributes(new G4VisAttributes(*PCathode));
  PC_log->SetSensitiveDetector(gg_sensitive);
  ZTrans += PC_box->GetZHalfLength();

// GEM -----------------------------------------

  std::string layerName[17] =
  {
     "GasGap1",
     "Gem1Copper1", "Gem1", "Gem1Copper2",
     "GasGap2",
     "Gem2Copper1", "Gem2", "Gem2Copper2",
     "GasGap3",
     "Gem3Copper1", "Gem3", "Gem3Copper2",
     "GasGap4",
     "ReadCopper1", "ReadoutBoard", "ReadCopper2",
     "FakeTop"
  };

 std::string layerNameLog[17];

  for (size_t i=0; i<17; i++) {
     layerNameLog[i]=layerName[i]+"Log";
  }

  G4Material* layerComp[17]=
 {
    ArCO2,                    //Drift Gap
    Copper,Kapton,Copper,      //gem1
    ArCO2,                    //Transfer I Gap
    Copper,Kapton,Copper,      //gem2
    ArCO2,                    //Transfer II Gap
    Copper,Kapton,Copper,      //gem3
    ArCO2,                    //Induction Gap
    Copper,FR4,Copper,      //Readout Board
    Vaccum                     //Fake
  };

  G4double gemSizeZ[17] =
  {
    2.*mm,                     //Drift Gap
    5.*um,50.*um,5.*um,         //gem1
    2.*mm,                     //Transfer I Gap
    5.*um,50.*um,5.*um,         //gem2
    2.*mm,                     //Transfer II Gap
    5.*um,50.*um,5.*um,        //gem3
    3.*mm,                     //Induction Gap
    35.*um,3.2*mm,35.*um,      //Readout Board
    0.1*mm                     //Fake
  };

  G4Box* SolidGem;
  G4LogicalVolume * LogicGem;
  for (G4int lyr=0; lyr<17; lyr++) {
    SolidGem = new G4Box(layerName[lyr], fGem_x, fGem_y, gemSizeZ[lyr]);
    LogicGem = new G4LogicalVolume(SolidGem, layerComp[lyr], layerNameLog[lyr]);
    LogicGem->SetVisAttributes(new G4VisAttributes(*gemAttributes));
    gemCollection.push_back(SolidGem);
    gemLogCollection.push_back(LogicGem);
  }

  PlaceGeometry(G4ThreeVector(0.,0.,0.), expHall_log, ZTrans);
  return expHall_phys;
}

void DetectorConstruction::PlaceGeometry(G4ThreeVector tlate, G4LogicalVolume* pMotherLogical, G4double Ztrans) {

  G4double ZTranslation = Ztrans;

  for (size_t i=0; i<gemCollection.size(); i++) {
     ZTranslation += gemCollection.at(i)->GetZHalfLength();
     G4ThreeVector position = tlate + G4ThreeVector(0,0,ZTranslation);
     new G4PVPlacement(0,
                       position,
                       gemLogCollection.at(i),
                       gemCollection.at(i)->GetName(),
                       pMotherLogical,
                       false,
                       i);
     ZTranslation += gemCollection.at(i)->GetZHalfLength();
  }
}
