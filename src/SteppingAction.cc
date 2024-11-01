//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExN03SteppingAction.cc,v 1.8 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "SteppingAction.hh"
#include "G4Material.hh"
#include "DetectorConstruction.hh"
#include "G4Track.hh"
#include "FontePosition.hh"
#include "G4UnitsTable.hh"
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif


G4int conta=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



SteppingAction::SteppingAction(EventAction* eventAction)
  : G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume1(0),
  fScoringVolume2(0),

  // extrapolar todos os aneis
  fScoringVolume3(0),
  fScoringVolume4(0),
  fScoringVolume5(0),
  fScoringVolume6(0),
  fScoringVolume7(0),
  fScoringVolume8(0),
  fScoringVolume9(0),
  fScoringVolume10(0),
  fScoringVolume11(0),
  fScoringVolume12(0)


{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
     if (!fScoringVolume1 || !fScoringVolume2 || !fScoringVolume3 || !fScoringVolume4 || !fScoringVolume5 || !fScoringVolume6 || !fScoringVolume7|| !fScoringVolume8 || !fScoringVolume9 || !fScoringVolume10 || !fScoringVolume11 || !fScoringVolume12) { 
          const DetectorConstruction* detectorConstruction
          = static_cast<const DetectorConstruction*>
               (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
          
          fScoringVolume1 = detectorConstruction->GetScoringVolume1();
          fScoringVolume2 = detectorConstruction->GetScoringVolume2();
          fScoringVolume3 = detectorConstruction->GetScoringVolume3();
          fScoringVolume4 = detectorConstruction->GetScoringVolume4();
          fScoringVolume5 = detectorConstruction->GetScoringVolume5();
          fScoringVolume6 = detectorConstruction->GetScoringVolume6();
          fScoringVolume7 = detectorConstruction->GetScoringVolume7();
          fScoringVolume8 = detectorConstruction->GetScoringVolume8();
          fScoringVolume9 = detectorConstruction->GetScoringVolume9();
          fScoringVolume10 = detectorConstruction->GetScoringVolume10();
          fScoringVolume11 = detectorConstruction->GetScoringVolume11();
          fScoringVolume12 = detectorConstruction->GetScoringVolume12();
     }
     // get volume of the current step
          G4LogicalVolume* volume 
               = step->GetPreStepPoint()->GetTouchableHandle()
                    ->GetVolume()->GetLogicalVolume();

     // check if we are in scoring volume
     if (volume != fScoringVolume1 && volume != fScoringVolume2 && volume != fScoringVolume3  && volume != fScoringVolume4 && volume != fScoringVolume5 && volume != fScoringVolume6 && volume != fScoringVolume7 && volume != fScoringVolume8 && volume != fScoringVolume9 && volume != fScoringVolume10 && volume != fScoringVolume11 && volume != fScoringVolume12) return;

     // collect energy deposited in this step
     G4double edepStep = step->GetTotalEnergyDeposit();
     

     if (volume == fScoringVolume1) {
        fEventAction->AddEdepAnel1(edepStep);  // Acumula energia no primeiro anel
     } else if (volume == fScoringVolume2) {
        fEventAction->AddEdepAnel2(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume3) {
        fEventAction->AddEdepAnel3(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume4) {
        fEventAction->AddEdepAnel4(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume5) {
        fEventAction->AddEdepAnel5(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume6) {
        fEventAction->AddEdepAnel6(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume7) {
        fEventAction->AddEdepAnel7(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume8) {
        fEventAction->AddEdepAnel8(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume9) {
        fEventAction->AddEdepAnel9(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume10) {
        fEventAction->AddEdepAnel10(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume11) {
        fEventAction->AddEdepAnel11(edepStep);  // Acumula energia no segundo anel
     } else if (volume == fScoringVolume12) {
        fEventAction->AddEdepAnel12(edepStep);  // Acumula energia no segundo anel
     }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......






