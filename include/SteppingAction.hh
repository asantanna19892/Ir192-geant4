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
// $Id: ExN03SteppingAction.hh,v 1.8 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "globals.hh"

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4VUserDetectorConstruction.hh"

class EventAction;
class G4LogicalVolume;
class DetectorConstruction;
class FontePosition;

//class EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(EventAction* eventAction);
    virtual ~SteppingAction();

    void UserSteppingAction(const G4Step*);
    G4double Myrmy[101];
  


private:
  DetectorConstruction* detector;
  FontePosition* PosInicialFonte;

  EventAction*  fEventAction;
  //G4LogicalVolume* fScoringVolume;
  G4LogicalVolume* fScoringVolume1;
  G4LogicalVolume* fScoringVolume2;
  G4LogicalVolume* fScoringVolume3;
  G4LogicalVolume* fScoringVolume4;
  G4LogicalVolume* fScoringVolume5;
  G4LogicalVolume* fScoringVolume6;
  G4LogicalVolume* fScoringVolume7;
  G4LogicalVolume* fScoringVolume8;
  G4LogicalVolume* fScoringVolume9;
  G4LogicalVolume* fScoringVolume10;
  G4LogicalVolume* fScoringVolume11;
  G4LogicalVolume* fScoringVolume12;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
