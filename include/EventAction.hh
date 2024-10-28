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
//
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;

/// Event action class
///

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    // Métodos para adicionar energia nos dois anéis
    void AddEdepAnel1(G4double edep) { fEdepAnel1 += edep; }
    void AddEdepAnel2(G4double edep) { fEdepAnel2 += edep; }
    void AddEdepAnel3(G4double edep) { fEdepAnel3 += edep; }
    void AddEdepAnel4(G4double edep) { fEdepAnel4 += edep; }
    void AddEdepAnel5(G4double edep) { fEdepAnel5 += edep; }
    void AddEdepAnel6(G4double edep) { fEdepAnel6 += edep; }
    void AddEdepAnel7(G4double edep) { fEdepAnel7 += edep; }
    void AddEdepAnel8(G4double edep) { fEdepAnel8 += edep; }
    void AddEdepAnel9(G4double edep) { fEdepAnel9 += edep; }
    void AddEdepAnel10(G4double edep) { fEdepAnel10 += edep; }
    void AddEdepAnel11(G4double edep) { fEdepAnel11 += edep; }
    void AddEdepAnel12(G4double edep) { fEdepAnel12 += edep; }
  private:
    RunAction* fRunAction;
    // Energia depositada em cada anel
    G4double fEdepAnel1;
    G4double fEdepAnel2;
    G4double fEdepAnel3;
    G4double fEdepAnel4;
    G4double fEdepAnel5;
    G4double fEdepAnel6;
    G4double fEdepAnel7;
    G4double fEdepAnel8;
    G4double fEdepAnel9;
    G4double fEdepAnel10;
    G4double fEdepAnel11;
    G4double fEdepAnel12;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
