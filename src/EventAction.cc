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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdepAnel1(0.),
  fEdepAnel2(0.),
  fEdepAnel3(0.),
  fEdepAnel4(0.),
  fEdepAnel5(0.),
  fEdepAnel6(0.),
  fEdepAnel7(0.),
  fEdepAnel8(0.),
  fEdepAnel9(0.),
  fEdepAnel10(0.),
  fEdepAnel11(0.),
  fEdepAnel12(0.)
  

{
  
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{    
  // Inicializar a energia depositada em ambos os anéis
  fEdepAnel1 = 0.;
  fEdepAnel2 = 0.;
  fEdepAnel3 = 0.;
  fEdepAnel4 = 0.;
  fEdepAnel5 = 0.;
  fEdepAnel6 = 0.;
  fEdepAnel7 = 0.;
  fEdepAnel8 = 0.;
  fEdepAnel9 = 0.;
  fEdepAnel10 = 0.;
  fEdepAnel11 = 0.;
  fEdepAnel12 = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{   
  // Acumular estatísticas separadamente para os dois anéis

  fRunAction->AddEdepAnel1(fEdepAnel1);  // Acumular energia do Anel 1
  fRunAction->AddEdepAnel2(fEdepAnel2);  // Acumular energia do Anel 2
  fRunAction->AddEdepAnel3(fEdepAnel3);  // Acumular energia do Anel 1
  fRunAction->AddEdepAnel4(fEdepAnel4);  // Acumular energia do Anel 2
  fRunAction->AddEdepAnel5(fEdepAnel5);  // Acumular energia do Anel 1
  fRunAction->AddEdepAnel6(fEdepAnel6);  // Acumular energia do Anel 2
  fRunAction->AddEdepAnel7(fEdepAnel7);  // Acumular energia do Anel 1
  fRunAction->AddEdepAnel8(fEdepAnel8);  // Acumular energia do Anel 2
  fRunAction->AddEdepAnel9(fEdepAnel9);  // Acumular energia do Anel 1
  fRunAction->AddEdepAnel10(fEdepAnel10);  // Acumular energia do Anel 2
  fRunAction->AddEdepAnel11(fEdepAnel11);  // Acumular energia do Anel 1
  fRunAction->AddEdepAnel12(fEdepAnel12);  // Acumular energia do Anel 2
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

