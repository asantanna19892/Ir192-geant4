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
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    // Métodos para adicionar energia nos dois anéis
    void AddEdepAnel1(G4double edep);
    void AddEdepAnel2(G4double edep);
    void AddEdepAnel3(G4double edep);
    void AddEdepAnel4(G4double edep);
    void AddEdepAnel5(G4double edep);
    void AddEdepAnel6(G4double edep);
    void AddEdepAnel7(G4double edep);
    void AddEdepAnel8(G4double edep);
    void AddEdepAnel9(G4double edep);
    void AddEdepAnel10(G4double edep);
    void AddEdepAnel11(G4double edep);
    void AddEdepAnel12(G4double edep);

  private:
    G4double fEdepAnel1; // Energia acumulada no primeiro anel
    G4double fEdep2Anel1; // Variância do primeiro anel
    
    G4double fEdepAnel2;  // Nova variável para o segundo anel
    G4double fEdep2Anel2; // Nova variável para o segundo anel

    G4double fEdepAnel3;  // Nova variável para o segundo anel
    G4double fEdep2Anel3; // Nova variável para o segundo anel

    G4double fEdepAnel4;  // Nova variável para o segundo anel
    G4double fEdep2Anel4; // Nova variável para o segundo anel

    G4double fEdepAnel5;  // Nova variável para o segundo anel
    G4double fEdep2Anel5; // Nova variável para o segundo anel

    G4double fEdepAnel6;  // Nova variável para o segundo anel
    G4double fEdep2Anel6; // Nova variável para o segundo anel

    G4double fEdepAnel7;  // Nova variável para o segundo anel
    G4double fEdep2Anel7; // Nova variável para o segundo anel

    G4double fEdepAnel8;  // Nova variável para o segundo anel
    G4double fEdep2Anel8; // Nova variável para o segundo anel

    G4double fEdepAnel9;  // Nova variável para o segundo anel
    G4double fEdep2Anel9; // Nova variável para o segundo anel

    G4double fEdepAnel10;  // Nova variável para o segundo anel
    G4double fEdep2Anel10; // Nova variável para o segundo anel

    G4double fEdepAnel11;  // Nova variável para o segundo anel
    G4double fEdep2Anel11; // Nova variável para o segundo anel

    G4double fEdepAnel12;  // Nova variável para o segundo anel
    G4double fEdep2Anel12; // Nova variável para o segundo anel

};

#endif


