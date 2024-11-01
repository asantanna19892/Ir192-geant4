#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4Run.hh"

RunAction::RunAction()
  : G4UserRunAction(),
    fEdepAnel1(0.), fEdep2Anel1(0.),
    fEdepAnel2(0.), fEdep2Anel2(0.),
    fEdepAnel3(0.), fEdep2Anel3(0.),
    fEdepAnel4(0.), fEdep2Anel4(0.),
    fEdepAnel5(0.), fEdep2Anel5(0.),
    fEdepAnel6(0.), fEdep2Anel6(0.),
    fEdepAnel7(0.), fEdep2Anel7(0.),
    fEdepAnel8(0.), fEdep2Anel8(0.),
    fEdepAnel9(0.), fEdep2Anel9(0.),
    fEdepAnel10(0.), fEdep2Anel10(0.),
    fEdepAnel11(0.), fEdep2Anel11(0.),
    fEdepAnel12(0.), fEdep2Anel12(0.)
{
  // Adiciona novas unidades para dose
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy", "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy", "Dose", microgray);
  new G4UnitDefinition("nanogray", "nanoGy", "Dose", nanogray);
  new G4UnitDefinition("picogray", "picoGy", "Dose", picogray); 
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run*)
{
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  fEdepAnel1 = 0.;
  fEdep2Anel1 = 0.;
  fEdepAnel2 = 0.;
  fEdep2Anel2 = 0.;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Calcular a dose para anel 1
  G4double edepAnel1 = fEdepAnel1;
  G4double edep2Anel1 = fEdep2Anel1;
  G4double rmsAnel1 = edep2Anel1 - edepAnel1*edepAnel1/nofEvents;
  if (rmsAnel1 > 0.) rmsAnel1 = std::sqrt(rmsAnel1); else rmsAnel1 = 0.;

  // Calcular a dose para anel 2
  G4double edepAnel2 = fEdepAnel2;
  G4double edep2Anel2 = fEdep2Anel2;
  G4double rmsAnel2 = edep2Anel2 - edepAnel2*edepAnel2/nofEvents;
  if (rmsAnel2 > 0.) rmsAnel2 = std::sqrt(rmsAnel2); else rmsAnel2 = 0.;

  // Calcular a dose para anel 3
  G4double edepAnel3 = fEdepAnel3;
  G4double edep2Anel3 = fEdep2Anel3;
  G4double rmsAnel3 = edep2Anel3 - edepAnel3*edepAnel3/nofEvents;
  if (rmsAnel3 > 0.) rmsAnel3 = std::sqrt(rmsAnel3); else rmsAnel3 = 0.;

  // Calcular a dose para anel 4
  G4double edepAnel4 = fEdepAnel4;
  G4double edep2Anel4 = fEdep2Anel4;
  G4double rmsAnel4 = edep2Anel4 - edepAnel4*edepAnel4/nofEvents;
  if (rmsAnel4 > 0.) rmsAnel4 = std::sqrt(rmsAnel4); else rmsAnel4 = 0.;

  // Calcular a dose para anel 5
  G4double edepAnel5 = fEdepAnel5;
  G4double edep2Anel5 = fEdep2Anel5;
  G4double rmsAnel5 = edep2Anel5 - edepAnel5*edepAnel5/nofEvents;
  if (rmsAnel5 > 0.) rmsAnel5 = std::sqrt(rmsAnel5); else rmsAnel5 = 0.;

  // Calcular a dose para anel 6
  G4double edepAnel6 = fEdepAnel6;
  G4double edep2Anel6 = fEdep2Anel6;
  G4double rmsAnel6 = edep2Anel6 - edepAnel6*edepAnel6/nofEvents;
  if (rmsAnel6 > 0.) rmsAnel6 = std::sqrt(rmsAnel6); else rmsAnel6 = 0.;

  // Calcular a dose para anel 7
  G4double edepAnel7 = fEdepAnel7;
  G4double edep2Anel7 = fEdep2Anel7;
  G4double rmsAnel7 = edep2Anel7 - edepAnel7*edepAnel7/nofEvents;
  if (rmsAnel7 > 0.) rmsAnel7 = std::sqrt(rmsAnel7); else rmsAnel7 = 0.;

  // Calcular a dose para anel 8
  G4double edepAnel8 = fEdepAnel8;
  G4double edep2Anel8 = fEdep2Anel8;
  G4double rmsAnel8 = edep2Anel8 - edepAnel8*edepAnel8/nofEvents;
  if (rmsAnel8 > 0.) rmsAnel8 = std::sqrt(rmsAnel8); else rmsAnel8 = 0.;

  // Calcular a dose para anel 9
  G4double edepAnel9 = fEdepAnel9;
  G4double edep2Anel9 = fEdep2Anel9;
  G4double rmsAnel9 = edep2Anel9 - edepAnel9*edepAnel9/nofEvents;
  if (rmsAnel9 > 0.) rmsAnel9 = std::sqrt(rmsAnel9); else rmsAnel9 = 0.;

  // Calcular a dose para anel 10
  G4double edepAnel10 = fEdepAnel10;
  G4double edep2Anel10 = fEdep2Anel10;
  G4double rmsAnel10 = edep2Anel10 - edepAnel10*edepAnel10/nofEvents;
  if (rmsAnel10 > 0.) rmsAnel10 = std::sqrt(rmsAnel10); else rmsAnel10 = 0.;

  // Calcular a dose para anel 11
  G4double edepAnel11 = fEdepAnel11;
  G4double edep2Anel11 = fEdep2Anel11;
  G4double rmsAnel11 = edep2Anel11 - edepAnel11*edepAnel11/nofEvents;
  if (rmsAnel11 > 0.) rmsAnel11 = std::sqrt(rmsAnel11); else rmsAnel11 = 0.;

  // Calcular a dose para anel 12
  G4double edepAnel12 = fEdepAnel12;
  G4double edep2Anel12 = fEdep2Anel12;
  G4double rmsAnel12 = edep2Anel12 - edepAnel12*edepAnel12/nofEvents;
  if (rmsAnel12 > 0.) rmsAnel12 = std::sqrt(rmsAnel12); else rmsAnel12 = 0.;

  // Obter a massa do volume de scoring
  const DetectorConstruction* detectorConstruction
    = static_cast<const DetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double massAnel1 = detectorConstruction->GetScoringVolume1()->GetMass();
  G4double massAnel2 = detectorConstruction->GetScoringVolume2()->GetMass();
  G4double massAnel3 = detectorConstruction->GetScoringVolume3()->GetMass();
  G4double massAnel4 = detectorConstruction->GetScoringVolume4()->GetMass();
  G4double massAnel5 = detectorConstruction->GetScoringVolume5()->GetMass();
  G4double massAnel6 = detectorConstruction->GetScoringVolume6()->GetMass();
  G4double massAnel7 = detectorConstruction->GetScoringVolume7()->GetMass();
  G4double massAnel8 = detectorConstruction->GetScoringVolume8()->GetMass();
  G4double massAnel9 = detectorConstruction->GetScoringVolume9()->GetMass();
  G4double massAnel10 = detectorConstruction->GetScoringVolume10()->GetMass();
  G4double massAnel11 = detectorConstruction->GetScoringVolume11()->GetMass();
  G4double massAnel12 = detectorConstruction->GetScoringVolume12()->GetMass();

  
  
  G4double doseAnel1 = edepAnel1/massAnel1;
  G4double rmsDoseAnel1 = rmsAnel1/massAnel1;

  G4double doseAnel2 = edepAnel2/massAnel2; 
  G4double rmsDoseAnel2 = rmsAnel2/massAnel2;

  G4double doseAnel3 = edepAnel1/massAnel3;
  G4double rmsDoseAnel3 = rmsAnel1/massAnel3;

  G4double doseAnel4 = edepAnel4/massAnel4;
  G4double rmsDoseAnel4 = rmsAnel4/massAnel4;

  G4double doseAnel5 = edepAnel5/massAnel5;
  G4double rmsDoseAnel5 = rmsAnel5/massAnel5;

  G4double doseAnel6 = edepAnel6/massAnel6; 
  G4double rmsDoseAnel6 = rmsAnel6/massAnel6;


  G4double doseAnel7 = edepAnel7/massAnel7;
  G4double rmsDoseAnel7 = rmsAnel7/massAnel7;

  G4double doseAnel8 = edepAnel8/massAnel8; 
  G4double rmsDoseAnel8 = rmsAnel8/massAnel8;

  G4double doseAnel9 = edepAnel9/massAnel9;
  G4double rmsDoseAnel9 = rmsAnel9/massAnel9;

  G4double doseAnel10 = edepAnel10/massAnel10;
  G4double rmsDoseAnel10 = rmsAnel10/massAnel10;

  G4double doseAnel11 = edepAnel11/massAnel11;
  G4double rmsDoseAnel11 = rmsAnel11/massAnel11;

  G4double doseAnel12 = edepAnel12/massAnel12; 
  G4double rmsDoseAnel12 = rmsAnel12/massAnel12;

  const PrimaryGeneratorAction* generatorAction
    = static_cast<const PrimaryGeneratorAction*>
    (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
    {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
  // Imprimir as doses calculadas para ambos os anéis
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " events."
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 1): " 
     << G4BestUnit(doseAnel1,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel1,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 2): " 
     << G4BestUnit(doseAnel2,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel2,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 3): " 
     << G4BestUnit(doseAnel3,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel3,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 4): " 
     << G4BestUnit(doseAnel4,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel4,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 5): " 
     << G4BestUnit(doseAnel5,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel5,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 6): " 
     << G4BestUnit(doseAnel6,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel6,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 7): " 
     << G4BestUnit(doseAnel7,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel7,"Dose")
     << G4endl     
     << " Cumulated dose per run, in scoring volume (Anel 8): " 
     << G4BestUnit(doseAnel8,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel8,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 9): " 
     << G4BestUnit(doseAnel9,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel9,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 10): " 
     << G4BestUnit(doseAnel10,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel10,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 11): " 
     << G4BestUnit(doseAnel11,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel11,"Dose")
     << G4endl
     << " Cumulated dose per run, in scoring volume (Anel 12): " 
     << G4BestUnit(doseAnel12,"Dose") << " rms = " << G4BestUnit(rmsDoseAnel12,"Dose")
     << G4endl;
     
}

void RunAction::AddEdepAnel1(G4double edep)
{
  fEdepAnel1 += edep;
  fEdep2Anel1 += edep*edep;
}

void RunAction::AddEdepAnel2(G4double edep)
{
  fEdepAnel2 += edep;
  fEdep2Anel2 += edep*edep;
}


void RunAction::AddEdepAnel3(G4double edep)
{
  fEdepAnel3 += edep;
  fEdep2Anel3 += edep*edep;
}

void RunAction::AddEdepAnel4(G4double edep)
{
  fEdepAnel4 += edep;
  fEdep2Anel4 += edep*edep;
}

void RunAction::AddEdepAnel5(G4double edep)
{
  fEdepAnel5 += edep;
  fEdep2Anel5 += edep*edep;
}

void RunAction::AddEdepAnel6(G4double edep)
{
  fEdepAnel6 += edep;
  fEdep2Anel6 += edep*edep;
}

void RunAction::AddEdepAnel7(G4double edep)
{
  fEdepAnel7 += edep;
  fEdep2Anel7 += edep*edep;
}

void RunAction::AddEdepAnel8(G4double edep)
{
  fEdepAnel8 += edep;
  fEdep2Anel8 += edep*edep;
}


void RunAction::AddEdepAnel9(G4double edep)
{
  fEdepAnel9 += edep;
  fEdep2Anel9 += edep*edep;
}

void RunAction::AddEdepAnel10(G4double edep)
{
  fEdepAnel10 += edep;
  fEdep2Anel10 += edep*edep;
}

void RunAction::AddEdepAnel11(G4double edep)
{
  fEdepAnel11 += edep;
  fEdep2Anel11 += edep*edep;
}

void RunAction::AddEdepAnel12(G4double edep)
{
  fEdepAnel12 += edep;
  fEdep2Anel12 += edep*edep;
}
