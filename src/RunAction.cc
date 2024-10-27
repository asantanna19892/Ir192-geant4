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
    fEdepAnel2(0.), fEdep2Anel2(0.)
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

}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Calcular a dose para o primeiro anel
  G4double edepAnel1 = fEdepAnel1;
  G4double edep2Anel1 = fEdep2Anel1;
  G4double rmsAnel1 = edep2Anel1 - edepAnel1*edepAnel1/nofEvents;
  if (rmsAnel1 > 0.) rmsAnel1 = std::sqrt(rmsAnel1); else rmsAnel1 = 0.;

  // Calcular a dose para o segundo anel
  G4double edepAnel2 = fEdepAnel2;
  G4double edep2Anel2 = fEdep2Anel2;
  G4double rmsAnel2 = edep2Anel2 - edepAnel2*edepAnel2/nofEvents;
  if (rmsAnel2 > 0.) rmsAnel2 = std::sqrt(rmsAnel2); else rmsAnel2 = 0.;

  // Obter a massa do volume de scoring
  const DetectorConstruction* detectorConstruction
    = static_cast<const DetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double massAnel1 = detectorConstruction->GetScoringVolume1()->GetMass();
  G4double massAnel2 = detectorConstruction->GetScoringVolume2()->GetMass();
  
  G4double doseAnel1 = edepAnel1/massAnel1;
  G4double rmsDoseAnel1 = rmsAnel1/massAnel1;

  G4double doseAnel2 = edepAnel2/massAnel2; // Assumindo a mesma massa para simplificar
  G4double rmsDoseAnel2 = rmsAnel2/massAnel2;
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

