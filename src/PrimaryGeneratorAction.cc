#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "FontePosition.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::PrimaryGeneratorAction(): 
    G4VUserPrimaryGeneratorAction(), fParticleGun(0)

{
G4int nofParticles = 1;

 fParticleGun = new G4ParticleGun(nofParticles);
 G4String particleName;
    //Definindo a tabela de particulas
G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma");
//G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");

fParticleGun->SetParticleDefinition(particle);
//fParticleGun->SetParticleEnergy(1*MeV);
fParticleGun->SetParticleEnergy(370*keV);
fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.5,1));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
delete fParticleGun;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
// this function is called at the begining of event

G4ThreeVector posIni = PosInicialFonte->GetPosition();

//*********************************************************************
    srand ( time(NULL) ); 
    //Definindo a posição das particulas
    //sorteio da posição     
    G4double x,y,z;
    G4double raio=0.5*mm;
    do{
        x = (G4UniformRand()-0.5)*2*raio;
        y = (G4UniformRand()-0.5)*2*raio;
      }while(x*x+y*y > raio*raio);
        z = (G4UniformRand()-0.5)*1.3*mm;

	  G4ThreeVector position(x,y,z);
      fParticleGun->SetParticlePosition(position+posIni);

   
//**********************************************************************
    srand ( time(NULL) );
    //Definindo o momento das particulas
    //sorteio da direçãoo i, j, k
    G4double a,b,c;
    G4double n;
    do{
       a = (G4UniformRand()-0.5)/0.5;
       b = (G4UniformRand()-0.5)/0.5; 
       c = (G4UniformRand()-0.5)/0.5;
       n = a*a+b*b+c*c;
      }while(n>1 || n == 0.0);
       n = sqrt(n);
       a /= n;//normalizando os valores para i
       b /= n;//       "     "     "     "   j
       c /= n;//       "     "     "     "   k 

      G4ThreeVector direction(a,b,c);
      fParticleGun->SetParticleMomentumDirection(direction); 


    // G4double primaryParticleEnergy;
    // G4int i = 0;

    // G4double E1[] = {0*keV,8*keV,9*keV,10*keV,11*keV,12*keV,13*keV,14*keV,62*keV,63*keV,64*keV,65*keV,66*keV,67*keV,72*keV,73*keV,74*keV,75*keV,76*keV,77*keV,78*keV,79*keV,137*keV,202*keV,206*keV,284*keV,296*keV,309*keV,317*keV,375*keV,417*keV,469*keV,485*keV,490*keV,589*keV,605*keV,613*keV,885*keV};

    // G4double soma[] = {0.0,0.003506,0.007012,0.010518,0.014024,0.017530,0.021036,0.024542,0.032102,0.039662,0.047222,0.054781,0.062341,0.069901,0.071431,0.072960,0.074490,0.076020,0.077549,0.079079,0.080608,0.082138,0.082904,0.084956,0.099046,0.100172,0.222245,0.349395,0.699747,0.702798,0.705607,0.907863,0.921234,0.923041,0.941997,0.976525,0.998781,1.000000};

    // G4double r = G4UniformRand();
 
    //     while (r>=soma[i]){
    //       i++;
    //     }

    //   G4double e = E1[i];
    //   primaryParticleEnergy = e;
   //   fParticleGun->SetParticleEnergy(primaryParticleEnergy);


fParticleGun->GeneratePrimaryVertex(anEvent);
}
