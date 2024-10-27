#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class FontePosition;
class G4ParticleGun;
class G4Event;
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction();
  // methods
  virtual void GeneratePrimaries(G4Event*);

    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
  

private:
  // data members
  G4ParticleGun* fParticleGun; //pointer a to G4 service class
  FontePosition* PosInicialFonte;
};
#endif
