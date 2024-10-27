#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class FontePosition;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    FontePosition* PosInicialFonte;
    G4LogicalVolume* GetScoringVolume1() const { return fScoringVolume1; }
    G4LogicalVolume* GetScoringVolume2() const { return fScoringVolume2; }
    //G4VPhysicalVolume* GetScoringVolume() const { return fScoringVolume; }
  private:
    G4LogicalVolume* fScoringVolume1;
    G4LogicalVolume* fScoringVolume2;
};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


