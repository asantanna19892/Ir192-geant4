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
    G4LogicalVolume* GetScoringVolume1() const { return fScoringVolume1; };
    G4LogicalVolume* GetScoringVolume2() const { return fScoringVolume2; };
    G4LogicalVolume* GetScoringVolume3() const { return fScoringVolume3; };
    G4LogicalVolume* GetScoringVolume4() const { return fScoringVolume4; };
    G4LogicalVolume* GetScoringVolume5() const { return fScoringVolume5; };
    G4LogicalVolume* GetScoringVolume6() const { return fScoringVolume6; };
    G4LogicalVolume* GetScoringVolume7() const { return fScoringVolume7; };
    G4LogicalVolume* GetScoringVolume8() const { return fScoringVolume8; };
    G4LogicalVolume* GetScoringVolume9() const { return fScoringVolume9; };
    G4LogicalVolume* GetScoringVolume10() const { return fScoringVolume10; };
    G4LogicalVolume* GetScoringVolume11() const { return fScoringVolume11; };
    G4LogicalVolume* GetScoringVolume12() const { return fScoringVolume12; }
    //G4VPhysicalVolume* GetScoringVolume() const { return fScoringVolume; }
  private:
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


