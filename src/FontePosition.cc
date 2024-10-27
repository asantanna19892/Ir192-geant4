#include "globals.hh"
#include "G4RunManager.hh"
#include "FontePosition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
FontePosition::FontePosition()
{;}


FontePosition::~FontePosition()
{;}

G4ThreeVector FontePosition::GetPosition()
{
//G4ThreeVector posFonte(0,0,0);9,-30,285

G4double X,Y,Z;
G4double I,J,K;


//  converte XYZ para IJK

/*

X=9*mm;
Y=-30*mm;
Z=285*mm;





I=(192/2)+(X/3.6)    +1;
J=(96/2)+(X/3.6)     +1;
K=(498/2)+(X/3.6)  +1;

G4cout << I <<" , "<< J <<" , "<< K << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<G4endl;

// Converte de IJK para XYZ;99 , 51 , 252

I=99;
J=51;
K=252;


X= +(I*3.6)+1.8 -((192*3.6)/2) ;
Y= -((96*3.6/2))+(J*3.6)+1.8 ;
Z= -((498*3.6)/2)+(K*3.6)+1.8 ;


G4cout << X <<" , "<< Y <<" , "<< Z << "##############################################################"<<G4endl;
*/


G4ThreeVector posFonte (0*mm,0*mm,0*mm);
//G4ThreeVector posFonte (12.6*mm,36*mm,-7.2*mm);

return posFonte;
}

