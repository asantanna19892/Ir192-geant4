
#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Sphere.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "FontePosition.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume1(0),
  fScoringVolume2(0) 
{ }

DetectorConstruction::~DetectorConstruction()
{
  delete PosInicialFonte;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  //------------------------------------------------------ materials
  G4String symbol;
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4int ncomponents, natoms;
  G4double  fractionmass;

  //
  // define Elements
  //

  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  //G4Element* Si = new G4Element("Silicon",symbol="Si" , z= 14., a= 28.09*g/mole);

  G4Element*Na   = new G4Element("Sodium"  ,"Na" , 11,23.00*g/mole); 
  G4Element*Mg   = new G4Element("Magnesium"  ,"Mg" , 12,24.30*g/mole);
  G4Element*P    = new G4Element("Phosphorus"  ,"P" , 15,31.00*g/mole);
  G4Element*S    = new G4Element("Sulfur"  ,"S" , 16,32.10*g/mole);
  G4Element*Cl   = new G4Element("Chlorine"  ,"Cl" , 17,35.50*g/mole);
  G4Element*K    = new G4Element("Potassium"  ,"K" , 19,39.10*g/mole);
  G4Element*Ca   = new G4Element("Calcium"  ,"Ca" , 20,40.10*g/mole);
  G4Element*Fe   = new G4Element("Iron"  ,"Fe" , 26,55.90*g/mole);
  G4Element*Zn   = new G4Element("Zinc"  ,"Zn" , 30,65.50*g/mole);
  G4Element*Rb   = new G4Element("Rubidium"  ,"Rb" , 37,85.50*g/mole);
  G4Element*Sr   = new G4Element("Strontiom"  ,"Sr" , 38,87.60*g/mole);
  G4Element*Zr   = new G4Element("Zirconium"  ,"Zr" , 20,91.22*g/mole);
  G4Element*Pb   = new G4Element("Lead"  ,"Pb" , 82,207.20*g/mole);
  G4Element*I     =new G4Element("Iodine"  ,"I" , 53,127.90*g/mole);

  G4Material*H2O= 
  new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

 //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  G4Material* Air = 
  new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);


  //
  // examples of gas in non STP conditions
  //

  G4Material* CO2 = 
  new G4Material("CarbonicGas", density= 27.*mg/cm3, ncomponents=2,
                                kStateGas, 325.*kelvin, 50.*atmosphere);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);
  
  G4Material* steam = 
  new G4Material("WaterSteam", density= 0.3*mg/cm3, ncomponents=1,
                              kStateGas, 500.*kelvin, 2.*atmosphere);
  steam->AddMaterial(H2O, fractionmass=1.);


  //****************************************************
  // Iridium(medical Physical, vol 25, No 10, out 1998)
  //****************************************************


  G4Material* Ir192 = new 
  G4Material("Iridium", z=77, a=191.96260*g/mole, density=22.42*g/cm3);


  G4Material* W184 = new 
  G4Material("Tungstenio", z=74, a=184.0*g/mole, density=19.25*g/cm3);


  G4Material* I125 = new 
  G4Material("Iodine", z=53, a=124.9046242*g/mole, density=19.25*g/cm3);


  G4Material* Ti = new 
  G4Material ("Titanio", z=22, a=47.9*g/mole, density=4.54*g/cm3);


  G4Material*polystyrene = new
  G4Material("polystyrene",density = 1.06*g/cm3,ncomponents= 2);
  polystyrene->AddElement(H, 0.077);
  polystyrene->AddElement(C, 0.923);


  //********
  // Stainless Stell 
  //********

  // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
  
  
  G4Element* Mn  = new G4Element("Manganese","Mn",z = 25.,a  =  54.94*g/mole);
   
  G4Element* Si  = new G4Element("Silicon","Si",z = 14.,a = 28.09*g/mole);
 
  G4Element* Cr  = new G4Element("Chromium","Cr",z = 24.,a = 52.00*g/mole);
  
  G4Element* Ni  = new G4Element("Nickel","Ni",z = 28.,a = 58.70*g/mole);
 
  //G4Element* Fe  = new G4Element("Iron","Fe",z = 26.,a = 55.85*g/mole);


  G4Material* steelAISI321 = new G4Material("Stainless steel",density = 7.9*g/cm3,ncomponents= 5); 
  steelAISI321->AddElement(Mn, 0.02);
  steelAISI321->AddElement(Si, 0.01);
  steelAISI321->AddElement(Cr, 0.18);
  steelAISI321->AddElement(Ni, 0.10);
  steelAISI321->AddElement(Fe, 0.69);

  G4Material* steelAISI301 = new G4Material("Stainless steel",density = 5.6*g/cm3,ncomponents=5); //effective density according to Ballester (2001).
  steelAISI301->AddElement(Mn, 0.02);
  steelAISI301->AddElement(Si, 0.02);
  steelAISI301->AddElement(Cr, 0.18);
  steelAISI301->AddElement(Ni, 0.09);
  steelAISI301->AddElement(Fe, 0.69);


  //////////////////////////
  // examples of vacuum
  /////////////////////////



  G4Material* Vacuum =
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,

                            kStateGas, 2.73*kelvin, 3.e-18*pascal);

  G4Material* beam = 
  new G4Material("Beam", density= 1.e-5*g/cm3, ncomponents=1,
                        kStateGas, STP_Temperature, 2.e-2*bar);
  beam->AddMaterial(Air, fractionmass=1.);

  G4ThreeVector posIni = PosInicialFonte->GetPosition();
  G4double posFonteZ = posIni.z();
  G4double posFonteY = posIni.y();
  G4double posFonteX = posIni.x();  
  //     
  // World
  //
  

  G4double world_sizeXY = 0.5*m;
  G4double world_sizeZ  = 0.5*m;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        H2O,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  // G4double pRMin = 12*cm;
  // G4double pRMax = 20*cm;
  // G4double pDz = 30*cm;
  // G4double pSPhi = 0;
  // G4double pDPhi = 1.5*pi;
  // G4ThreeVector pLowNorm;
  // pLowNorm = G4ThreeVector(0,0,-1);
  // G4ThreeVector pHighNorm = G4ThreeVector(0,0,1);

  // G4CutTubs* solidPipe = new G4CutTubs( "solid_pipe", pRMin, pRMax, pDz, pSPhi, pDPhi, pLowNorm, pHighNorm );
  
  // G4LogicalVolume* logicPipe = new G4LogicalVolume(solidPipe, H2O, "logic_pipe");


  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////FONTE DE IRIDIO 192////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  ///// Feita por Sane Simone /////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////


  /////////////////////////
  // Encapsulamento 
  /////////////////////////


  G4double raio_int = 0.*mm; //Dimensoes
  G4double raio_ext = 0.80*mm;
  G4double altura = 2.35*mm;//(4.5mm/2)
  G4double angulo_inicial = 0*rad;
  G4double angulo_final = (2*pi)*rad;

  G4Tubs* blind = new G4Tubs("blind", raio_int, raio_ext, altura, angulo_inicial, angulo_final); //Formato

  G4LogicalVolume* blind_log = new G4LogicalVolume (blind, steelAISI321,"blind_log", 0, 0, 0); //Composicao

  G4PVPlacement* blind_fis = new G4PVPlacement (0,G4ThreeVector (posIni.x(),posIni.y(),posIni.z()+1.40*mm), blind_log,"blind",logicWorld, false, 0); //Posicao


  ////////////////////
  // Cable
  ///////////////////


  //---------
  // Interno
  //---------

  G4double raio_intC = 0.*mm; //Dimensoes
  G4double raio_extC = 0.55*mm; //(1.1mm de diametro)
  G4double alturaC = 0.65*mm;//(1.3mm/2)

  G4Tubs* cabo = new G4Tubs("cabo", raio_intC, raio_extC, alturaC, angulo_inicial, angulo_final);//Formato

  G4LogicalVolume* cabo_log = new G4LogicalVolume(cabo, steelAISI301,"cabo_log", 0, 0, 0);//Composicao

  G4PVPlacement* cabo_fis = new G4PVPlacement(0,G4ThreeVector(0,0,3.1*mm), cabo_log, "cabo",logicWorld , false, 0);//Posicao



  //----------
  // Externo
  //---------


  G4double raio_intCE = 0.*mm; //Dimensoes
  G4double raio_extCE = 0.55*mm;
  G4double alturaCE = 30.0*mm;//(1.3mm/2)

  G4Tubs* caboE = new G4Tubs("caboE", raio_intCE, raio_extCE, alturaCE, angulo_inicial, angulo_final);//Formato

  G4LogicalVolume* caboE_log = new G4LogicalVolume(caboE, steelAISI301,"caboE_log", 0, 0, 0);//Composicao

  G4PVPlacement* caboE_fis = new G4PVPlacement(0,G4ThreeVector(posIni.x(),posIni.y(),posIni.z()+33.75*mm), caboE_log, "caboE", logicWorld, false, 0);//Posicao


  //-------------
  // Cavidade
  //-------------

  G4double raio_intcavit = 0*mm; //Dimensoes
  G4double raio_extcavit = 0.6*mm;
  G4double alturacavit = 0.70*mm;//(3.6mm/2)
  G4double ang_inicialcavit = 0*rad;
  G4double ang_finalcavit = (2*pi)*rad;

  G4Tubs* cavit = new G4Tubs("cavit", raio_intcavit, raio_extcavit,alturacavit,ang_inicialcavit, ang_finalcavit); //Formato
  G4LogicalVolume* cavit_log = new G4LogicalVolume(cavit, Air, "cavit_log",0,0,0); //Composicao
  G4PVPlacement* cavit_fis = new G4PVPlacement(0,G4ThreeVector(0,0,-1.45*mm), cavit_log, "cavit",blind_log , true, 0); //Posicao


  //-------------
  // Ir-192 core
  //-------------


  G4double raio_intcore = 0*mm; //Dimensoes
  G4double raio_extcore = 0.5*mm;
  G4double alturacore = 0.65*mm;
  G4double ang_inicialcore = 0*rad;
  G4double ang_finalcore = (2*pi)*rad;

  G4Tubs* font = new G4Tubs("font", raio_intcore, raio_extcore,alturacore,ang_inicialcore, ang_finalcore); //Formato
  G4LogicalVolume* font_log = new G4LogicalVolume(font, Ir192, "font_log",0,0,0); //Composicao
  G4PVPlacement* font_fis = new G4PVPlacement(0,G4ThreeVector(0,0,0.05*mm), font_log, "font",cavit_log, true, 0); //Posicao


  //-------------
  // Semiesfera
  //-------------


  G4double pRmin = 0*mm;
  G4double pRmax = 0.8*mm;
  G4double pSPhi = 0*rad;
  G4double pDPhi = (2*pi)*rad;
  G4double pSTheta = (0.5*pi)*rad;
  G4double pDTheta = (pi)*rad;


  G4Sphere* sesf = new G4Sphere("sesf", pRmin, pRmax, pSPhi, pDPhi, pSTheta, pDTheta);    
  G4LogicalVolume* sesf_log = new G4LogicalVolume(sesf, steelAISI321, "sesf_log",0,0,0);
  G4PVPlacement* sesf_fis = new G4PVPlacement(0,G4ThreeVector(posIni.x(),posIni.y(),posIni.z()-0.95*mm), sesf_log, "sesf", logicWorld, false, 0); 

  // Primeiro anel
  G4double alfa1 = 180*deg;
  G4RotationMatrix* xRot1 = new G4RotationMatrix;
  xRot1->rotateX(alfa1); 
  G4double pRmin1 = 1.00*cm;
  G4double pRmax1 = 1.10*cm;
  pSPhi = 0*rad;
  pDPhi = (2*pi)*rad;
  pSTheta = (50)*deg;
  pDTheta = (60)*deg;

  G4Sphere* esf1 = new G4Sphere("esf1", pRmin1, pRmax1, pSPhi, pDPhi, pSTheta, pDTheta);
  G4LogicalVolume* esf_log1 = new G4LogicalVolume(esf1, H2O, "esf_log1", 0, 0, 0);
  G4PVPlacement* esf_fis1 = new G4PVPlacement(xRot1, G4ThreeVector(0, 0, 0.0*mm), esf_log1, "esf1", logicWorld, false, 0);

  // Segundo anel
  G4double alfa2 = 180*deg;
  G4RotationMatrix* xRot2 = new G4RotationMatrix;
  xRot2->rotateX(alfa2); 
  G4double pRmin2 = 2.00*cm;  // Diferente do primeiro anel
  G4double pRmax2 = 2.10*cm;

  G4Sphere* esf2 = new G4Sphere("esf2", pRmin2, pRmax2, pSPhi, pDPhi, pSTheta, pDTheta);
  G4LogicalVolume* esf_log2 = new G4LogicalVolume(esf2, H2O, "esf_log2", 0, 0, 0);
  G4PVPlacement* esf_fis2 = new G4PVPlacement(xRot2, G4ThreeVector(0, 0, 0.0*mm), esf_log2, "esf2", logicWorld, false, 0);

  // Definindo volumes de pontuação
  fScoringVolume1 = esf_log1;
  fScoringVolume2 = esf_log2;


  // //-------------
  // //Mudando a cor
  // //-------------


  G4VisAttributes* blindVisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    blind_log->SetVisAttributes(blindVisAtt);

  G4VisAttributes* sesfAtt= new G4VisAttributes(G4Colour(0.0,0.0,5.0));
    sesf_log->SetVisAttributes(sesfAtt);
    
  G4VisAttributes* caboVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));//white
    cabo_log->SetVisAttributes(caboVisAtt);  

  G4VisAttributes* FontVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    font_log->SetVisAttributes(FontVisAtt);

  G4VisAttributes* esfAtt1 = new G4VisAttributes(G4Colour(3.0,0.0,5.0));
    esf_log1->SetVisAttributes(esfAtt1);

  G4VisAttributes* esfAtt2 = new G4VisAttributes(G4Colour(3.0,0.0,6.0));
    esf_log2->SetVisAttributes(esfAtt2);

  // definir material da fonte como vacuo ou ar e emitir elétrons da fonte e transformar detector em uma esfera onde elétrons deeriam parar. OK
  // depositar toda energia no esfera, sabendo E eletron e qnts elétrons emitidos, sabemos energia depositada na casca esférica. precisamos saber a massa da casca. OK
  // 4pi* rao cubo* /3 * densidade = m
  // 4pi* rao cubo* /3 * 7.9*g/cm3 = m
  // 4pi* (0,31ao cubo-0,3ao cubo)* /3 *7.9 g/cm3 = m
  // 10000 elétrons com E = 1 MeV
  // Cumulated dose per run, in scoring volume : 15.5679352128739 picoGy  rms = 0.01620041806157135 picoGy 

  
  return physWorld;
}
