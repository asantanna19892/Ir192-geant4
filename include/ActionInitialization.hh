#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization();
  ~ActionInitialization();
  virtual void Build() const;
};

#endif