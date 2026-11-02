#ifndef MAGFIELD_UNIFORM_FIELD_FAST_MESSENGER_HH
#define MAGFIELD_UNIFORM_FIELD_FAST_MESSENGER_HH

#include "G4UImessenger.hh"

class G4UIcmdWithABool;
class G4UIcommand;
class G4UIdirectory;

class UniformFieldFastMessenger final : public G4UImessenger {
public:
  static UniformFieldFastMessenger* Instance();

  G4bool IsEnabled() const { return fEnabled; }
  void SetNewValue(G4UIcommand* command, G4String value) override;

private:
  UniformFieldFastMessenger();

  G4UIdirectory* fDirectory = nullptr;
  G4UIcmdWithABool* fEnableCommand = nullptr;
  G4bool fEnabled = false;
};

#endif
