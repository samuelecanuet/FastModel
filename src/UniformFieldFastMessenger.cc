#include "UniformFieldFastMessenger.hh"

#include "G4ApplicationState.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"

UniformFieldFastMessenger* UniformFieldFastMessenger::Instance() {
  static auto* configuration = new UniformFieldFastMessenger;
  return configuration;
}

UniformFieldFastMessenger::UniformFieldFastMessenger() {
  fDirectory = new G4UIdirectory("/fastModel/");
  fDirectory->SetGuidance("Fast simulation model configuration.");

  fEnableCommand = new G4UIcmdWithABool("/fastModel/enable", this);
  fEnableCommand->SetGuidance(
      "Enable or disable the uniform-field fast simulation model.");
  fEnableCommand->SetParameterName("enabled", false);
  fEnableCommand->SetDefaultValue(false);
  fEnableCommand->AvailableForStates(G4State_PreInit);
}

void UniformFieldFastMessenger::SetNewValue(G4UIcommand* command,
                                            G4String value) {
  if (command == fEnableCommand) {
    fEnabled = fEnableCommand->GetNewBoolValue(value);
  }
}
