#ifndef MAGFIELD_UNIFORM_FIELD_FAST_SIMULATION_MODEL_HH
#define MAGFIELD_UNIFORM_FIELD_FAST_SIMULATION_MODEL_HH

#include "G4VFastSimulationModel.hh"

class G4FastTrack;
class G4Material;
class G4Region;
class G4Track;
class UniformFieldFastGeometry;

class UniformFieldFastSimulationModel final : public G4VFastSimulationModel {
public:
  explicit UniformFieldFastSimulationModel(
      const UniformFieldFastGeometry* fastGeometry);

  // Validate geometry overlap policy, then attach the model if enabled by the
  // macro configuration.
  static void ConfigureAndAttach(
      const UniformFieldFastGeometry* fastGeometry);
  static void ActivateForChargedParticles();

  G4bool IsApplicable(const G4ParticleDefinition& particle) override;
  G4bool ModelTrigger(const G4FastTrack& fastTrack) override;
  void DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) override;

private:
  G4bool GetTargetSurface(const G4FastTrack& fastTrack,
                          G4double& targetZ) const;
  G4bool GetUniformMagneticFieldZ(G4double& magneticFieldZ) const;
  G4bool IsOrbitContainedInEnvelope(const G4FastTrack& fastTrack,
                                    G4double magneticFieldZ) const;

  const G4Material* fRequiredMaterial;
};

#endif
