#ifndef MAGFIELD_UNIFORM_FIELD_FAST_PROPAGATOR_HH
#define MAGFIELD_UNIFORM_FIELD_FAST_PROPAGATOR_HH

#include "G4ThreeVector.hh"

struct UniformFieldTrajectoryState {
  G4ThreeVector position;
  G4ThreeVector momentum;
  G4double pathLength = 0.0;
  G4double timeOfFlight = 0.0;
  G4double properTimeOfFlight = 0.0;
};

class UniformFieldFastPropagator final {
public:
  // Propagate an initial state exactly through a uniform Bz field to a z plane.
  // Position is returned in Geant4 length units and momentum in Geant4 units.
  static G4bool PropagateToZ(const G4ThreeVector& initialPosition,
                            const G4ThreeVector& initialMomentum,
                            G4double charge, G4double mass,
                            G4double kineticEnergy, G4double magneticFieldZ,
                            G4double targetZ,
                            UniformFieldTrajectoryState& finalState);
};

#endif
