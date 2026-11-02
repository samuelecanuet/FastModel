#include "UniformFieldFastPropagator.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

G4bool UniformFieldFastPropagator::PropagateToZ(
    const G4ThreeVector& initialPosition,
    const G4ThreeVector& initialMomentum, G4double charge, G4double mass,
    G4double kineticEnergy, G4double magneticFieldZ, G4double targetZ,
    UniformFieldTrajectoryState& finalState) {
  const G4double pz = initialMomentum.z();
  const G4double momentumMagnitude = initialMomentum.mag();
  const G4double deltaZ = targetZ - initialPosition.z();

  if (momentumMagnitude == 0.0 || std::abs(pz) < 1.0e-15 * MeV ||
      deltaZ / pz <= 0.0) {
    return false;
  }

  const G4double pzGeV = pz / GeV;
  const G4double pxGeV = initialMomentum.x() / GeV;
  const G4double pyGeV = initialMomentum.y() / GeV;
  const G4double curvaturePerLength =
      0.299792458 * (charge / eplus) * (magneticFieldZ / tesla) /
      pzGeV / m;
  const G4double angle = curvaturePerLength * deltaZ;

  G4double deltaX = 0.0;
  G4double deltaY = 0.0;
  G4double sinAngle = 0.0;
  G4double cosAngle = 1.0;
  if (std::abs(angle) < 1.0e-6) {
    deltaX = initialMomentum.x() / pz * deltaZ +
             0.5 * initialMomentum.y() / pz * curvaturePerLength *
                 deltaZ * deltaZ;
    deltaY = initialMomentum.y() / pz * deltaZ -
             0.5 * initialMomentum.x() / pz * curvaturePerLength *
                 deltaZ * deltaZ;
    sinAngle = angle;
    cosAngle = 1.0 - 0.5 * angle * angle;
  } else if (curvaturePerLength != 0.0) {
    sinAngle = std::sin(angle);
    cosAngle = std::cos(angle);
    const G4double oneMinusCos =
        2.0 * std::sin(0.5 * angle) * std::sin(0.5 * angle);
    const G4double denominator = pzGeV * curvaturePerLength;
    deltaX = (pxGeV * sinAngle + pyGeV * oneMinusCos) / denominator;
    deltaY = (pyGeV * sinAngle - pxGeV * oneMinusCos) / denominator;
  } else {
    deltaX = initialMomentum.x() / pz * deltaZ;
    deltaY = initialMomentum.y() / pz * deltaZ;
  }

  const G4double finalPx = initialMomentum.x() * cosAngle +
                           initialMomentum.y() * sinAngle;
  const G4double finalPy = initialMomentum.y() * cosAngle -
                           initialMomentum.x() * sinAngle;
  finalState.position = G4ThreeVector(initialPosition.x() + deltaX,
                                      initialPosition.y() + deltaY, targetZ);
  finalState.momentum = G4ThreeVector(finalPx, finalPy, pz);
  finalState.pathLength = deltaZ * momentumMagnitude / pz;
  finalState.timeOfFlight =
      finalState.pathLength * (kineticEnergy + mass) /
      (momentumMagnitude * c_light);
  finalState.properTimeOfFlight =
      finalState.pathLength * mass / (momentumMagnitude * c_light);
  return true;
}
