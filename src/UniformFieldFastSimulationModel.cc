#include "UniformFieldFastSimulationModel.hh"

#include "UniformFieldFastGeometry.hh"
#include "UniformFieldFastMessenger.hh"
#include "UniformFieldFastPropagator.hh"

#include "G4AffineTransform.hh"
#include "G4Box.hh"
#include "G4Exception.hh"
#include "G4FastSimulationHelper.hh"
#include "G4FastStep.hh"
#include "G4FastTrack.hh"
#include "G4GeometryTolerance.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4Region.hh"
#include "G4Tubs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4UniformMagField.hh"

#include <cmath>
#include <limits>
#include <utility>

namespace {
constexpr G4double kSurfaceClearance = 0.01 * mm;

G4bool IsSupportedEnvelopeSolid(const G4VSolid* solid) {
  if (dynamic_cast<const G4Box*>(solid) != nullptr) {
    return true;
  }

  const auto* tube = dynamic_cast<const G4Tubs*>(solid);
  if (tube == nullptr || tube->GetInnerRadius() != 0.0) {
    return false;
  }
  const auto angularTolerance =
      G4GeometryTolerance::GetInstance()->GetAngularTolerance();
  return tube->GetDeltaPhiAngle() >= twopi - angularTolerance;
}

G4Region* GetRegion(const UniformFieldFastGeometry* fastGeometry) {
  if (fastGeometry == nullptr) {
    G4ExceptionDescription description;
    description << "The UniformFieldFastSimulationModel requires a non-null "
                   "UniformFieldFastGeometry pointer.";
    G4Exception("UniformFieldFastSimulationModel::UniformFieldFastSimulationModel",
                "UniformFieldFastSim004", FatalException, description);
    return nullptr;
  }

  auto* region = fastGeometry->GetRegion();
  if (region == nullptr) {
    G4ExceptionDescription description;
    description << "The supplied UniformFieldFastGeometry has no region.";
    G4Exception("UniformFieldFastSimulationModel::UniformFieldFastSimulationModel",
                "UniformFieldFastSim009", FatalException, description);
  }
  return region;
}
}

UniformFieldFastSimulationModel::UniformFieldFastSimulationModel(
    const UniformFieldFastGeometry* fastGeometry)
    : G4VFastSimulationModel("FastMagneticField", GetRegion(fastGeometry)),
      fRequiredMaterial(
          G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic")) {
  auto* region = fastGeometry == nullptr ? nullptr : fastGeometry->GetRegion();
  if (region == nullptr || fRequiredMaterial == nullptr) {
    G4ExceptionDescription description;
    description << "A region and the G4_Galactic material are required.";
    G4Exception("UniformFieldFastSimulationModel::UniformFieldFastSimulationModel",
                "UniformFieldFastSim001", FatalException, description);
    return;
  }

  if (region->GetNumberOfRootVolumes() == 0) {
    G4ExceptionDescription description;
    description << "The fast-simulation region must contain at least one root "
                   "logical volume.";
    G4Exception("UniformFieldFastSimulationModel::UniformFieldFastSimulationModel",
                "UniformFieldFastSim002", FatalException, description);
    return;
  }

  auto rootVolume = region->GetRootLogicalVolumeIterator();
  for (std::size_t i = 0; i < region->GetNumberOfRootVolumes(); ++i, ++rootVolume) {
    const auto* logicalVolume = *rootVolume;
    if (logicalVolume == nullptr ||
        !IsSupportedEnvelopeSolid(logicalVolume->GetSolid())) {
      G4ExceptionDescription description;
      description << "Fast-simulation envelope root volumes must use G4Box "
                     "or solid, full-phi G4Tubs solids (inner radius zero). "
                     "Unsupported root volume: "
                  << (logicalVolume == nullptr ? "<null>"
                                                : logicalVolume->GetName());
      G4Exception("UniformFieldFastSimulationModel::UniformFieldFastSimulationModel",
                  "UniformFieldFastSim003", FatalException, description);
      return;
    }
  }

  G4double magneticFieldZ = 0.0;
  GetUniformMagneticFieldZ(magneticFieldZ);
}

void UniformFieldFastSimulationModel::ConfigureAndAttach(
    const UniformFieldFastGeometry* fastGeometry) {
  const auto* vacuumMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  UniformFieldFastGeometry::CheckNonVacuumOverlaps(vacuumMaterial);
  UniformFieldFastGeometry::CheckRegisteredOverlaps();

  if (!UniformFieldFastMessenger::Instance()->IsEnabled()) {
    return;
  }

  if (fastGeometry == nullptr) {
    G4ExceptionDescription description;
    description << "Fast mode is enabled, but no UniformFieldFastGeometry was "
                   "provided.";
    G4Exception("UniformFieldFastSimulationModel::ConfigureAndAttach",
                "UniformFieldFastSim010", FatalException, description);
    return;
  }

  new UniformFieldFastSimulationModel(fastGeometry);
}

void UniformFieldFastSimulationModel::ActivateForChargedParticles() {
  if (!UniformFieldFastMessenger::Instance()->IsEnabled()) {
    return;
  }

  auto* particleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    auto* particle = particleIterator->value();
    if (particle->GetPDGCharge() == 0.0) {
      continue;
    }
    auto* processManager = particle->GetProcessManager();
    if (processManager != nullptr) {
      G4FastSimulationHelper::ActivateFastSimulation(processManager);
    }
  }
}

G4bool UniformFieldFastSimulationModel::IsApplicable(
    const G4ParticleDefinition& particle) {
  return std::abs(particle.GetPDGCharge()) > 0.0;
}

G4bool UniformFieldFastSimulationModel::GetUniformMagneticFieldZ(
    G4double& magneticFieldZ) const {
  auto* fieldManager =
      G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if (fieldManager == nullptr) {
    G4ExceptionDescription description;
    description << "The global field manager is not available. "
                   "FastMagneticField requires a global, uniform magnetic "
                   "field aligned with z.";
    G4Exception("UniformFieldFastSimulationModel::GetUniformMagneticFieldZ",
                "UniformFieldFastSim005", FatalException, description);
    return false;
  }

  const auto* detectorField = fieldManager->GetDetectorField();
  if (detectorField == nullptr) {
    G4ExceptionDescription description;
    description << "No detector field is installed on the global field "
                   "manager. FastMagneticField cannot propagate without a "
                   "field.";
    G4Exception("UniformFieldFastSimulationModel::GetUniformMagneticFieldZ",
                "UniformFieldFastSim006", FatalException, description);
    return false;
  }

  const auto* uniformField =
      dynamic_cast<const G4UniformMagField*>(detectorField);
  if (uniformField == nullptr) {
    G4ExceptionDescription description;
    description << "The installed global detector field is not a "
                   "G4UniformMagField. FastMagneticField only supports a "
                   "uniform magnetic field aligned with z.";
    G4Exception("UniformFieldFastSimulationModel::GetUniformMagneticFieldZ",
                "UniformFieldFastSim007", FatalException, description);
    return false;
  }

  const auto fieldValue = uniformField->GetConstantFieldValue();
  constexpr G4double transverseFieldTolerance = 1.0e-12 * tesla;
  if (std::abs(fieldValue.x()) > transverseFieldTolerance ||
      std::abs(fieldValue.y()) > transverseFieldTolerance) {
    G4ExceptionDescription description;
    description << "The installed uniform magnetic field must be aligned "
                   "with z; received B = ("
                << fieldValue.x() / tesla << ", "
                << fieldValue.y() / tesla << ", "
                << fieldValue.z() / tesla << ") tesla.";
    G4Exception("UniformFieldFastSimulationModel::GetUniformMagneticFieldZ",
                "UniformFieldFastSim008", FatalException, description);
    return false;
  }
  magneticFieldZ = fieldValue.z();
  return true;
}

G4bool UniformFieldFastSimulationModel::GetTargetSurface(
    const G4FastTrack& fastTrack, G4double& targetZ) const {
  const auto* track = fastTrack.GetPrimaryTrack();
  const auto* solid = fastTrack.GetEnvelopeSolid();
  const auto* globalToLocal = fastTrack.GetAffineTransformation();
  const auto* localToGlobal = fastTrack.GetInverseAffineTransformation();
  if (track == nullptr || solid == nullptr || globalToLocal == nullptr ||
      localToGlobal == nullptr) {
    return false;
  }

  G4double halfLengthZ = 0.0;
  if (const auto* box = dynamic_cast<const G4Box*>(solid)) {
    halfLengthZ = box->GetZHalfLength();
  } else if (const auto* tube = dynamic_cast<const G4Tubs*>(solid)) {
    halfLengthZ = tube->GetZHalfLength();
  } else {
    return false;
  }

  const auto localGlobalZ =
      globalToLocal->TransformAxis(G4ThreeVector(0.0, 0.0, 1.0)).unit();
  constexpr G4double axisTolerance = 1.0e-10;
  if (std::abs(localGlobalZ.x()) > axisTolerance ||
      std::abs(localGlobalZ.y()) > axisTolerance) {
    return false;
  }

  const G4double pz = track->GetMomentum().z();
  if (std::abs(pz) < 1.0e-15 * MeV) {
    return false;
  }

  const G4double directionZ = pz > 0.0 ? 1.0 : -1.0;
  const G4double currentZ = track->GetPosition().z();
  const G4double lowerZ = localToGlobal
                              ->TransformPoint(G4ThreeVector(0.0, 0.0,
                                                            -halfLengthZ))
                              .z();
  const G4double upperZ = localToGlobal
                              ->TransformPoint(G4ThreeVector(0.0, 0.0,
                                                            halfLengthZ))
                              .z();
  const G4double surfaces[] = {lowerZ, upperZ};
  G4double nearestForwardDistance = std::numeric_limits<G4double>::max();
  G4bool foundSurface = false;

  for (const auto surfaceZ : surfaces) {
    const G4double candidateZ =
        surfaceZ - directionZ * kSurfaceClearance;
    const G4double deltaZ = candidateZ - currentZ;
    const G4double distanceAlongTrack = deltaZ / pz;
    if (distanceAlongTrack <= 0.0) {
      continue;
    }
    if (std::abs(deltaZ) < nearestForwardDistance) {
      nearestForwardDistance = std::abs(deltaZ);
      targetZ = candidateZ;
      foundSurface = true;
    }
  }
  return foundSurface;
}

G4bool UniformFieldFastSimulationModel::ModelTrigger(
    const G4FastTrack& fastTrack) {
  const auto* track = fastTrack.GetPrimaryTrack();
  if (track == nullptr ||
      track->GetMaterial() != fRequiredMaterial) {
    return false;
  }

  G4double magneticFieldZ = 0.0;
  if (!GetUniformMagneticFieldZ(magneticFieldZ)) {
    return false;
  }

  // IsApplicable only receives the particle definition. State-dependent
  // checks belong here, where the initial position and momentum are available.
  if (!IsOrbitContainedInEnvelope(fastTrack, magneticFieldZ)) {
    return false;
  }

  G4double targetZ = 0.0;
  if (!GetTargetSurface(fastTrack, targetZ)) {
    return false;
  }

  UniformFieldTrajectoryState targetState;
  return UniformFieldFastPropagator::PropagateToZ(
      track->GetPosition(), track->GetMomentum(),
      track->GetDefinition()->GetPDGCharge(),
      track->GetDefinition()->GetPDGMass(), track->GetKineticEnergy(),
      magneticFieldZ, targetZ, targetState);
}

G4bool UniformFieldFastSimulationModel::IsOrbitContainedInEnvelope(
    const G4FastTrack& fastTrack, G4double magneticFieldZ) const {
  const auto* track = fastTrack.GetPrimaryTrack();
  const auto* solid = fastTrack.GetEnvelopeSolid();
  const auto* box = dynamic_cast<const G4Box*>(solid);
  const auto* tube = dynamic_cast<const G4Tubs*>(solid);
  const auto* transform = fastTrack.GetAffineTransformation();
  if (track == nullptr || (!box && !tube) || transform == nullptr) {
    return false;
  }

  const auto localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  const auto localMomentum = fastTrack.GetPrimaryTrackLocalMomentum();
  // A track entering the fast region can start its first step exactly on the
  // envelope surface. Accept that state; GetTargetSurface below still requires
  // a forward z face, so a track leaving the envelope will not trigger here.
  if (solid->Inside(localPosition) == kOutside) {
    return false;
  }

  // The analytic propagator assumes B is parallel to global z. It can be
  // bounded in the envelope frame only when that axis is also a local volume axis.
  const auto localGlobalZ =
      transform->TransformAxis(G4ThreeVector(0.0, 0.0, 1.0)).unit();
  constexpr G4double axisTolerance = 1.0e-10;
  if (std::abs(localGlobalZ.x()) > axisTolerance ||
      std::abs(localGlobalZ.y()) > axisTolerance) {
    return false;
  }

  const auto momentum = localMomentum.mag();
  if (momentum <= 0.0) {
    return false;
  }

  const auto localFieldZ =
      magneticFieldZ * localGlobalZ.z();
  const auto charge = track->GetDefinition()->GetPDGCharge();

  G4double targetZ = 0.0;
  UniformFieldTrajectoryState targetState;
    if (!GetTargetSurface(fastTrack, targetZ) ||
      !UniformFieldFastPropagator::PropagateToZ(
          track->GetPosition(), track->GetMomentum(), charge,
          track->GetDefinition()->GetPDGMass(), track->GetKineticEnergy(),
          magneticFieldZ, targetZ, targetState) ||
      solid->Inside(transform->TransformPoint(targetState.position)) ==
          kOutside) {
    return false;
  }

  if (std::abs(localFieldZ) < 1.0e-15 * tesla || charge == 0.0) {
    // With no curvature, the accepted z-directed segment is straight.
    return true;
  }

  const auto angularRate =
      0.299792458 * (charge / eplus) * (localFieldZ / tesla) /
      (momentum / GeV) / m;
  if (angularRate == 0.0) {
    return false;
  }

  // In a uniform Bz field the transverse path is a circle. Bound the entire
  // orbit using its center and Larmor radius, not just the fast-step endpoint.
  const auto centerX =
      localPosition.x() + localMomentum.y() / (momentum * angularRate);
  const auto centerY =
      localPosition.y() - localMomentum.x() / (momentum * angularRate);
  const auto transverseMomentum =
      std::hypot(localMomentum.x(), localMomentum.y());
  const auto larmorRadius =
      transverseMomentum / (momentum * std::abs(angularRate));
  const auto margin = kSurfaceClearance +
                      G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  if (box != nullptr) {
    const auto availableHalfX = box->GetXHalfLength() - margin;
    const auto availableHalfY = box->GetYHalfLength() - margin;
    return availableHalfX > 0.0 && availableHalfY > 0.0 &&
           std::abs(centerX) + larmorRadius <= availableHalfX &&
           std::abs(centerY) + larmorRadius <= availableHalfY;
  }

  // A solid, full-phi G4Tubs volume is rotationally symmetric about z. Require
  // the complete orbit to fit within its outer radius.
  const auto angularTolerance =
      G4GeometryTolerance::GetInstance()->GetAngularTolerance();
  if (tube->GetDeltaPhiAngle() < twopi - angularTolerance) {
    return false;
  }
  const auto orbitCenterRadius = std::hypot(centerX, centerY);
  const auto maximumOrbitRadius = orbitCenterRadius + larmorRadius;
  const auto maximumTubeRadius = tube->GetOuterRadius() - margin;
  return maximumTubeRadius > 0.0 && maximumOrbitRadius <= maximumTubeRadius;
}

void UniformFieldFastSimulationModel::DoIt(const G4FastTrack& fastTrack,
                                           G4FastStep& fastStep) {
  const auto* track = fastTrack.GetPrimaryTrack();
  if (track == nullptr) {
    return;
  }

  G4double magneticFieldZ = 0.0;
  if (!GetUniformMagneticFieldZ(magneticFieldZ)) {
    return;
  }

  G4double targetZ = 0.0;
  UniformFieldTrajectoryState finalState;
  if (!GetTargetSurface(fastTrack, targetZ) ||
      !UniformFieldFastPropagator::PropagateToZ(
          track->GetPosition(), track->GetMomentum(),
          track->GetDefinition()->GetPDGCharge(),
          track->GetDefinition()->GetPDGMass(), track->GetKineticEnergy(),
          magneticFieldZ, targetZ, finalState)) {
    return;
  }

  fastStep.ProposePrimaryTrackFinalPosition(finalState.position, false);
  fastStep.ProposePrimaryTrackFinalMomentumDirection(
      finalState.momentum.unit(), false);
  fastStep.ProposePrimaryTrackFinalKineticEnergy(track->GetKineticEnergy());
  fastStep.ProposePrimaryTrackFinalTime(track->GetGlobalTime() +
                                        finalState.timeOfFlight);
  fastStep.ProposePrimaryTrackFinalProperTime(track->GetProperTime() +
                                              finalState.properTimeOfFlight);
  fastStep.ProposePrimaryTrackPathLength(finalState.pathLength);
}
