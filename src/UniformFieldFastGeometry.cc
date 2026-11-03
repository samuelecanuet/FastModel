#include "UniformFieldFastGeometry.hh"

#include "G4Box.hh"
#include "G4AffineTransform.hh"
#include "G4Exception.hh"
#include "G4GeometryTolerance.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Navigator.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysicalConstants.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4Tubs.hh"
#include "G4TransportationManager.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"

#include <algorithm>
#include <cstddef>
#include <mutex>
#include <unordered_set>
#include <vector>

namespace {
void CollectWorldTransforms(
    const G4VPhysicalVolume* currentVolume,
    const G4VPhysicalVolume* targetVolume,
    const G4AffineTransform& currentToWorld,
    std::vector<G4AffineTransform>& worldTransforms) {
  if (currentVolume == nullptr) {
    return;
  }

  if (currentVolume == targetVolume) {
    worldTransforms.push_back(currentToWorld);
  }

  const auto* logicalVolume = currentVolume->GetLogicalVolume();
  if (logicalVolume == nullptr) {
    return;
  }

  for (G4int daughterIndex = 0;
       daughterIndex < logicalVolume->GetNoDaughters(); ++daughterIndex) {
    const auto* daughter = logicalVolume->GetDaughter(daughterIndex);
    if (daughter == nullptr) {
      continue;
    }

    const G4AffineTransform daughterToMother(
        daughter->GetObjectRotation(), daughter->GetObjectTranslation());
    CollectWorldTransforms(daughter, targetVolume,
                           currentToWorld * daughterToMother,
                           worldTransforms);
  }
}

std::vector<G4AffineTransform> GetWorldTransforms(
    const G4VPhysicalVolume* volume,
    const G4VPhysicalVolume* worldVolume) {
  if (volume == nullptr || worldVolume == nullptr) {
    return {};
  }

  const G4AffineTransform worldToGlobal(worldVolume->GetObjectRotation(),
                                        worldVolume->GetObjectTranslation());
  std::vector<G4AffineTransform> worldTransforms;
  CollectWorldTransforms(worldVolume, volume, worldToGlobal, worldTransforms);
  return worldTransforms;
}

struct FastVolumeInstance {
  const G4VPhysicalVolume* volume;
  G4AffineTransform localToWorld;
};

G4bool HaveSameTransform(const G4AffineTransform& first,
                         const G4AffineTransform& second) {
  const G4ThreeVector probes[] = {
      G4ThreeVector(), G4ThreeVector(1.0 * mm, 0.0, 0.0),
      G4ThreeVector(0.0, 1.0 * mm, 0.0),
      G4ThreeVector(0.0, 0.0, 1.0 * mm)};
  const auto tolerance =
      G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  for (const auto& probe : probes) {
    if ((first.TransformPoint(probe) - second.TransformPoint(probe)).mag() >
        tolerance) {
      return false;
    }
  }
  return true;
}

G4bool DoFastVolumeInstancesOverlap(const FastVolumeInstance& first,
                                    const FastVolumeInstance& second) {
  const auto* firstSolid = first.volume->GetLogicalVolume()->GetSolid();
  const auto* secondSolid = second.volume->GetLogicalVolume()->GetSolid();
  if (HaveSameTransform(first.localToWorld, second.localToWorld) &&
      firstSolid->Inside(G4ThreeVector()) == kInside &&
      secondSolid->Inside(G4ThreeVector()) == kInside) {
    return true;
  }

  constexpr G4int samplesPerDirection = 10000;
  for (G4int sample = 0; sample < samplesPerDirection; ++sample) {
    const auto firstSurfacePoint = first.localToWorld.TransformPoint(
        firstSolid->GetPointOnSurface());
    const auto firstPointInSecond = second.localToWorld.InverseTransformPoint(
        firstSurfacePoint);
    if (secondSolid->Inside(firstPointInSecond) == kInside) {
      return true;
    }

    const auto secondSurfacePoint = second.localToWorld.TransformPoint(
        secondSolid->GetPointOnSurface());
    const auto secondPointInFirst = first.localToWorld.InverseTransformPoint(
        secondSurfacePoint);
    if (firstSolid->Inside(secondPointInFirst) == kInside) {
      return true;
    }
  }
  return false;
}

std::vector<G4VPhysicalVolume*>& RegisteredFastVolumes() {
  static std::vector<G4VPhysicalVolume*> fastVolumes;
  return fastVolumes;
}

std::mutex& FastVolumeRegistryMutex() {
  static std::mutex mutex;
  return mutex;
}

void RegisterFastVolume(G4VPhysicalVolume* fastVolume) {
  const std::lock_guard<std::mutex> lock(FastVolumeRegistryMutex());
  auto& fastVolumes = RegisteredFastVolumes();
  const auto* volumeStore = G4PhysicalVolumeStore::GetInstance();
  fastVolumes.erase(
      std::remove_if(fastVolumes.begin(), fastVolumes.end(),
                     [volumeStore](const G4VPhysicalVolume* volume) {
                       return std::find(volumeStore->begin(),
                                        volumeStore->end(), volume) ==
                              volumeStore->end();
                     }),
      fastVolumes.end());
  if (std::find(fastVolumes.begin(), fastVolumes.end(), fastVolume) ==
      fastVolumes.end()) {
    fastVolumes.push_back(fastVolume);
  }
}
}

UniformFieldFastGeometry::UniformFieldFastGeometry(
    const G4String& name, G4LogicalVolume* mother,
    const G4ThreeVector& position, const G4ThreeVector& boxHalfLengths) {
  if (boxHalfLengths.x() <= 0.0 || boxHalfLengths.y() <= 0.0 ||
      boxHalfLengths.z() <= 0.0) {
    G4Exception("UniformFieldFastGeometry::UniformFieldFastGeometry",
                "UniformFieldFastGeometry001", FatalException,
                "Box half-lengths must be positive.");
    return;
  }

  Build(name, mother, position,
        new G4Box(name + "_Solid", boxHalfLengths.x(), boxHalfLengths.y(),
                  boxHalfLengths.z()));
}

UniformFieldFastGeometry::UniformFieldFastGeometry(
    const G4String& name, G4LogicalVolume* mother,
    const G4ThreeVector& position, G4double outerRadius,
    G4double halfLengthZ) {
  if (outerRadius <= 0.0 || halfLengthZ <= 0.0) {
    G4Exception("UniformFieldFastGeometry::UniformFieldFastGeometry",
                "UniformFieldFastGeometry002", FatalException,
                "A positive outer radius and z half-length are required for "
                "a solid fast tube.");
    return;
  }

  Build(name, mother, position,
        new G4Tubs(name + "_Solid", 0.0, outerRadius, halfLengthZ, 0.0,
                   twopi));
}

void UniformFieldFastGeometry::Build(const G4String& name,
                                     G4LogicalVolume* mother,
                                     const G4ThreeVector& position,
                                     G4VSolid* solid) {
  if (name.empty() || mother == nullptr || solid == nullptr) {
    G4ExceptionDescription description;
    description << "A non-empty name, mother logical volume, and solid are "
                   "required to build fast geometry.";
    G4Exception("UniformFieldFastGeometry::Build",
                "UniformFieldFastGeometry003", FatalException, description);
    return;
  }

  auto* vacuum =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  if (vacuum == nullptr) {
    G4Exception("UniformFieldFastGeometry::Build",
                "UniformFieldFastGeometry008", FatalException,
                "Could not find or build the required G4_Galactic material.");
    return;
  }

  fRegionName = name + "_Region";
  fSolid = solid;
  fLogicalVolume = new G4LogicalVolume(fSolid, vacuum, name + "_Logical");
  fPhysicalVolume = new G4PVPlacement(nullptr, position, fLogicalVolume, name,
                                      mother, false, 0, false);
  fRegion = G4RegionStore::GetInstance()->FindOrCreateRegion(fRegionName);
  fRegion->AddRootLogicalVolume(fLogicalVolume);
  RegisterFastVolume(fPhysicalVolume);
}

void UniformFieldFastGeometry::CheckForOverlaps(
    const std::vector<G4VPhysicalVolume*>& fastVolumes) {
  auto* transportationManager =
      G4TransportationManager::GetTransportationManager();
  auto* navigator = transportationManager == nullptr
                        ? nullptr
                        : transportationManager->GetNavigatorForTracking();
  const auto* worldVolume =
      navigator == nullptr ? nullptr : navigator->GetWorldVolume();
  if (worldVolume == nullptr) {
    G4Exception("UniformFieldFastGeometry::CheckForOverlaps",
                "UniformFieldFastGeometry004", FatalException,
                "Cannot check fast-volume overlaps before the tracking world "
                "is available.");
    return;
  }

  // A geometry reinitialization can leave old placements in Geant4's stores.
  // Keep only registry entries reachable from the current world, then make
  // one record per active placement instance (a logical mother can be reused).
  std::vector<FastVolumeInstance> activeInstances;
  std::unordered_set<const G4VPhysicalVolume*> seenVolumes;
  for (const auto* volume : fastVolumes) {
    if (volume == nullptr || !seenVolumes.insert(volume).second) {
      continue;
    }

    for (const auto& transform : GetWorldTransforms(volume, worldVolume)) {
      activeInstances.push_back({volume, transform});
    }
  }

  for (std::size_t firstIndex = 0; firstIndex < activeInstances.size();
       ++firstIndex) {
    for (std::size_t secondIndex = firstIndex + 1;
         secondIndex < activeInstances.size(); ++secondIndex) {
      const auto& first = activeInstances[firstIndex];
      const auto& second = activeInstances[secondIndex];
      if (!DoFastVolumeInstancesOverlap(first, second)) {
        continue;
      }

      G4ExceptionDescription description;
      description << "Fast-simulation volumes '" << first.volume->GetName()
                  << "' and '" << second.volume->GetName()
                  << "' overlap. Fast volumes must not overlap, including "
                     "when both use vacuum material.";
      G4Exception("UniformFieldFastGeometry::CheckForOverlaps",
                  "UniformFieldFastGeometry005", FatalException, description);
    }
  }
}

void UniformFieldFastGeometry::CheckRegisteredOverlaps() {
  // ConstructSDandField is called for each worker in MT runs. The master owns
  // the shared geometry check; worker navigators and stores must not race over
  // the process-wide registry or repeat its randomized surface sampling.
  if (G4Threading::IsWorkerThread()) {
    return;
  }

  std::vector<G4VPhysicalVolume*> fastVolumes;
  {
    const std::lock_guard<std::mutex> lock(FastVolumeRegistryMutex());
    fastVolumes = RegisteredFastVolumes();
  }
  CheckForOverlaps(fastVolumes);
}
