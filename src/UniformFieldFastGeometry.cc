#include "UniformFieldFastGeometry.hh"

#include "G4Box.hh"
#include "G4AffineTransform.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysicalConstants.hh"
#include "G4Region.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace {
std::vector<G4AffineTransform> GetWorldTransforms(
    const G4VPhysicalVolume* volume, std::size_t depth = 0) {
  const auto* volumeStore = G4PhysicalVolumeStore::GetInstance();
  if (volume == nullptr || depth > volumeStore->size()) {
    return {};
  }

  const G4AffineTransform localToMother(volume->GetObjectRotation(),
                                        volume->GetObjectTranslation());
  const auto* motherLogical = volume->GetMotherLogical();
  if (motherLogical == nullptr) {
    return {localToMother};
  }

  std::vector<G4AffineTransform> worldTransforms;
  for (const auto* possibleMother : *volumeStore) {
    if (possibleMother->GetLogicalVolume() != motherLogical) {
      continue;
    }
    for (const auto& motherTransform :
         GetWorldTransforms(possibleMother, depth + 1)) {
      worldTransforms.push_back(motherTransform * localToMother);
    }
  }
  return worldTransforms;
}

G4bool DoFastVolumesOverlap(const G4VPhysicalVolume* first,
                            const G4VPhysicalVolume* second) {
  const auto* firstSolid = first->GetLogicalVolume()->GetSolid();
  const auto* secondSolid = second->GetLogicalVolume()->GetSolid();
  const auto firstTransforms = GetWorldTransforms(first);
  const auto secondTransforms = GetWorldTransforms(second);
  if (firstTransforms.empty() || secondTransforms.empty()) {
    G4ExceptionDescription description;
    description << "Could not resolve the placement hierarchy while checking "
                   "fast-volume overlap between '"
                << first->GetName() << "' and '" << second->GetName() << "'.";
    G4Exception("UniformFieldFastGeometry::CheckForOverlaps",
                "UniformFieldFastGeometry004", FatalException, description);
    return false;
  }

  constexpr G4int samplesPerDirection = 10000;
  for (const auto& firstTransform : firstTransforms) {
    for (const auto& secondTransform : secondTransforms) {
      for (G4int sample = 0; sample < samplesPerDirection; ++sample) {
        const auto firstSurfacePoint =
            firstTransform.TransformPoint(firstSolid->GetPointOnSurface());
        const auto firstPointInSecond =
            secondTransform.InverseTransformPoint(firstSurfacePoint);
        if (secondSolid->Inside(firstPointInSecond) == kInside) {
          return true;
        }

        const auto secondSurfacePoint =
            secondTransform.TransformPoint(secondSolid->GetPointOnSurface());
        const auto secondPointInFirst =
            firstTransform.InverseTransformPoint(secondSurfacePoint);
        if (firstSolid->Inside(secondPointInFirst) == kInside) {
          return true;
        }
      }
    }
  }
  return false;
}

std::vector<G4VPhysicalVolume*>& RegisteredFastVolumes() {
  static std::vector<G4VPhysicalVolume*> fastVolumes;
  return fastVolumes;
}

void RegisterFastVolume(G4VPhysicalVolume* fastVolume) {
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
  fastVolumes.push_back(fastVolume);
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
  fRegion = new G4Region(fRegionName);
  fRegion->AddRootLogicalVolume(fLogicalVolume);
  RegisterFastVolume(fPhysicalVolume);
}

void UniformFieldFastGeometry::CheckForOverlaps(
    const std::vector<G4VPhysicalVolume*>& fastVolumes) {
  for (std::size_t firstIndex = 0; firstIndex < fastVolumes.size();
       ++firstIndex) {
    for (std::size_t secondIndex = firstIndex + 1;
         secondIndex < fastVolumes.size(); ++secondIndex) {
      const auto* first = fastVolumes[firstIndex];
      const auto* second = fastVolumes[secondIndex];
      if (first == nullptr || second == nullptr ||
          !DoFastVolumesOverlap(first, second)) {
        continue;
      }

      G4ExceptionDescription description;
      description << "Fast-simulation volumes '" << first->GetName()
                  << "' and '" << second->GetName()
                  << "' overlap. Fast volumes must not overlap, including "
                     "when both use vacuum material.";
      G4Exception("UniformFieldFastGeometry::CheckForOverlaps",
                  "UniformFieldFastGeometry005", FatalException, description);
    }
  }
}

void UniformFieldFastGeometry::CheckRegisteredOverlaps() {
  CheckForOverlaps(RegisteredFastVolumes());
}

void UniformFieldFastGeometry::CheckNonVacuumOverlaps(
    const G4Material* vacuumMaterial) {
  if (vacuumMaterial == nullptr) {
    G4Exception("UniformFieldFastGeometry::CheckNonVacuumOverlaps",
                "UniformFieldFastGeometry006", FatalException,
                "A vacuum material is required to check geometry overlaps.");
    return;
  }

  for (auto* volume : *G4PhysicalVolumeStore::GetInstance()) {
    if (volume == nullptr || !volume->CheckOverlaps(10000, 0.0, false)) {
      continue;
    }

    const auto* logicalVolume = volume->GetLogicalVolume();
    const auto* material =
        logicalVolume == nullptr ? nullptr : logicalVolume->GetMaterial();
    if (material == vacuumMaterial) {
      continue;
    }

    G4ExceptionDescription description;
    description << "Geometry overlap detected for physical volume '"
                << volume->GetName() << "' with material '"
                << (material == nullptr ? "<null>" : material->GetName())
                << "'. Overlaps are allowed only for vacuum volumes.";
    G4Exception("UniformFieldFastGeometry::CheckNonVacuumOverlaps",
                "UniformFieldFastGeometry007", FatalException, description);
  }
}
