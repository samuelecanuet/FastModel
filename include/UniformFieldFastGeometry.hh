#ifndef MAGFIELD_UNIFORM_FIELD_FAST_GEOMETRY_HH
#define MAGFIELD_UNIFORM_FIELD_FAST_GEOMETRY_HH

#include "G4String.hh"
#include "G4ThreeVector.hh"

#include <vector>

class G4LogicalVolume;
class G4Material;
class G4Region;
class G4VSolid;
class G4VPhysicalVolume;

// Convenience builder for an axis-aligned fast-simulation box or solid,
// full-phi tube (inner radius zero).
// The region name is <name>_Region, while the solid and logical-volume names
// are <name>_Solid and <name>_Logical.
class UniformFieldFastGeometry final {
public:
  UniformFieldFastGeometry(const G4String& name, G4LogicalVolume* mother,
                           const G4ThreeVector& position,
                           const G4ThreeVector& boxHalfLengths);

  UniformFieldFastGeometry(const G4String& name, G4LogicalVolume* mother,
                           const G4ThreeVector& position,
                           G4double outerRadius, G4double halfLengthZ);

  const G4String& GetRegionName() const { return fRegionName; }
  G4VSolid* GetSolid() const { return fSolid; }
  G4LogicalVolume* GetLogicalVolume() const { return fLogicalVolume; }
  G4VPhysicalVolume* GetPhysicalVolume() const { return fPhysicalVolume; }
  G4Region* GetRegion() const { return fRegion; }

  // Check all fast geometries created on this thread after construction ends.
  static void CheckRegisteredOverlaps();

  // Check an explicit collection of fast-volume placements.
  static void CheckForOverlaps(
      const std::vector<G4VPhysicalVolume*>& fastVolumes);

  // Check every placement in the Geant4 physical-volume store. Overlaps are
  // allowed only when the overlapping placement uses vacuumMaterial.
  static void CheckNonVacuumOverlaps(const G4Material* vacuumMaterial);

private:
  void Build(const G4String& name, G4LogicalVolume* mother,
             const G4ThreeVector& position, G4VSolid* solid);

  G4String fRegionName;
  G4VSolid* fSolid = nullptr;
  G4LogicalVolume* fLogicalVolume = nullptr;
  G4VPhysicalVolume* fPhysicalVolume = nullptr;
  G4Region* fRegion = nullptr;
};

#endif
