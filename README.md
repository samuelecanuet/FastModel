# Uniform-field fast model

This directory contains a reusable Geant4 component for analytically
propagating charged tracks through a vacuum fast-simulation volume in a
uniform magnetic field along z. It is built as the `FastModel::UniformField`
CMake target.

## What it does

`UniformFieldFastGeometry` creates a solid, logical volume, placement, and
Geant4 region. `UniformFieldFastSimulationModel` attaches the model to that
region when fast mode is enabled. `UniformFieldFastPropagator` computes the
helix state at a z plane, including position, momentum, path length, global
time, and proper time. The fast step ends 0.01 mm inside the next z face of
the fast volume; Geant4 continues ordinary transport after that point.

The model applies to every charged particle, including secondaries. At each eligible track
state, it checks whether the orbit fits inside the fast volume. If the
conditions are not met, `ModelTrigger` declines the fast step and Geant4 uses
ordinary transport.

## Requirements and limits

- The global field manager must contain a `G4UniformMagField` aligned with z.
- The model triggers only for tracks in `G4_Galactic` material.
- Fast-region root volumes must be `G4Box` or solid, full-phi `G4Tubs` with
  inner radius zero, and their z axis must be aligned to the field. Hollow or
  partial-phi tubes are rejected when the model is attached.
- The analytic propagator requires nonzero forward or backward z momentum
  toward a z face. Tracks that fail its checks use ordinary transport.
- Fast volumes must not overlap each other.
  General geometry overlaps are allowed only for placements made of the
  supplied vacuum material. The overlap checks run when
  `ConfigureAndAttach` is called, whether fast mode is enabled or not.
- **Physics is skipped during the fast step.** `DoIt` analytically updates the
  charged particle's state and time but does not run the normal Geant4 physics
  processes. The model advances proper time but does not sample a decay or
  generate its daughters inside the fast step. An unstable particle may only be 
  checked for decay after the shortcut, so a decay that should occur inside the
   volume will not have its correct position or daughters. Do not use this model 
  where decay in flight matters unless decay handling is added.

## Create a fast volume

Include `UniformFieldFastGeometry.hh` and keep the helper alive as part of the
detector construction. The box constructor takes **half-lengths**. The tube
constructor takes only an outer radius and z half-length. It creates a solid,
full-phi `G4Tubs` with inner radius zero; the API does not accept an inner
radius or phi span. Both constructors create their logical volume with
`G4_Galactic`, resolved internally through `G4NistManager`; no material
argument is required or accepted.

```cpp
#include "UniformFieldFastGeometry.hh"

#include "G4SystemOfUnits.hh"

// In DetectorConstruction, after creating the world logical volume:
fFastGeometry = std::make_unique<UniformFieldFastGeometry>(
    "FastVacuum", worldLogical, G4ThreeVector(),
    G4ThreeVector(100.0 * mm, 100.0 * mm, 150.0 * mm));
```

The base name `FastVacuum` creates `FastVacuum_Solid`,
`FastVacuum_Logical`, a physical volume named `FastVacuum`, and the region
`FastVacuum_Region`. The geometry helper supports axis-aligned placements.

For a cylindrical fast volume, use the tube constructor:

```cpp
fFastGeometry = std::make_unique<UniformFieldFastGeometry>(
    "FastVacuum", worldLogical, G4ThreeVector(),
    100.0 * mm, 150.0 * mm);  // outer radius, z half-length
```

## Attach and activate it

After constructing the geometry and installing the global magnetic field,
call `ConfigureAndAttach` from `ConstructSDandField`:

```cpp
#include "UniformFieldFastSimulationModel.hh"

void DetectorConstruction::ConstructSDandField() {
  // Set up the global field manager first.
  // ...
  UniformFieldFastSimulationModel::ConfigureAndAttach(fFastGeometry.get());
}
```

The physics list must register Geant4's fast-simulation process for the
charged particles. Call this after the base physics list has constructed its
processes:

```cpp
void PhysicsList::ConstructProcess() {
  G4VModularPhysicsList::ConstructProcess();
  UniformFieldFastSimulationModel::ActivateForChargedParticles();
}
```

Instantiate the messenger before executing UI macros, for example in
`main`:

```cpp
#include "UniformFieldFastMessenger.hh"

int main(int argc, char** argv) {
  UniformFieldFastMessenger::Instance();
  // Create the run manager and execute the selected macro.
}
```

Fast mode is off by default. Enable it before `/run/initialize` in a macro:

```text
/fastModel/enable true
/run/initialize
```

Use `/fastModel/enable false` to keep ordinary transport. The command is
available in Geant4's `PreInit` state.

## Add it to another CMake project

After finding Geant4, add this directory and link the exported target:

```cmake
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_subdirectory(path/to/CODEX_Magfield/FastModel)
target_link_libraries(my_simulation PRIVATE FastModel::UniformField)
```

The target exports the `FastModel/include` directory and requires C++17. The
component uses Geant4 only.

## Files

- `UniformFieldFastGeometry`build and validate fast-region geometry.
- `UniformFieldFastMessenger`provide `/fastModel/enable`.
- `UniformFieldFastSimulationModel` implement Geant4 model attachment, trigger checks, and fast-step updates.
- `UniformFieldFastPropagator` provide the analytic helix propagation calculation.
