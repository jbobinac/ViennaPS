#pragma once

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>

#include <psDomain.hpp>

/**
  Creates a plane in z(3D)/y(2D) direction.
*/
template <class NumericType, int D> class psMakePlane {
  using LSPtrType = psSmartPointer<lsDomain<NumericType, D>>;
  using PSPtrType = psSmartPointer<psDomain<NumericType, D>>;

public:
  PSPtrType domain = nullptr;

  NumericType gridDelta = .25;
  NumericType xExtent = 10;
  NumericType yExtent = 10;
  NumericType height = 0.;
  bool periodicBoundary = false;

  psMakePlane(PSPtrType passedDomain) : domain(passedDomain) {}

  psMakePlane(PSPtrType passedDomain, const NumericType passedGridDelta,
              const NumericType passedXExtent, const NumericType passedYExtent,
              const NumericType passedHeight, const bool passedPeriodic = false)
      : domain(passedDomain), gridDelta(passedGridDelta),
        xExtent(passedXExtent), yExtent(passedYExtent), height(passedHeight),
        periodicBoundary(passedPeriodic) {}

  void apply() {
    domain->clear();
    double bounds[2 * D];
    bounds[0] = -xExtent / 2.;
    bounds[1] = xExtent / 2.;

    if constexpr (D == 3) {
      bounds[2] = -yExtent / 2.;
      bounds[3] = yExtent / 2.;
      bounds[4] = -gridDelta;
      bounds[5] = gridDelta;
    } else {
      bounds[2] = -gridDelta;
      bounds[3] = gridDelta;
    }

    typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];

    for (int i = 0; i < D - 1; i++) {
      if (periodicBoundary) {
        boundaryCons[i] =
            lsDomain<NumericType, D>::BoundaryType::PERIODIC_BOUNDARY;
      } else {
        boundaryCons[i] =
            lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
      }
    }
    boundaryCons[D - 1] =
        lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

    auto substrate = LSPtrType::New(bounds, boundaryCons, gridDelta);
    NumericType normal[D] = {0.};
    NumericType origin[D] = {0.};
    normal[D - 1] = 1.;
    origin[D - 1] = height;
    lsMakeGeometry<NumericType, D>(
        substrate, lsSmartPointer<lsPlane<NumericType, D>>::New(origin, normal))
        .apply();
  }
};