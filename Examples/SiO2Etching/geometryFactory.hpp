#ifndef GEOMFAC_HPP
#define GEOMFAC_HPP

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>

#define PI 3.141592

template <class T, int D> class MakeMask {
  using LSPtrType = lsSmartPointer<lsDomain<T, D>>;
  LSPtrType mask;

  std::array<T, 3> maskOrigin = {};
  T maskBaseRadius = 0.;
  T maskTopRadius = 0.;
  T maskHeight = 1.2;

public:
  MakeMask(LSPtrType passedMask) : mask(passedMask) {}

  void setMaskOrigin(std::array<T, 3> &origin) { maskOrigin = origin; }

  void setMaskBaseRadius(T radius) { maskBaseRadius = radius; }

  void setMaskTopRadius(T radius) { maskTopRadius = radius; }

  void setMaskHeight(T height) { maskHeight = height; }

  void apply() {
    auto &grid = mask->getGrid();
    auto &boundaryCons = grid.getBoundaryConditions();
    auto gridDelta = grid.getGridDelta();
    T extent = grid.getGridExtent(0) * gridDelta;

    // create mask
    {
      T normal[3] = {static_cast<T>(0.),
                     (D == 2) ? static_cast<T>(1.) : static_cast<T>(0.),
                     (D == 3) ? static_cast<T>(1.) : static_cast<T>(0.)};
      T origin[3] = {static_cast<T>(0.),
                     (D == 2) ? maskHeight : static_cast<T>(0.),
                     (D == 3) ? maskHeight : static_cast<T>(0.)};

      lsMakeGeometry<T, D>(mask,
                           lsSmartPointer<lsPlane<T, D>>::New(origin, normal))
          .apply();
      normal[D - 1] = static_cast<T>(-1.0);
      origin[D - 1] = static_cast<T>(0.);
      auto maskBottom = LSPtrType::New(grid);
      lsMakeGeometry<T, D>(maskBottom,
                           lsSmartPointer<lsPlane<T, D>>::New(origin, normal))
          .apply();
      lsBooleanOperation<T, D>(mask, maskBottom,
                               lsBooleanOperationEnum::INTERSECT)
          .apply();
  {
    std::cout << "Writing mask" << std::endl;
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<T, D>(mask, mesh).apply();
    lsVTKWriter(mesh, lsFileFormatEnum::VTP, "mask.vtp").apply();
  }
      // Create sidewalls
       auto maskHole = LSPtrType::New(grid);
       if( D == 2 )
      {
        double angleRad = 1.5499660;

        normal[0] = -1.;
        normal[1] = -1./tan(angleRad);
        origin[0] = -maskTopRadius;
        origin[1] = maskHeight;

        lsMakeGeometry<T, D>(maskHole,
                            lsSmartPointer<lsPlane<T, D>>::New(origin, normal))
            .apply();
        normal[0] = 1.;
        normal[1] = static_cast<T>(-1./tan(angleRad));
        origin[0] = static_cast<T>(maskTopRadius);
        auto maskRight = LSPtrType::New(grid);
        lsMakeGeometry<T, D>(maskRight,
                            lsSmartPointer<lsPlane<T, D>>::New(origin, normal))
            .apply();
        lsBooleanOperation<T, D>(maskHole, maskRight,
                                lsBooleanOperationEnum::INTERSECT)
            .apply();
      }
      else{
        T origin[D] = {0., 0., 0.};
        T maskOrigin[D] = {0., 0., 0.};
        T axisDirection[D] = {0., 0., 1.};

        lsMakeGeometry<T, D>(maskHole,
                    lsSmartPointer<lsPlane<T, D>>::New(maskOrigin, normal));

        auto cylinder = lsSmartPointer<lsCylinder<T, D>>::New(
            origin, axisDirection, maskHeight, maskBaseRadius, maskTopRadius);
        lsMakeGeometry<T, D>(maskHole, cylinder).apply();
      }

      lsBooleanOperation<T, D>(mask, maskHole,
                               lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
          .apply();
    }
  }
};

template <typename T, int D> class MakeLayer {
  using LSPtrType = lsSmartPointer<lsDomain<T, D>>;
  LSPtrType mask;

public:
  MakeLayer(lsSmartPointer<lsDomain<T, D>> passedMask) : mask(passedMask) {}

  LSPtrType apply() {
    auto &grid = mask->getGrid();
    auto materialLayer = LSPtrType::New(grid);
    T origin[D] = {0};
    origin[D - 1] = 0.;
    T normal[D] = {0};
    normal[D - 1] = 1.;
  
    auto plane = lsSmartPointer<lsPlane<T, D>>::New(origin, normal);
    lsMakeGeometry<T, D>(materialLayer, plane).apply();
    return materialLayer;
  }
};

#endif