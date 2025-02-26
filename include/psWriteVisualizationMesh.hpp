#pragma once

#include <lsToSurfaceMesh.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include <psDomain.hpp>

template <class NumericType, int D> class psWriteVisualizationMesh {
public:
  psWriteVisualizationMesh() {}
  psWriteVisualizationMesh(
      psSmartPointer<psDomain<NumericType, D>> passedDomain,
      std::string passedFileName)
      : domain(passedDomain), fileName(passedFileName) {}

  void apply() {
    lsWriteVisualizationMesh<NumericType, D> visMesh;
    visMesh.setFileName(fileName);
    int i = 0;
    for (auto ls : *domain->getLevelSets()) {
      visMesh.insertNextLevelSet(ls);
      auto mesh = psSmartPointer<lsMesh<NumericType>>::New();
      lsToSurfaceMesh<NumericType, D>(ls, mesh).apply();
      psVTKWriter<NumericType>(mesh, "visMesh_" + std::to_string(i++) + ".vtp")
          .apply();
    }
    visMesh.apply();
  }

  void setFileName(std::string passedFileName) { fileName = passedFileName; }

  void setDomain(psSmartPointer<psDomain<NumericType, D>> passedDomain) {
    domain = passedDomain;
  }

private:
  psSmartPointer<psDomain<NumericType, D>> domain;
  std::string fileName;
};