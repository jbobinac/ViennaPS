#include <Geometries/psMakeHole.hpp>
#include <psMakeHole.hpp>
#include <psProcess.hpp>
#include <psToSurfaceMesh.hpp>
#include <psVTKWriter.hpp>

#include "SF6O2EtchingCompact.hpp"

int main(int argc, char *argv[]) {
  using NumericType = double;
  constexpr int D = 3;

  // Parameters
  int p = 25;
  int y = 44;
  int ionEnergy = 100;
  int t = 150;
  int mH = 12;

  if (argc > 1) {
    p = std::stoi(argv[1]);
    y = std::stoi(argv[2]);
    // ionEnergy = std::stoi(argv[3]);
    // t = std::stoi(argv[4]);
    // mH = std::stoi(argv[5]);
    }

  // Generate the initial geometry. SF6O2 Compact model is currently only
  // applicable to the geometry below. TO DO: extend the model to capture
  // geometry variations
  auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
  psMakeHole<NumericType, D>(geometry, 0.02 /* grid delta */, 1 /*x extent*/,
                             1 /*y extent*/, 0.175 /*hole radius*/,
                             mH/10. /* mask height*/, 1.193 /*taper angle*/,
                             true /*create mask*/)
      .apply();

  SF6O2CompactEtching<NumericType, D> model(
      "dimensionsP&y35ptGrid.csv" /*filename*/, p /*chamber pressure*/,
      y /*chamber oxygen percentage*/, ionEnergy /*ionEnergy*/,
      t /* etching time*/, mH /* mask height*/);

  psProcess<NumericType, D> process;
  process.setDomain(geometry);
  process.setProcessModel(model.getProcessModel());

  auto mesh = psSmartPointer<lsMesh<NumericType>>::New();
  psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
  psVTKWriter<NumericType>(mesh, "initial.vtp").apply();

  process.apply();

  psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
  // psVTKWriter<NumericType>(mesh,
  //                          "CMHole_p_" + std::to_string(p) + "_y_" +
  //                              std::to_string(y) + "_V_" +
  //                              std::to_string(int(ionEnergy)) + "_mH_" +
  //                              std::to_string(int(mH)) + "_" +
  //                              std::to_string(t) + ".vtp")
  //     .apply();
  psVTKWriter<NumericType>(mesh,
                          "CMHole_p_" + std::to_string(p) + "_y_" +
                              std::to_string(y) + ".vtp")
    .apply();

  return EXIT_SUCCESS;
}
