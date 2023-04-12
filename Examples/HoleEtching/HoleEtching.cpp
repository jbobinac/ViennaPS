#include <Geometries/psMakeHole.hpp>
#include <SF6O2Etching.hpp>
#include <psProcess.hpp>
#include <psToSurfaceMesh.hpp>
#include <psUtils.hpp>

#include "Parameters.hpp"
#include <psVTKWriter.hpp>

#include "geometryFactory.hpp"

int main(int argc, char *argv[]) {
  using NumericType = double;
  constexpr int D = 3;

  // Parse the parameters
  int P, y;

  Parameters<NumericType> params;
  if (argc > 1) {
    auto config = psUtils::readConfigFile(argv[1]);
    if (config.empty()) {
      std::cerr << "Empty config provided" << std::endl;
      return -1;
    }
    params.fromMap(config);
  }
  params.print();
  auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
  psMakeHole<NumericType, D>(
      geometry, params.gridDelta /* grid delta */, params.xExtent /*x extent*/,
      params.yExtent /*y extent*/, params.holeRadius /*hole radius*/,
      params.maskHeight /* mask height*/,
      params.taperAngle /* tapering angle in degrees */, 0 /* base height */,
      false /* periodic boundary */, true /*create mask*/)
      .apply();


  SF6O2Etching<NumericType, D> model(params.totalIonFlux /*ion flux*/,
                                     params.totalEtchantFlux /*etchant flux*/,
                                     params.totalOxygenFlux /*oxygen flux*/,
                                     params.ionEnergy /*min ion energy (eV)*/,
                                     params.A_O /*oxy sputter yield*/,
                                     params.A_SiO /*oxide sputter yield*/,
                                     0 /*mask material ID*/);

  // {auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
  // lsToSurfaceMesh<NumericType, D>(geometry->getLevelSets()->back(), mesh).apply();
  // auto meshSubstrate = lsSmartPointer<lsMesh<NumericType>>::New();
  // lsToSurfaceMesh<NumericType, D>(geometry->getLevelSets()->at(0), meshSubstrate).apply();
  // lsVTKWriter<NumericType>(mesh, lsFileFormatEnum::VTP ,"meshTest.vtp").apply();
  // lsVTKWriter<NumericType>(meshSubstrate, lsFileFormatEnum::VTP ,"substrateTest.vtp").apply();}

  psProcess<NumericType, D> process;
  process.setDomain(geometry);
  process.setProcessModel(model.getProcessModel());
  process.setMaxCoverageInitIterations(10);
  process.setNumberOfRaysPerPoint(50);

  auto mesh = psSmartPointer<lsMesh<NumericType>>::New();
  psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
  psVTKWriter<NumericType>(mesh, "initial.vtp").apply();

  static constexpr int intermediateSteps = 3;
  for (size_t counter = 1; counter <= intermediateSteps; counter++) {
    std::cout << "Step: " << counter << "/" << intermediateSteps
              << std::endl;
    auto duration = params.processTime / intermediateSteps;
    process.setProcessDuration(duration);

    process.apply();

    psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
    psVTKWriter<NumericType>(
        mesh, "./Examples/HoleEtching/MaskFaceting/HoleEtch_y_" + std::to_string(params.y) 
        + "_tA_" + std::to_string(int(std::round(params.taperAngle *10))) +
                  "_mH_" + std::to_string(int(std::round(params.maskHeight *10))) + "_" +
                  std::to_string(int(std::round(counter * duration))) + ".vtp")
        .apply();

    // psVTKWriter<NumericType>(
    //     mesh, "./ValidationFilesP&yOut/HoleEtch_F_" + std::to_string(params.totalEtchantFlux) + "_O_" +
    //               std::to_string(params.totalOxygenFlux) + "_V_" +
    //               std::to_string(int(params.ionEnergy)) + "_AO_" +
    //               std::to_string(params.A_O) + "_IF_" +
    //               std::to_string(params.totalIonFlux) +
    //               "_" + std::to_string(int(counter * duration)) + ".vtp")
    //     .apply();
  }

  return EXIT_SUCCESS;
}
