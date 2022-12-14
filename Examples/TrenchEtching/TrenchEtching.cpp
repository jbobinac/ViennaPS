#include <Geometries/psMakeHole.hpp>
#include <SF6O2Etching.hpp>
#include <psConfigParser.hpp>
#include <psProcess.hpp>
#include <psToSurfaceMesh.hpp>
#include <psVTKWriter.hpp>

int main(int argc, char *argv[]) {
  using NumericType = double;
  constexpr int D = 3;

  // Parse the parameters
  psProcessParameters<NumericType> params;
  if (argc > 1) {
    psConfigParser<NumericType> parser(argv[1]);
    parser.apply();
    params = parser.getParameters();
  }
  params.print();
  auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
  psMakeTrench<NumericType, D>(geometry, params.gridDelta, params.xExtent,
                               params.yExtent, params.trenchWidth,
                               params.trenchHeight, params.taperAngle).apply();

  SF6O2Etching<NumericType, D>
      model(params.totalIonFlux /*ion flux*/,
            params.totalEtchantFlux /*etchant flux*/,
            params.totalOxygenFlux /*oxygen flux*/,
            params.ionEnergy /*min ion energy (eV)*/,
            params.A_O /*oxy sputter yield*/,
            params.A_SiO /*oxide sputter yield*/, 0 /*mask material ID*/);

  psProcess<NumericType, D> process;
  process.setDomain(geometry);
  process.setProcessModel(model.getProcessModel());
  process.setMaxCoverageInitIterations(10);
  process.setNumberOfRaysPerPoint(20);

  auto mesh = psSmartPointer<lsMesh<NumericType>>::New();
  psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
  psVTKWriter<NumericType>(mesh, "initialTrench.vtp").apply();

  static constexpr int intermediateSteps = 3;
  for (size_t counter = 1; counter <= intermediateSteps; counter++) {
    std::cout << "Step: " << counter << "/" << intermediateSteps << std::endl;
    auto duration = params.processTime / intermediateSteps;
    process.setProcessDuration(duration);

    process.apply();

    psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
    psVTKWriter<NumericType>(
        mesh, "./Examples/TrenchEtching/MaskFaceting/TrenchEtch_mH_" // + std::to_string(params.y) 
        //+ "_tA_" + std::to_string(int(std::round(params.taperAngle *10))) + "_mH_"
                  + std::to_string(int(std::round(params.trenchHeight * 1000))) + "_" +
                  std::to_string(int(std::round(counter * duration))) + ".vtp")
        .apply();

    // psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
    // psVTKWriter<NumericType>(
    //     mesh, "./Output5D/HoleEtch_p_" + std::to_string(int(params.P * 10)) +
    //               "_y_" + std::to_string(int(params.y)) + "_V_" +
    //               std::to_string(int(params.ionEnergy)) + "_mH_" +
    //               std::to_string(int(params.maskHeight * 10)) + "_" +
    //               std::to_string(int(counter * duration)) + ".vtp")
    //     .apply();

    // psVTKWriter<NumericType>(
    //     mesh, "./ValidationFilesP&yOut/HoleEtch_F_" +
    //     std::to_string(params.totalEtchantFlux) + "_O_" +
    //               std::to_string(params.totalOxygenFlux) + "_V_" +
    //               std::to_string(int(params.ionEnergy)) + "_AO_" +
    //               std::to_string(params.A_O) + "_IF_" +
    //               std::to_string(params.totalIonFlux) +
    //               "_" + std::to_string(int(counter * duration)) + ".vtp")
    //     .apply();
  }

  return EXIT_SUCCESS;
}
