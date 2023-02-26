#include <Geometries/psMakeHole.hpp>
#include <SiO2Etching.hpp>
#include <psConfigParser.hpp>
#include <psProcess.hpp>
#include <psToSurfaceMesh.hpp>
#include <psVTKWriter.hpp>

#include "geometryFactory.hpp"

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
  psMakeHole<NumericType, D>(
      geometry, params.gridDelta /* grid delta */, params.xExtent /*x extent*/,
      params.yExtent /*y extent*/, params.holeRadius /*hole radius*/,
      params.maskHeight /* mask height*/,
      params.taperAngle /* tapering angle in degrees */, true /*create mask*/)
      .apply();



  SiO2Etching<NumericType, D> model(params.totalIonFlux /*ion flux*/,
                                     params.totalEtchantFlux /*etchant flux*/,
                                     params.totalPolymerFlux /*polymer flux*/,
                                     params.ionEnergy /*min ion energy (eV)*/,
                                     0 /*mask material ID*/);

  // {auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
  // lsToSurfaceMesh<NumericType, D>(geometry->getLevelSets()->back(),
  // mesh).apply(); auto meshSubstrate =
  // lsSmartPointer<lsMesh<NumericType>>::New(); lsToSurfaceMesh<NumericType,
  // D>(geometry->getLevelSets()->at(0), meshSubstrate).apply();
  // lsVTKWriter<NumericType>(mesh, lsFileFormatEnum::VTP
  // ,"meshTest.vtp").apply(); lsVTKWriter<NumericType>(meshSubstrate,
  // lsFileFormatEnum::VTP ,"substrateTest.vtp").apply();}

  psProcess<NumericType, D> process;
  process.setDomain(geometry);
  process.setProcessModel(model.getProcessModel());
  process.setMaxCoverageInitIterations(10);
  process.setNumberOfRaysPerPoint(50);

  auto mesh = psSmartPointer<lsMesh<NumericType>>::New();
  psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
  psVTKWriter<NumericType>(mesh, "initial.vtp").apply();
  static constexpr int intermediateSteps = 2;
  for (size_t counter = 1; counter <= intermediateSteps; counter++) {
    std::cout << "Step: " << counter << "/" << intermediateSteps << std::endl;
    auto duration = params.processTime / intermediateSteps;

    process.setProcessDuration(duration);

    process.apply();

    psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();

    psVTKWriter<NumericType>(
        mesh, "./SiO2Etch_" +
                  std::to_string(int(std::round(counter * duration))) + ".vtp")
        .apply();
  }

  return EXIT_SUCCESS;
}
