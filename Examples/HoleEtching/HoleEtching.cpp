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

  auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
  psMakeHole<NumericType, D>(
      geometry, params.gridDelta /* grid delta */, params.xExtent /*x extent*/,
      params.yExtent /*y extent*/, params.holeRadius /*hole radius*/,
      params.maskHeight /* mask height*/,
      params.taperAngle /* tapering angle in degrees */, true /*create mask*/)
      .apply();

  SF6O2Etching<NumericType, D> model(params.totalIonFlux /*ion flux*/,
                                     params.totalEtchantFlux /*etchant flux*/,
                                     params.totalOxygenFlux /*oxygen flux*/,
                                     params.ionEnergy /*min ion energy (eV)*/,
                                     params.A_O /*oxy sputter yield*/,
                                     0 /*mask material ID*/);

  psProcess<NumericType, D> process;
  process.setDomain(geometry);
  process.setProcessModel(model.getProcessModel());
  process.setMaxCoverageInitIterations(10);
  process.setNumberOfRaysPerPoint(1000);
  process.setProcessDuration(params.processTime);

  auto mesh = psSmartPointer<lsMesh<NumericType>>::New();
  psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
  psVTKWriter<NumericType>(mesh, "initial.vtp").apply();

  process.apply();

  psToSurfaceMesh<NumericType, D>(geometry, mesh).apply();
  psVTKWriter<NumericType>(mesh, "final.vtp").apply();

  return EXIT_SUCCESS;
}
