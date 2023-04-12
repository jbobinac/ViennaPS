#include <lsVTKWriter.hpp>
#include <psCSVDataSource.hpp>
#include <psCSVWriter.hpp>
#include <psDataScaler.hpp>
#include <psMakeHole.hpp>
#include <psNearestNeighborsInterpolation.hpp>
#include <psProcessModel.hpp>
#include <psRectilinearGridInterpolation.hpp>
#include <psToSurfaceMesh.hpp>
#include <psVTKWriter.hpp>

template <typename NumericType, int D>
class SF6O2CompactGeomModel : public psGeometricModel<NumericType, D> {
  using psGeometricModel<NumericType, D>::domain;
  NumericType pressure;
  NumericType chamberOxyPercentage;
  NumericType ionEnergy;
  NumericType t;
  NumericType maskHeight;

public:
  SF6O2CompactGeomModel(std::string passedFilename,
                        const NumericType passedPressure,
                        const NumericType passedChamberOxyPercentage,
                        const NumericType passedIonEnergy,
                        const NumericType passedMaskHeight,
                        const NumericType passedTime)
      : pressure(passedPressure),
        chamberOxyPercentage(passedChamberOxyPercentage),
        ionEnergy(passedIonEnergy), maskHeight(passedMaskHeight), t(passedTime) {}

  void apply() {

    static constexpr int InputDim = 2;
    static constexpr int TargetDim = 15;

    static constexpr int DataDim = InputDim + TargetDim;

    auto dataSource =
        psSmartPointer<psCSVDataSource<NumericType, DataDim>>::New();
    dataSource->setFilename("dimensionsP&y35ptGrid.csv");

    psRectilinearGridInterpolation<NumericType, InputDim, TargetDim>
        interpolation;

    interpolation.setDataSource(dataSource);
    interpolation.initialize();
    std::array<NumericType, InputDim> inputParams{
        pressure, chamberOxyPercentage};//, ionEnergy, maskHeight, t};

    auto [value, isInside] = interpolation.estimate(inputParams);
    for (auto a :value)
      std::cout << a << " ";
    std::cout << std::endl;

    // Construct a final hole geometry from parameters
    auto cloud = lsSmartPointer<lsPointCloud<double, D>>::New();
    auto pointMesh = lsSmartPointer<lsMesh<>>::New();

    constexpr double limit = 2 * M_PI - 1e-6;
    unsigned numPoints = 30;
    double smallAngle = 2.0 * M_PI / double(numPoints);
    auto sampleDepths = dataSource->getPositionalParameters();
    // create and insert points at base
    unsigned count = 0;
    for (int i = 1; i < TargetDim; ++i) {
      for (int j = 0; j < numPoints; ++j) {
        double angle = 2.0 * M_PI * j / double(numPoints);
        NumericType x = value[i] * std::cos(angle);
        NumericType y = value[i] * std::sin(angle);
        NumericType z = -value[0] * sampleDepths[i - 1];
        cloud->insertNextPoint(hrleVectorType<double, D>(x, y, z));
        // Add the same points to the pointMesh for debugging purposes
        pointMesh->insertNextNode(std::array<double, D>{x, y, z});
        pointMesh->insertNextVertex(std::array<unsigned, 1>{count++});
      }
    }
    cloud->insertNextPoint(hrleVectorType<double, D>(0, 0, -value[0]));
    pointMesh->insertNextNode(std::array<double, D>{0, 0, -value[0]});
    pointMesh->insertNextVertex(std::array<unsigned, 1>{count++});

    // Create the etched shape from point cloud
    auto hull = lsSmartPointer<lsMesh<>>::New();
    lsConvexHull<double, D>(hull, cloud).apply();
    lsVTKWriter<double>(pointMesh, lsFileFormatEnum::VTP,
                        "pointMesh_" + std::to_string(1) + ".vtp")
        .apply();

    // Create a new level set on the same grid
    auto etchedMaterial = lsSmartPointer<lsDomain<NumericType, D>>::New(
        domain->getLevelSets()->back()->getGrid());

    lsFromSurfaceMesh<NumericType, D>(etchedMaterial, hull, false).apply();

    // Boolean operation
    auto levelSets = domain->getLevelSets();
    lsBooleanOperation<NumericType, D>(
        levelSets->back(), etchedMaterial,
        lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

    for (unsigned i = 0; i < levelSets->size() - 1; ++i)
      lsBooleanOperation<NumericType, D>(levelSets->at(i), levelSets->back(),
                                         lsBooleanOperationEnum::INTERSECT)
          .apply();
  }
};

template <typename NumericType, int D> class SF6O2CompactEtching {
  psSmartPointer<psProcessModel<NumericType, D>> processModel = nullptr;

public:
  SF6O2CompactEtching(std::string filename, const NumericType pressure,
                      const NumericType chamberOxyPercentage,
                      const NumericType ionEnergy = 100, const NumericType time = 100, const NumericType maskHeight = 12) {
    processModel = psSmartPointer<psProcessModel<NumericType, D>>::New();

    auto geomModel = psSmartPointer<SF6O2CompactGeomModel<NumericType, D>>::New(
        filename, pressure, chamberOxyPercentage, ionEnergy, maskHeight, time);

    processModel->setGeometricModel(geomModel);
    processModel->setProcessName("SF6O2EtchingCM");
  }

  psSmartPointer<psProcessModel<NumericType, D>> getProcessModel() {
    return processModel;
  }
};
