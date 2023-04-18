#include <algorithm>
#include <array>
#include <filesystem>
#include <map>
#include <set>
#include <numeric>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsMesh.hpp>
#include <lsSmartPointer.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>
#include "CSVWriter.hpp"
#include "FilenameUtils.hpp"
// Auto-generated header which includes data directory path
#include "/home/bobinac/Documents/ViennaTools/ViennaPS/build/Examples/SF6O2EtchingCompact/parameters.hpp"


template <class T> struct Dimensions {
  T depth;
  std::vector<T> widths;
};

int main(int argc, const char *argv[]) {
  using NumericType = double;
  static constexpr int D = 2;
  // set the eps to half the grid delta as a filter to find points at half of max depth
  static constexpr NumericType eps = 2.e-2;

  // const std::string filenameRegex =
  //     "\\bHoleEtch_p_([0-9]+)_y_([0-9]+)_V_([0-9]+)_mH_([0-9]+)_([0-9]+)\\.vtp\\b";
    const std::string filenameRegex =
      "\\bHoleEtch_p_([0-9]+)_[O|y]_([0-9]+)_100\\.vtp\\b";

  std::string dataDir = params::defaultDataDir;

  if (argc > 1) {
    dataDir = argv[1];
  }

  std::map<std::tuple<NumericType, NumericType>, Dimensions<NumericType>> //, NumericType, NumericType, NumericType
      dimensionData;
  CSVWriter<NumericType> dimWriter("P&yValidation.csv");

  //dimWriter.writeLine("# P,y,ionEnergy,mH,t,maxDepth,widths");
  dimWriter.writeLine("# P,y,maxDepth,widths");
  std::vector<NumericType> sampleDepths{0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99};
  dimWriter.writePositionalData(sampleDepths);
  int count =0;
  for (auto dirItem : std::filesystem::directory_iterator(dataDir)) {
    //std::cout << count++ << std::endl;
    // Check whether the file is a regular file
    if (!dirItem.is_regular_file())
      continue;

    // Check whether the file is named according to our internal disk mesh
    // naming convention
    auto filePath = std::filesystem::path(dirItem);
    auto filename = std::string(filePath.filename());

    // auto [pressure, oxygenFraction, ionEnergy, t, maskHeight] =
    //     extractParameters(filename, filenameRegex);

    auto [pressure, oxygenFraction] =
        extractParameters(filename, filenameRegex);

    std::cout << pressure << "   " << oxygenFraction << std::endl; //"  " << ionEnergy << "  " << t << "  " << maskHeight << std::endl;

  
    if (pressure.empty() || oxygenFraction.empty())// || ionEnergy.empty() || t.empty() || maskHeight.empty()) 
      continue;
    NumericType P, y, ionE, T, mH;
    // Correct the pressure parsing error for 17.5 and 32.5
    P = std::stod(pressure);
    // if (P == 32 || P == 17)
    //   P += 0.5;
    y = std::stod(oxygenFraction);
    // ionE = std::stod(ionEnergy);
    // T = std::stod(t);
    // mH = std::stod(maskHeight);
// std::cout << 3 << std::endl;

     Dimensions<NumericType> dimensions;

    // Load the mesh file
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    try {
      lsVTKReader<NumericType>(mesh, filePath.generic_string()).apply();
    } catch (std::exception e) {
      std::cout << e.what() << std::endl;
      continue;
    }

    const auto &nodes = mesh->getNodes();

    std::vector<NumericType> xpos;
    std::vector<NumericType> ypos;
    std::vector<NumericType> zpos;
    xpos.reserve(nodes.size());
    ypos.reserve(nodes.size());
    zpos.reserve(nodes.size());

    std::transform(nodes.cbegin(), nodes.cend(), std::back_inserter(xpos),
                   [&](const auto &n) { return n[0]; });
    std::transform(nodes.cbegin(), nodes.cend(), std::back_inserter(ypos),
                   [&](const auto &n) { return n[1]; });
    std::transform(nodes.cbegin(), nodes.cend(), std::back_inserter(zpos),
                   [&](const auto &n) { return n[2]; });

    // bottom height
    auto maxDepth = *std::min_element(zpos.begin(), zpos.end());

    //std::cout << "- max depth: " << maxDepth << std::endl;

    dimensions.depth = std::fabs(maxDepth);

    for(int j = 0; j < sampleDepths.size(); ++j){
      // width at specific percentage of the total depth +/- gridDelta/2 and take an average radius
      std::vector<NumericType> radii2;
      NumericType maxpos = -1;
      unsigned count = 0;
      for (unsigned i = 0; i < zpos.size(); ++i) {
        const auto &a = xpos[i];
        const auto &b = ypos[i];
        const auto &p = zpos[i];
        if (p < (maxDepth * sampleDepths[j]) + eps && p > (maxDepth * sampleDepths[j] - eps) ) { 
          radii2.push_back(a*a + b*b);
          //count++;
        }
      }
      dimensions.widths.push_back(std::sqrt(std::accumulate(radii2.begin(), radii2.end(), 0.0) / radii2.size()));
    }

     dimensionData[std::tuple{y, P}] = //, ionE, T, mH
         dimensions;
//std::cout << 4 << std::endl;
    std::vector<NumericType> temp = {P, y, dimensions.depth}; //, ionE, T, mH
    for (auto i: dimensions.widths){
      temp.push_back(i);
    }
    //if(y == 58 || y== 53 || y == 60)
      dimWriter.writeRow(temp);
     count++;
   }
  std::cout << "Count is  " << count << std::endl;
}


// // Now store the extracted data in a Mesh
//   auto dataMesh = lsSmartPointer<lsMesh<>>::New();

//   std::set<NumericType> pressures;
//   std::set<NumericType> oxygenFractions;

//   std::vector<NumericType> depths;
//   std::vector<std::vector<NumericType>> widthData(sampleDepths.size());

//   for (const auto &[k, v] : dimensionData) {
//     depths.push_back(v.depth);
//     for (int w = 0; w < v.widths.size(); ++w){
//       widthData[w].push_back(v.widths.at(w));
//     }
//     dataMesh->insertNextNode(std::array<NumericType, 3>{k.first, k.second, 0.});
//     oxygenFractions.insert(k.first);
//     pressures.insert(k.second);
//   }

//   unsigned nYs = oxygenFractions.size();
//   unsigned nPressures = pressures.size();
//   for (unsigned j = 0; j < nYs - 1; ++j)
//     for (unsigned i = 0; i < nPressures - 1; ++i) {
//       unsigned index = j * nPressures + i;
//       dataMesh->insertNextTetra(std::array<unsigned, 4>{
//           index, index + nPressures, index + nPressures + 1, index + 1});
//     }

//   dataMesh->getPointData().insertNextScalarData(depths,
//                                                 "depths");     
//   dataMesh->getPointData().insertNextScalarData(widthData[0], "widths005");                                                                                      
//   dataMesh->getPointData().insertNextScalarData(widthData[1], "widths005");
//   dataMesh->getPointData().insertNextScalarData(widthData[2], "widths01");
//   dataMesh->getPointData().insertNextScalarData(widthData[3], "widths02");
//   dataMesh->getPointData().insertNextScalarData(widthData[4], "widths03");
//   dataMesh->getPointData().insertNextScalarData(widthData[5], "widths04");
//   dataMesh->getPointData().insertNextScalarData(widthData[6], "widths05");
//   dataMesh->getPointData().insertNextScalarData(widthData[7], "widths06");
//   dataMesh->getPointData().insertNextScalarData(widthData[8], "widths07");
//   dataMesh->getPointData().insertNextScalarData(widthData[9], "widths08");
//   dataMesh->getPointData().insertNextScalarData(widthData[10], "widths09");
//   dataMesh->getPointData().insertNextScalarData(widthData[11], "widths095");

//   lsVTKWriter<NumericType>(dataMesh, "dataMesh.vtu").apply();
//}
