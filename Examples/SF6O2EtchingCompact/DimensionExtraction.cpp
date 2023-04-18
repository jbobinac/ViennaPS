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
#include "CSVWriter.hpp"
#include "FilenameUtils.hpp"
// Auto-generated header which includes data directory path
#include "/home/bobinac/Documents/ViennaTools/ViennaPS/build/Examples/SF6O2EtchingCompact/parameters.hpp"

#include <queue>
#include <deque>
#include <iostream>

template <typename T, int MaxLen, typename Container=std::deque<T>>
class FixedQueue : public std::queue<T, Container> {
public:
    void push(const T& value) {
        if (this->size() == MaxLen) {
           this->c.pop_front();
        }
        std::queue<T, Container>::push(value);
    }
};


template <class T> struct Dimensions {
  T depth;
  T width;
};

int main(int argc, const char *argv[]) {
  using NumericType = double;
  static constexpr int D = 2;
  // set the eps to half the grid delta as a filter to find points at half of max depth
  static constexpr NumericType eps = 2.e-2;

  // Extract from mask height variation study
  // const std::string filenameRegex =
  //     "\\bHoleEtch_tA_([0-9]+)_m[y|Y]_([0-9]+)_([0-9]+)\\.vtp\\b";

  // const std::string filenameRegex =
  //     "\\bTrenchEtch_y_([0-9]+)_tA_([0-9]+)_mY_0_150\\.vtp\\b";

  const std::string filenameRegex =
    "\\bHoleEtch_p_([0-9]+)_[O|y]_([0-9]+)_100\\.vtp\\b";

  std::string dataDir = "data/"; //params::defaultDataDir;

  if (argc > 1) {
    dataDir = argv[1];
  }

  std::map<std::pair<NumericType, NumericType>, Dimensions<NumericType>>
      dimensionData;
  CSVWriter<NumericType> dimWriter("P&yTest.csv");
  int counter = 0;
  for (auto dirItem : std::filesystem::directory_iterator(dataDir)) {

    // Check whether the file is a regular file
    if (!dirItem.is_regular_file())
      continue;

    // Check whether the file is named according to our internal disk mesh
    // naming convention
    auto filePath = std::filesystem::path(dirItem);
    auto filename = std::string(filePath.filename());

    // auto [y, tapering, maskHeight] =
    //     extractParameters(filename, filenameRegex);

    auto [maskHeight, taperAngle] =
        extractParameters(filename, filenameRegex);

    // NumericType mY = std::stod(maskHeight) / 100.;
    // if (mY > 0.5)
    //   continue;

    counter++;
     
    //std::cout << tapering << "   " << maskHeight << "   " << timeSteps << std::endl;

    if (maskHeight.empty() || taperAngle.empty()) // || timeSteps.empty()
      continue;

    // NumericType taper;
    // // if (tapering.size() == 1)
    // //     tapering+= '0';
    // try {
    //     taper = std::stod(tapering) / 10.;
    // } catch (std::invalid_argument e) {
    //   continue;
    // } catch (std::out_of_range e) {
    //   continue;
    // }

    Dimensions<NumericType> dimensions;
   
    // Load the mesh file
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    try {
      lsVTKReader<NumericType>(mesh, filePath.generic_string()).apply();
    } catch (std::exception e) {
      std::cout << e.what() << std::endl;
      continue;
    }

    // std::cout << filename << "\n- tapering: " << taper
    //           << std::endl;

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

    // std::cout << "- max depth: " << maxDepth << std::endl;
    // std::cout << "- timeSteps: " << std::stod(timeSteps) << std::endl;
    // std::cout << "- mH: " << std::stod(maskHeight) << std::endl;

    dimensions.depth = std::fabs(maxDepth);

    // width at 1/2 depth - select points at 1/2 depth +/- gridDelta/2 and take an average radius
    std::vector<NumericType> radii2;

    // Set up the vector to get 20 largest radii
    // std::vector<NumericType> maxRadii(10, -1);
    // std::vector<NumericType> maxRadiiPos(10, 1);

    //NumericType maxRadiusDepth = 1;
    NumericType maxpos = -1;
    unsigned count = 0;
    for (unsigned i = 0; i < zpos.size(); ++i) {
      const auto &a = xpos[i];
      const auto &b = ypos[i];
      const auto &p = zpos[i];
      NumericType r = a*a + b*b;
      if (p < (maxDepth / 2) + eps && p > (maxDepth/2 - eps) ) { 
        radii2.push_back(r);
        count++;
      }
      // Find the maximum radius and its location
      // if (p < 0 && p > maxDepth * 0.6){
      //   if (r > maxRadii.back()){
      //     maxRadii.insert(maxRadii.begin(), r);
      //     maxRadiiPos.insert(maxRadiiPos.begin(), p);
      //     maxRadii.pop_back();
      //     maxRadiiPos.pop_back();
      //   }
      // }
    }
    dimensions.width = std::sqrt(std::accumulate(radii2.begin(), radii2.end(), 0.0) / radii2.size());
  
    // NumericType maxRadius = std::sqrt(std::accumulate(maxRadii.begin(), maxRadii.end(), 0.0) / maxRadii.size());
    // NumericType maxRadiusRelPos = std::accumulate(maxRadiiPos.begin(), maxRadiiPos.end(), 0.0) / maxRadiiPos.size();
    // maxRadiusRelPos = std::fabs(maxRadiusRelPos);// / dimensions.depth;

    // dimensionData[std::pair{params.spacerWidth - 2 * iso, sticking}] =
    //     dimensions;
    // dimWriter.writeRow(
    //     std::vector<NumericType>{ std::stod(maskHeight) /100., std::stod(taperAngle), maxRadius, maxRadiusRelPos, dimensions.depth});
    
    dimWriter.writeRow(
        std::vector<NumericType>{ std::stod(taperAngle) /10., std::stod(maskHeight), dimensions.depth, dimensions.width});
         // taper,   std::stod(maskHeight) / 100., dimensions.width
  }
  std::cout << counter << std::endl;
}
