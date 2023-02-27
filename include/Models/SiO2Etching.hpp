#pragma once

#include <csTracingParticle.hpp>

#include <rayParticle.hpp>
#include <rayReflection.hpp>
#include <rayUtil.hpp>

#include <psProcessModel.hpp>
#include <psSmartPointer.hpp>
#include <psSurfaceModel.hpp>
#include <psToSurfaceMesh.hpp>
#include <psVelocityField.hpp>

template <typename NumericType, int D>
class SiO2SurfaceModel : public psSurfaceModel<NumericType> {
  using psSurfaceModel<NumericType>::Coverages;
  //using psSurfaceModel<NumericType>::domain;

  NumericType totalIonFlux = 1e16;
  NumericType totalEtchantFlux = 2.5e17;
  NumericType totalPolymerFlux = 1.0e17;       // 5.0e16;
  static constexpr NumericType inv_rho_p =
      5.0e-23; // in (atoms/cm³)⁻¹ (rho Si)
  static constexpr NumericType inv_rho_SiO2 = 1. / (2.6e22);
  static constexpr double k_ev = 2.;
  static constexpr double k_ei = 2.;

    // Temperature
  static constexpr NumericType T = 300; // K
  static constexpr NumericType E_SiO2 = 0.168;
  static constexpr NumericType K_SiO2 = 0.0027;
  static constexpr NumericType kB = 0.000086173324; // m2 kg s-2 K-1

public:
  SiO2SurfaceModel(const double ionFlux, const double etchantFlux,
                    const double polymerFlux)
      : totalIonFlux(ionFlux), totalEtchantFlux(etchantFlux),
        totalPolymerFlux(polymerFlux) {}

  void initializeCoverages(unsigned numGeometryPoints) override {
    if (Coverages == nullptr) {
      Coverages = psSmartPointer<psPointData<NumericType>>::New();
    } else {
      Coverages->clear();
    }
    std::vector<NumericType> cov(numGeometryPoints);
    Coverages->insertNextScalarData(cov, "activeSitesPolymerCov");
    Coverages->insertNextScalarData(cov, "polymerCoverage");
    Coverages->insertNextScalarData(cov, "etchantCoverage");
  }

  psSmartPointer<std::vector<NumericType>>
  calculateVelocities(psSmartPointer<psPointData<NumericType>> Rates,
                      const std::vector<NumericType> &materialIds) override {

    updateCoverages(Rates);

    const auto numPoints = Rates->getScalarData(0)->size();
    std::vector<NumericType> surfaceRates(numPoints, 0.);

    const auto ionSputteringRate = Rates->getScalarData("ionSputteringRate");
    const auto ionEnhancedEtchantRate =
        Rates->getScalarData("ionEnhancedEtchantRate");
    const auto ionEnhancedPolymerRate =
        Rates->getScalarData("ionEnhancedPolymerRate");
    const auto polymerRate = Rates->getScalarData("polymerRate");
    const auto dielectricEtchantRate =
        Rates->getScalarData("dielectricEtchantRate");
    const auto polymerEtchantRate = Rates->getScalarData("polymerEtchantRate");
    //const auto evaporationRate = Rates->getScalarData("evaporationRate");
        auto eToEvFluxConst = K_SiO2 * std::exp(- E_SiO2 / kB * T);


    const auto activeSitesPolymerCov =
        Coverages->getScalarData("activeSitesPolymerCov");
    const auto polymerCoverage = Coverages->getScalarData("polymerCoverage");
    const auto etchantCoverage = Coverages->getScalarData("etchantCoverage");

    for (size_t i = 0; i < numPoints; ++i) {
      
       if (polymerCoverage->at(i) > 1.){ //Deposition
          auto ER_p = (ionEnhancedPolymerRate->at(i) * totalIonFlux * activeSitesPolymerCov->at(i)) * inv_rho_p * 1e4 ;
          auto DR_p = polymerRate->at(i) * totalPolymerFlux * inv_rho_p * 1e4 ;
          if(materialIds[i] >= 1){
            surfaceRates[i] = (DR_p - ER_p); 
          }
          else
            surfaceRates[i] = 0.;
        }
        else{ // Etching
            if(materialIds[i] >= 1){
              surfaceRates[i] = - inv_rho_SiO2 * 1e4 * (ionEnhancedEtchantRate->at(i) * totalIonFlux * etchantCoverage->at(i) +
                          ionSputteringRate->at(i) * totalIonFlux * (1. - etchantCoverage->at(i)) +
                           eToEvFluxConst * dielectricEtchantRate->at(i) * totalEtchantFlux * etchantCoverage->at(i));; 
            }
            else
              surfaceRates[i] = 0.;
        }
    }

    return psSmartPointer<std::vector<NumericType>>::New(surfaceRates);
  }

  void
  updateCoverages(psSmartPointer<psPointData<NumericType>> Rates) override {
    // update coverages based on fluxes
    const auto numPoints = Rates->getScalarData(0)->size();

     const auto ionEnhancedEtchantRate = Rates->getScalarData("ionEnhancedEtchantRate");
    const auto ionEnhancedPolymerRate = Rates->getScalarData("ionEnhancedPolymerRate");
    const auto polymerRate = Rates->getScalarData("polymerRate");
    const auto dielectricEtchantRate = Rates->getScalarData("dielectricEtchantRate");
    const auto polymerEtchantRate = Rates->getScalarData("polymerEtchantRate");
    //const auto evaporationRate = Rates->getScalarData("evaporationRate");
    // Get the coefficient to convert etchant flux to evaporation flux
    auto eToEvFluxConst = K_SiO2 * std::exp(- E_SiO2 / kB * T);

    // Active sites on polymer coverage e,P
    auto activeSitesPolymerCov = Coverages->getScalarData("activeSitesPolymerCov");
    activeSitesPolymerCov->resize(numPoints);
    // Polymer coverage
    auto polymerCoverage = Coverages->getScalarData("polymerCoverage");
    polymerCoverage->resize(numPoints);
    // Etchant coverage
    auto etchantCoverage = Coverages->getScalarData("etchantCoverage");
    etchantCoverage->resize(numPoints);

    for (size_t i = 0; i < numPoints; ++i)
    {
      //Active sites on polymer coverage e,P
      activeSitesPolymerCov->at(i) = (polymerEtchantRate->at(i) * totalEtchantFlux == 0)?0. :
      polymerEtchantRate->at(i) * totalEtchantFlux /(polymerEtchantRate->at(i) * totalEtchantFlux + ionEnhancedPolymerRate->at(i) * totalIonFlux);

      //Polymer coverage
      polymerCoverage->at(i) = (polymerRate->at(i) * totalPolymerFlux == 0.)?0.:
                              ((activeSitesPolymerCov->at(i) == 0.) || (ionEnhancedPolymerRate->at(i) * totalIonFlux==0.))?1. :
                                polymerRate->at(i) * totalPolymerFlux /(ionEnhancedPolymerRate->at(i) * totalIonFlux * activeSitesPolymerCov->at(i));

      //Neutral coverage
      if (polymerCoverage->at(i) < 1.) {
        etchantCoverage->at(i) = (dielectricEtchantRate->at(i) * totalEtchantFlux == 0)?0:
                                  dielectricEtchantRate->at(i) * totalEtchantFlux * (1.- polymerCoverage->at(i)) / 
                                  (dielectricEtchantRate->at(i) * totalEtchantFlux + k_ei * ionEnhancedEtchantRate->at(i) * totalIonFlux
                                  + k_ev * eToEvFluxConst * dielectricEtchantRate->at(i) * totalEtchantFlux);
      } else {
        etchantCoverage->at(i) = 0.;
      }
            // if (i <= 10)
            //   std::cout << "   "
            //   << "P_rate = " << polymerRate->at(i) 
            //   << ", IEP_r = " << ionEnhancedPolymerRate->at(i)
            //   << ", IEEr = " << ionEnhancedEtchantRate->at(i)
            //   << ", diEtchtRate = " << dielectricEtchantRate->at(i)
            //   << ", polyERate = " << polymerEtchantRate->at(i)
            //   << ",  aSPCov = " << activeSitesPolymerCov->at(i) 
            //   << ",  Etch Coverage = " << etchantCoverage->at(i)
            //   << ",  Poly Coverage = " << polymerCoverage->at(i) << std::endl;
    }
  }
};

template <class T> class SiO2VelocityField : public psVelocityField<T> {
public:
  SiO2VelocityField(const int maskId) : maskMaterial(maskId) {}

  T getScalarVelocity(const std::array<T, 3> & /*coordinate*/, int material,
                      const std::array<T, 3> & /*normalVector*/,
                      unsigned long pointID) override {
    //    if (material != maskMaterial)
    return velocities->at(pointID);
    // else
    //   return 0.;
  }

  void setVelocities(psSmartPointer<std::vector<T>> passedVelocities) override {
    velocities = passedVelocities;
  }

private:
  psSmartPointer<std::vector<T>> velocities = nullptr;
  const int maskMaterial = 0;
};

template <typename NumericType, int D>
class SiO2Ion : public rayParticle<SiO2Ion<NumericType, D>, NumericType> {
public:
  SiO2Ion(NumericType passedEnergy = 100.)
      : minIonEnergy(passedEnergy) {}

  void surfaceCollision(NumericType rayWeight,
                        const rayTriple<NumericType> &rayDir,
                        const rayTriple<NumericType> &geomNormal,
                        const unsigned int primID, const int materialId,
                        rayTracingData<NumericType> &localData,
                        const rayTracingData<NumericType> *globalData,
                        rayRNG &Rng) override final {
                              
    // collect data for this hit
    assert(primID < localData.getVectorData(0).size() && "id out of bounds");
    assert(E >= 0 && "Negative energy ion");

    const auto cosTheta = -rayInternal::DotProduct(rayDir, geomNormal);

    assert(cosTheta >= 0 && "Hit backside of disc");
    assert(cosTheta <= 1 + 1e6 && "Error in calculating cos theta");

    //const double angle = std::acos(cosTheta);

    // Ion Enhanced
    f_e_ie = cosTheta;
    // Etchant sputtering
    f_e_sp = std::max((1 + B_sp * (1.- cosTheta * cosTheta)) * cosTheta, 0.);

    const double sqrtE = std::sqrt(E);

    const double Y_sp = Ae_sp * std::max(sqrtE - std::sqrt(Eth_e_sp), 0.) * f_e_sp;
    const double Ye_ie = Ae_ie * std::max(sqrtE - std::sqrt(Eth_e_e), 0.) * f_e_ie;
    const double Yp_ie = Ap_ie * std::max(sqrtE - std::sqrt(Eth_e_p), 0.) * f_e_ie;

    // sputtering yield Y_sp ionSputteringRate
    localData.getVectorData(0)[primID] += rayWeight * Y_sp;

    // ion enhanced etching yield Ye_ie ionEnhancedEtchantRate
    localData.getVectorData(1)[primID] += rayWeight * Ye_ie;

    // ion enhanced polymer yield Yp_ie ionEnhancedPolymerRate
    localData.getVectorData(2)[primID] += rayWeight * Yp_ie;
  }

  std::pair<NumericType, rayTriple<NumericType>>
  surfaceReflection(NumericType rayWeight, const rayTriple<NumericType> &rayDir,
                    const rayTriple<NumericType> &geomNormal,
                    const unsigned int primId, const int materialId,
                    const rayTracingData<NumericType> *globalData,
                    rayRNG &Rng) override final {
    auto cosTheta = -rayInternal::DotProduct(rayDir, geomNormal);

    assert(cosTheta >= 0 && "Hit backside of disc");
    assert(cosTheta <= 1 + 1e-6 && "Error in calculating cos theta");

    const NumericType incAngle =
        std::acos(std::max(std::min(cosTheta, 1.), 0.));

    NumericType Eref_peak = 0;

    // Small incident angles are reflected with the energy fraction centered at
    // 0
    if (incAngle >= inflectAngle) {
      Eref_peak =
          Eref_max *
          (1 - (1 - A) * std::pow((rayInternal::PI / 2 - incAngle) /
                                      (rayInternal::PI / 2 - inflectAngle),
                                  n_r));
    } else {
      Eref_peak = Eref_max * A * std::pow(incAngle / inflectAngle, n_l);
    }
    // Gaussian distribution around the Eref_peak scaled by the particle energy
    NumericType tempEnergy = Eref_peak * E;

    NumericType NewEnergy;

    std::uniform_real_distribution<NumericType> uniDist;

    do {
      const auto rand1 = uniDist(Rng);
      const auto rand2 = uniDist(Rng);
      NewEnergy = tempEnergy +
                  (std::min((E - tempEnergy), tempEnergy) + E * 0.05) *
                      (1 - 2. * rand1) * std::sqrt(std::fabs(std::log(rand2)));

    } while (NewEnergy > E || NewEnergy < 0.);

    // Set the flag to stop tracing if the energy is below the threshold
    if (NewEnergy > 1) {
      E = NewEnergy;
      auto direction = rayReflectionConedCosine<NumericType, D>(
          rayInternal::PI / 2. - std::min(incAngle, minAngle), rayDir,
          geomNormal, Rng);
      return std::pair<NumericType, rayTriple<NumericType>>{NumericType(0),
                                                            direction};
    } else {
      return std::pair<NumericType, rayTriple<NumericType>>{
          1., rayTriple<NumericType>{0., 0., 0.}};
    }
  }

  void initNew(rayRNG &RNG) override final {
    std::uniform_real_distribution<NumericType> uniDist;
    do {
      const auto rand = uniDist(RNG);

      E = minIonEnergy + deltaIonEnergy * rand;
    } while (E < 0);
  }

  int getRequiredLocalDataSize() const override final { return 3; }
  NumericType getSourceDistributionPower() const override final { return 1000.; }
  std::vector<std::string> getLocalDataLabels() const override final {
    return std::vector<std::string>{"ionSputteringRate", "ionEnhancedEtchantRate", "ionEnhancedPolymerRate"};
  }

private:
	static constexpr NumericType Ae_ie = 0.0361;
	static constexpr NumericType Ae_sp = 0.0139;
	static constexpr NumericType Ap_ie = 0.1444;

	static constexpr NumericType B_sp = 9.3;

  static constexpr NumericType Eth_e_e = 4.;
  static constexpr NumericType Eth_e_p = 4.;
  static constexpr NumericType Eth_e_sp = 18.;

  NumericType minIonEnergy = 50.;
  static constexpr NumericType deltaIonEnergy = 40.;

  static constexpr NumericType inflectAngle = 1.55334;
  static constexpr NumericType n_l = 10.;
  static constexpr NumericType n_r = 1.;
  static constexpr NumericType Eref_max = 1.;

  static constexpr NumericType minAngle = 1.3962634;

  static constexpr NumericType A =
      1. / (1. + (n_l / n_r) * (rayInternal::PI / (2 * inflectAngle) - 1.));

  NumericType f_e_ie;
  NumericType f_e_sp;
  NumericType E;
};

template <typename NumericType, int D>
class SiO2Etchant
    : public rayParticle<SiO2Etchant<NumericType, D>, NumericType> {
public:
  void surfaceCollision(NumericType rayWeight,
                        const rayTriple<NumericType> &rayDir,
                        const rayTriple<NumericType> &geomNormal,
                        const unsigned int primID, const int materialId,
                        rayTracingData<NumericType> &localData,
                        const rayTracingData<NumericType> *globalData,
                        rayRNG &Rng) override final {

    // etchant Rate
    if(materialId >= 1){
      // Neutrals coverages
      const auto &phi_ep = globalData->getVectorData(0)[primID];
      const auto &phi_e = globalData->getVectorData(2)[primID];

      // Obtain the sticking probability
      NumericType Se_e = gamma_ee * std::max(1. - phi_e, 0.);
      NumericType Se_p = gamma_ep * std::max(1. - phi_ep, 0.);

      localData.getVectorData(0)[primID] += rayWeight * Se_e ;
      localData.getVectorData(1)[primID] += rayWeight * Se_p ;
      // evaporation flux Jev evaporationFlux
      //localData.getVectorData(2)[primID] =  K_SiO2 * std::exp(- E_SiO2 / kB * T);
      }
  }
  std::pair<NumericType, rayTriple<NumericType>>
  surfaceReflection(NumericType rayWeight, const rayTriple<NumericType> &rayDir,
                    const rayTriple<NumericType> &geomNormal,
                    const unsigned int primID, const int materialId,
                    const rayTracingData<NumericType> *globalData,
                    rayRNG &Rng) override final {
    // Neutrals coverages
    const auto &phi_ep = globalData->getVectorData(0)[primID];
    const auto &phi_e = globalData->getVectorData(2)[primID];

    auto direction = rayReflectionDiffuse<NumericType, D>(geomNormal, Rng);
    NumericType weightToDrop = 0;
    if(materialId == 1)
      weightToDrop = gamma_ee * std::max(1. - phi_e, 0.);
    if(materialId == 2)
      weightToDrop = gamma_ep * std::max(1. - phi_ep, 0.);

    //auto direction = rayReflectionSpecular<NumericType>(rayDir, geomNormal);
    return std::pair<NumericType, rayTriple<NumericType>>{weightToDrop, direction};
  }
  void initNew(rayRNG &RNG) override final {}
  int getRequiredLocalDataSize() const override final { return 3; }
  NumericType getSourceDistributionPower() const override final { return 1.; }
  std::vector<std::string> getLocalDataLabels() const override final {
    return std::vector<std::string>{"dielectricEtchantRate", "polymerEtchantRate", "evaporationRate"};
  }

private:
  static constexpr NumericType gamma_ee = 0.9;
  static constexpr NumericType gamma_ep = 0.6;

};

template <typename NumericType, int D>
class SiO2Polymer
    : public rayParticle<SiO2Polymer<NumericType, D>, NumericType> {
public:
 void surfaceCollision(NumericType rayWeight,
                        const rayTriple<NumericType> &rayDir,
                        const rayTriple<NumericType> &geomNormal,
                        const unsigned int primID, const int materialId,
                        rayTracingData<NumericType> &localData,
                        const rayTracingData<NumericType> *globalData,
                        rayRNG &Rng) override final {

    // Active sites on polymer coverage
    const auto &phi_ep = globalData->getVectorData(0)[primID];
    // Polymer surface coverage
    const auto &phi_p = globalData->getVectorData(1)[primID];
    // Etchant surface coverage
    const auto &phi_e = globalData->getVectorData(2)[primID];

    const auto S_p = std::max(gamma_p * (1.- phi_ep * phi_p - phi_e ),0.);
    localData.getVectorData(0)[primID] += rayWeight * S_p;
  }

  std::pair<NumericType, rayTriple<NumericType>>
  surfaceReflection(NumericType rayWeight, const rayTriple<NumericType> &rayDir,
                    const rayTriple<NumericType> &geomNormal,
                    const unsigned int primID, const int materialId,
                    const rayTracingData<NumericType> *globalData,
                    rayRNG &Rng) override final {

    auto direction = rayReflectionSpecular<NumericType>(rayDir, geomNormal);

    // Active sites on polymer coverage
    const auto &phi_ep = globalData->getVectorData(0)[primID];
    // Polymer surface coverage
    const auto &phi_p = globalData->getVectorData(1)[primID];
    // Etchant surface coverage
    const auto &phi_e = globalData->getVectorData(2)[primID];

    const auto S_p = std::max(gamma_p * (1.- phi_ep * phi_p - phi_e ),0.);

    return std::pair<NumericType, rayTriple<NumericType>>{S_p, direction};
  }
  void initNew(rayRNG &RNG) override final {}
  int getRequiredLocalDataSize() const override final { return 1; }
  NumericType getSourceDistributionPower() const override final { return 1.; }
  std::vector<std::string> getLocalDataLabels() const override final {
    return std::vector<std::string>{"polymerRate"};
  }

private:
  static constexpr NumericType gamma_p = 0.26;
};

template <typename NumericType, int D> class SiO2Etching {
  psSmartPointer<psProcessModel<NumericType, D>> processModel = nullptr;

public:
  SiO2Etching(const double ionFlux, const double etchantFlux,
               const double polymerFlux, const NumericType ionEnergy = 100.,
               const int maskMaterial = 0) {
    processModel = psSmartPointer<psProcessModel<NumericType, D>>::New();

    // particles
    auto ion = std::make_unique<SiO2Ion<NumericType, D>>(
        ionEnergy);
    auto etchant = std::make_unique<SiO2Etchant<NumericType, D>>();
    auto polymer = std::make_unique<SiO2Polymer<NumericType, D>>();

    // surface model
    auto surfModel = psSmartPointer<SiO2SurfaceModel<NumericType, D>>::New(
        ionFlux, etchantFlux, polymerFlux);

    // velocity field
    auto velField =
        psSmartPointer<SiO2VelocityField<NumericType>>::New(maskMaterial);

    processModel->setSurfaceModel(surfModel);
    processModel->setVelocityField(velField);
    processModel->setProcessName("SiO2Etching");
    processModel->insertNextParticleType(ion);
    processModel->insertNextParticleType(etchant);
    processModel->insertNextParticleType(polymer);
  }

  psSmartPointer<psProcessModel<NumericType, D>> getProcessModel() {
    return processModel;
  }
};