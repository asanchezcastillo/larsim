////////////////////////////////////////////////////////////////////////
// Class:       PhotonHitConverter
// Plugin Type: producer (Unknown Unknown)
// File:        PhotonHitConverter_module.cc
//
// Generated at Tue Dec 13 05:30:41 2022 by Alejandro Castillo using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "nurandom/RandomUtils/NuRandomService.h"

#include "artg4tk/pluginDetectors/gdml/PhotonHit.hh"

// Random numbers
#include "CLHEP/Random/RandFlat.h"

#include <memory>
#include <string>
#include <vector>

namespace sim {
  class PhotonHitConverter;
}

class sim::PhotonHitConverter : public art::EDProducer {
public:
  explicit PhotonHitConverter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  PhotonHitConverter(PhotonHitConverter const&) = delete;
  PhotonHitConverter(PhotonHitConverter&&) = delete;
  PhotonHitConverter& operator=(PhotonHitConverter const&) = delete;
  PhotonHitConverter& operator=(PhotonHitConverter&&) = delete;
  // Required functions.
  void produce(art::Event& e) override;

private:
  // Declare member data here.
  const bool fUseLitePhotons;
  size_t nOpChannels;
  double fThickness;

  const std::vector<double> fAbsSpectrum;
  const std::vector<double> fAbsEnergy;
  const std::vector<double> fEmissionEnergy;
  const std::vector<double> fEmissionSpectrum;
  const std::vector<double> fQEEnergy;
  const std::vector<double> fQESpectrum;

  //PDS mapping
  opdet::sbndPDMapAlg pdsmap;

  //Random number generator
  CLHEP::HepRandomEngine& fScintTimeEngine;
  std::unique_ptr<CLHEP::RandFlat> fUniformGen;

  const double fScinPreScale;


  void initRand(CLHEP::HepRandomEngine&);
  bool isAbsorbed(double& , double& );
  double linearInterpolation(const std::vector<double>& , const std::vector<double>& , double& );
  double isEmitted(const std::vector<double>& , const std::vector<double>&);
  bool isDetected(const std::vector<double>& ,const std::vector<double>& , double& );
};

sim::PhotonHitConverter::PhotonHitConverter(fhicl::ParameterSet const& p)
  : EDProducer{p}, fUseLitePhotons(p.get<bool>("UseLitePhotons")), fThickness(p.get<double>("fThickness")), 
  fAbsSpectrum(p.get<std::vector<double>>("fAbsSpectrum")), fAbsEnergy(p.get<std::vector<double>>("fAbsEnergy")),
  fEmissionEnergy(p.get<std::vector<double>>("fEmissionEnergy")), 
  fEmissionSpectrum(p.get<std::vector<double>>("fEmissionSpectrum")), fQEEnergy(p.get<std::vector<double>>("fQEEnergy")), 
  fQESpectrum(p.get<std::vector<double>>("fQESpectrum")), 
  fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
        createEngine(0, "HepJamesRandom", "scinttime"),
        "HepJamesRandom",
        "scinttime",
        p,
        "SeedScintTime")), fScinPreScale(p.get<double>("fScintPreScale")) // ,
// More initializers here.
{
  // Call appropriate produces<>() functions here.
  if (fUseLitePhotons) {
    produces<std::vector<sim::SimPhotonsLite>>();
    produces<std::vector<sim::SimPhotonsLite>>("Reflected");
  }
  else {
    produces<std::vector<sim::SimPhotons>>();
    produces<std::vector<sim::SimPhotons>>("Reflected");
  }

  this->initRand(fScintTimeEngine);

  if (fAbsSpectrum.size() != fAbsEnergy.size()) {
  throw art::Exception(art::errors::Configuration)
    << "Error with TPB absorption spectrum, check the lenght of the data vector."
    << "\n";
  }

}

void sim::PhotonHitConverter::produce(art::Event& e)
{
  geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
  nOpChannels = geom.NOpDets();

  //SimPhotonsLite
  std::unique_ptr<std::vector<sim::SimPhotonsLite>> photLiteCol{
    new std::vector<sim::SimPhotonsLite>{}};
  std::unique_ptr<std::vector<sim::SimPhotonsLite>> photLiteCol_ref{
    new std::vector<sim::SimPhotonsLite>{}};
  auto& photonLiteCollection{*photLiteCol};
  auto& photonLiteCollection_ref{*photLiteCol_ref};
  //SimPhotons
  std::unique_ptr<std::vector<sim::SimPhotons>> photCol{new std::vector<sim::SimPhotons>{}};
  std::unique_ptr<std::vector<sim::SimPhotons>> photCol_ref{new std::vector<sim::SimPhotons>{}};
  auto& photonCollection{*photCol};
  auto& photonCollection_ref{*photCol_ref};

  if (fUseLitePhotons) { //SimPhotonsLite
    photonLiteCollection.resize(nOpChannels);
    photonLiteCollection_ref.resize(nOpChannels);
    for (unsigned int i = 0; i < nOpChannels; ++i) {
      photonLiteCollection[i].OpChannel = i;
      photonLiteCollection_ref[i].OpChannel = i;      
    }
  }
  else { //SimPhotons
    photonCollection.resize(nOpChannels);
    photonCollection_ref.resize(nOpChannels);
    for (unsigned int i = 0; i < nOpChannels; ++i) {
      photonCollection[i].fOpChannel = i;
      photonCollection_ref[i].fOpChannel = i;
    }
  }

  int reach_pmt_coated=0;
  int reach_pmt_uncoated=0;
  int reach_xarapuca=0;
  int detected_pmt_coated=0;
  int detected_pmt_uncoated=0;

  typedef std::vector<art::Handle<artg4tk::PhotonHitCollection>> HandleVector;
  auto allSims = e.getMany<artg4tk::PhotonHitCollection>();
  for (HandleVector::const_iterator i = allSims.begin(); i != allSims.end(); ++i) {
    const artg4tk::PhotonHitCollection& sims(**i);
    for (artg4tk::PhotonHitCollection::const_iterator j = sims.begin(); j != sims.end(); ++j) {
      const artg4tk::PhotonHit& photonHit = *j;
        if (fUseLitePhotons) { 
          auto time = static_cast<int>(photonHit.GetTime());
          auto channel = static_cast<unsigned int>(photonHit.GetID()); 
          std::string pd_type = pdsmap.pdType(channel);
          auto energy = photonHit.GetEdep()/CLHEP::eV;
          int processID = photonHit.GetProcessID();
          if(processID==0) //If is generated by Cerenkov
          {
            if(fScinPreScale<fUniformGen->fire()) continue; // ScintPreScale which is applied for scintillation ligth
            bool detected = isDetected(fQEEnergy, fQESpectrum, energy); // Get if the photon is detected. 
            if( pd_type=="pmt_coated")
            {
              reach_pmt_coated+=1;
              double fAbsLength = linearInterpolation(fAbsEnergy, fAbsSpectrum, energy);
              bool Absorbed = isAbsorbed(fThickness, fAbsLength); // The photons is absorbed with a probability depending on its energy and thickness of the TPB. 
              if(Absorbed) // If the photon is absorbed, then reemit it
              {
                energy = isEmitted(fEmissionEnergy, fEmissionSpectrum); // Change the energy of the photon to the reemited energy
                bool detected = isDetected(fQEEnergy, fQESpectrum, energy); // Check if the photon is detected with the WLSed energy.
                if(detected)
                {
                  detected_pmt_coated+=1;
                  ++photonLiteCollection[channel].DetectedPhotons[time];
                }
              }
              else // If the photon is not absorbed check if it's visible and if it's detected is true
              {
                if(detected)
                {
                  detected_pmt_coated+=1;
                  ++photonLiteCollection_ref[channel].DetectedPhotons[time];
                }       
              } 
            } // End PMT coated
            else if( pd_type=="pmt_uncoated"){
              reach_pmt_uncoated+=1;
              if(detected)
              {
                ++photonLiteCollection_ref[channel].DetectedPhotons[time];
                detected_pmt_uncoated+=1;
              } 
              else continue;
            } // End PMT uncoated
            else if(pd_type=="xarapuca_vis" || pd_type=="xarapuca_vuv"){ // If !=pmt_coated and !=pmt_uncoated the it is XArapuca. Leave it as it was.
              reach_xarapuca+=1;
              if (photonHit.GetEdep() > 6.19 * CLHEP::eV) {
                ++photonLiteCollection[channel].DetectedPhotons[time];
              }
            } // End if XArapuca
          } // End If Cerenkov
          else
          {
            if (photonHit.GetEdep() > 6.19 * CLHEP::eV) {
              ++photonLiteCollection[channel].DetectedPhotons[time];
            }
            else {
              ++photonLiteCollection_ref[channel].DetectedPhotons[time];
            }
          } 
        } // End SimPhotonsLite
        else {
          sim::OnePhoton photon;
          photon.SetInSD = false;
          photon.Energy = photonHit.GetEdep()/CLHEP::eV;
          auto time = photonHit.GetTime();
          photon.Time = time;
          auto channel = static_cast<unsigned int>(photonHit.GetID());
          std::string pd_type = pdsmap.pdType(channel);
          auto energy = photonHit.GetEdep()/CLHEP::eV;
          int processID = photonHit.GetProcessID();
          if(processID==0)
          {
            if(fScinPreScale<fUniformGen->fire()) continue; // ScintPreScale which is applied for scintillation ligth
            bool detected = isDetected(fQEEnergy, fQESpectrum, energy); // Get if the photon is detected.
            if( pd_type=="pmt_coated")
            {
              reach_pmt_coated+=1;
              double fAbsLength = linearInterpolation(fAbsEnergy, fAbsSpectrum, energy);
              bool Absorbed = isAbsorbed(fThickness, fAbsLength); // The photons is absorbed with a probability depending on its energy and thickness of the TPB. 
              if(Absorbed) // If the photon is absorbed, then reemit it
              {
                energy = isEmitted(fEmissionEnergy, fEmissionSpectrum); // Change the energy of the photon to the reemited energy
                bool detected = isDetected(fQEEnergy, fQESpectrum, energy); // Check if the photon is detected with the WLSed energy.
                if(detected)
                {
                  photonCollection[channel].insert(photonCollection[channel].end(), 1, photon);
                  detected_pmt_coated+=1;
                }
              } // End if absorbed
              else // If the photon is not absorbed check if it's visible and it it's detected is true
              {
                if(detected)
                {
                  photonCollection_ref[channel].insert(photonCollection[channel].end(), 1, photon);
                }       
              } 
            }
            else if( pd_type=="pmt_uncoated"){
              reach_pmt_uncoated+=1;
              if(detected)
              {
                photonCollection_ref[channel].insert(photonCollection[channel].end(), 1, photon);
                detected_pmt_uncoated+=1;
              } 
              else continue;
            } // End PMT uncoated 
            else if(pd_type=="xarapuca_vis" || pd_type=="xarapuca_vuv"){ // If !=pmt_coated and !=pmt_uncoated the it is XArapuca. Leave it as it was.
              reach_xarapuca+=1;
              if (photonHit.GetEdep() > 6.19 * CLHEP::eV) {
                  photonCollection[channel].insert(photonCollection[channel].end(), 1, photon);
              }
              else {
                  photonCollection_ref[channel].insert(photonCollection[channel].end(), 1, photon);
              }  
            } // End if XArapuca
          } // End if Cerenkov
          else
          {
            if (photon.Energy > 6.19 * CLHEP::eV) {
              photonCollection[channel].insert(photonCollection[channel].end(), 1, photon);
            }
            else {
              photonCollection_ref[channel].insert(photonCollection_ref[channel].end(), 1, photon);
            }
          }
        } // End SimPhotons   
    }
  } 

  if (fUseLitePhotons) {
    e.put(std::move(photLiteCol));
    e.put(std::move(photLiteCol_ref), "Reflected");
  }
  else {
    e.put(std::move(photCol));
    e.put(std::move(photCol_ref), "Reflected");
  }
}

bool sim::PhotonHitConverter::isAbsorbed(double& fThickness, double& fAbsLength)
{
  double AbsProbability= 1 - std::exp(-fAbsLength/fThickness);
  bool Absorbed = (AbsProbability<fUniformGen->fire());
  return Absorbed;
}

double sim::PhotonHitConverter::isEmitted(const std::vector<double>& fEmissionEnergy, const std::vector<double>& fEmissionSpectrum)
{
  double sum = 0.0;
  for (const double& probability : fEmissionSpectrum) {
      sum += probability;
  } 
  double randomValue = fUniformGen->fire();
  double cumulativeProbability = 0.0;
  for (size_t i = 0; i < fEmissionEnergy.size(); i++) {
      cumulativeProbability += fEmissionSpectrum.at(i) / sum;
      if (randomValue <= cumulativeProbability) {
          return fEmissionEnergy.at(i);
      }
  }
  return 0;
}

bool sim::PhotonHitConverter::isDetected(const std::vector<double>& fQEEnergy,const std::vector<double>& fQESpectrum, double& energy)
{
  double Probability = linearInterpolation(fQEEnergy, fQESpectrum, energy);
  bool detected = (Probability<fUniformGen->fire());
  return detected;
}

void sim::PhotonHitConverter::initRand(CLHEP::HepRandomEngine& engine)
{
  fUniformGen = std::make_unique<CLHEP::RandFlat>(engine);
}

double sim::PhotonHitConverter::linearInterpolation(const std::vector<double>& X, const std::vector<double>& Y, double& A)
{
    // Find the index where A falls in the X vector
    size_t index = 0;
    while (X[index] < A) {
        index++;
    }
    // Perform linear interpolation
    if (index == 0) {
        // A is smaller than the first value in X, so we can't interpolate
        return Y[0];
    } else {
        double x0 = X[index - 1];
        double x1 = X[index];
        double y0 = Y[index - 1];
        double y1 = Y[index];
        // Linear interpolation formula
        return y0 + (A - x0) * (y1 - y0) / (x1 - x0);
    }
}

DEFINE_ART_MODULE(sim::PhotonHitConverter)
