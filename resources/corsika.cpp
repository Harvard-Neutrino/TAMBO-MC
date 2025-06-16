/*
 * (c) Copyright 2018 CORSIKA Project, corsika-project@lists.kit.edu
 *
 * This software is distributed under the terms of the GNU General Public
 * Licence version 3 (GPL Version 3). See file LICENSE for a full version of
 * the license.
 */

/* clang-format off */
// InteractionCounter used boost/histogram, which
// fails if boost/type_traits have been included before. Thus, we have
// to include it first...
#include <corsika/framework/process/InteractionCounter.hpp>
/* clang-format on */
#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/process/DynamicInteractionProcess.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>
#include <corsika/framework/utility/SaveBoostHistogram.hpp>

#include <corsika/modules/writers/EnergyLossWriter.hpp>
#include <corsika/modules/writers/LongitudinalWriter.hpp>
//  #include <corsika/modules/writers/PrimaryWriter.hpp>
#include <corsika/modules/writers/SubWriter.hpp>
#include <corsika/output/OutputManager.hpp>

#include <corsika/media/CORSIKA7Atmospheres.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/media/GeomagneticModel.hpp>
#include <corsika/media/GladstoneDaleRefractiveIndex.hpp>
#include <corsika/media/HomogeneousMedium.hpp>
#include <corsika/media/IMagneticFieldModel.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/MediumPropertyModel.hpp>
#include <corsika/media/NuclearComposition.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/UniformMagneticField.hpp>

#include <corsika/modules/BetheBlochPDG.hpp>
#include <corsika/modules/Epos.hpp>
#include <corsika/modules/LongitudinalProfile.hpp>
#include <corsika/modules/ObservationPlane.hpp>
#include <corsika/modules/PROPOSAL.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/Pythia8.hpp>
#include <corsika/modules/QGSJetII.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Sophia.hpp>
#include <corsika/modules/StackInspector.hpp>
#include <corsika/modules/thinning/EMThinning.hpp>

//added in 
#include <corsika/modules/tracking/TrackingStraight.hpp>

//for ICRC2023
#ifdef WITH_FLUKA
#include <corsika/modules/FLUKA.hpp>
#else
#include <corsika/modules/UrQMD.hpp>
#endif

#include <corsika/modules/UrQMD.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

//added in 
#include <corsika/stack/VectorStack.hpp>

#include <boost/filesystem.hpp>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>

#include <cstdlib>
#include <iomanip>
#include <limits>
#include <string>

using namespace corsika;
using namespace std;

using EnvironmentInterface =
    IRefractiveIndexModel<IMediumPropertyModel<IMagneticFieldModel<IMediumModel>>>;
using EnvType = Environment<EnvironmentInterface>;
using StackType = setup::Stack<EnvType>;

//using StackActive = setup::Stack<EnvType>;
using StackView = StackType::stack_view_type;

using TrackingType = setup::Tracking;
using Particle = StackType::particle_type;

long registerRandomStreams(long seed) {
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("qgsjet");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  RNGManager<>::getInstance().registerRandomStream("sophia");
  RNGManager<>::getInstance().registerRandomStream("epos");
  RNGManager<>::getInstance().registerRandomStream("pythia");
  RNGManager<>::getInstance().registerRandomStream("urqmd");
  RNGManager<>::getInstance().registerRandomStream("fluka");
  RNGManager<>::getInstance().registerRandomStream("proposal");
  RNGManager<>::getInstance().registerRandomStream("thinning");
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
    CORSIKA_LOG_INFO("random seed (auto) {}", seed);
  } else {
    CORSIKA_LOG_INFO("random seed {}", seed);
  }
  RNGManager<>::getInstance().setSeed(seed);
  return seed;
}

template <typename T>
using MyExtraEnv =
    GladstoneDaleRefractiveIndex<MediumPropertyModel<UniformMagneticField<T>>>;

int main(int argc, char** argv) {

  // the main command line description
  CLI::App app{"Simulate standard (upgoing) showers with CORSIKA 8."};

  // some options that we want to fill in
  int A, Z, nevent = 0;

  // the following section adds the options to the parser

  // we start by definining a sub-group for the primary ID
  auto opt_Z = app.add_option("-Z", Z, "Atomic number for primary")
                   ->check(CLI::Range(0, 26))
                   ->group("Primary");
  auto opt_A = app.add_option("-A", A, "Atomic mass number for primary")
                   ->needs(opt_Z)
                   ->check(CLI::Range(1, 58))
                   ->group("Primary");
  app.add_option("-p,--pdg", "PDG code for primary.")
      ->excludes(opt_A)
      ->excludes(opt_Z)
      ->group("Primary");
  // the remainding options
  app.add_option("-E,--energy", "Primary energy in GeV")
      ->required()
      ->check(CLI::PositiveNumber)
      ->group("Primary");
  app.add_option("-z,--zenith", "Primary zenith angle (deg)")
      ->default_val(0.)
      ->check(CLI::Range(0., 180.))
      ->group("Primary");
  app.add_option("-a,--azimuth", "Primary azimuth angle (deg)")
      ->default_val(0.)
      ->check(CLI::Range(0., 360.))
      ->group("Primary");
   app.add_option("--photoncut",
                 "Min. kin. energy of photons "
                 "in tracking (GeV)")
      ->default_val(0.5e-3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--emcut",
                 "Min. kin. energy of electrons and "
                 "positrons in tracking (GeV)")
      ->default_val(0.5e-3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--hadcut", "Min. kin. energy of hadrons in tracking (GeV)")
      ->default_val(0.3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--mucut", "Min. kin. energy of muons in tracking (GeV)")
      ->default_val(0.3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--cutoff-height",
                "Height above which no particles will be recorded (in km)")
      ->default_val(5.)
      ->check(CLI::Range(0.001,1.e2))
      ->group("Config");
  app.add_option("-O,--observation-height",
                 "Height above earth radius of the observation level (in km)")
      ->default_val(0.)
      ->check(CLI::Range(-1.e3, 1.e5))
      ->group("Config");
  app.add_option("-N,--nevent", nevent, "The number of events/showers to run.")
      ->default_val(1)
      ->check(CLI::PositiveNumber)
      ->group("Library/Output");
  app.add_option("-f,--filename", "Filename for output library.")
      ->required()
      ->default_val("corsika_library")
      ->check(CLI::NonexistentPath)
      ->group("Library/Output");
  app.add_option("-s,--seed", "The random number seed.")
      ->default_val(0)
      ->check(CLI::NonNegativeNumber)
      ->group("Misc.");
  bool force_interaction = false;
  app.add_flag("--force-interaction", force_interaction,
               "Force the location of the first interaction.")
      ->group("Misc.");
  app.add_option("-v,--verbosity", "Verbosity level: warn, info, debug, trace.")
      ->default_val("info")
      ->check(CLI::IsMember({"warn", "info", "debug", "trace"}))
      ->group("Misc.");
  app.add_option("-M,--hadronModel", "High-energy hadronic interaction model")
      ->default_val("SIBYLL-2.3d")
      ->check(CLI::IsMember({"SIBYLL-2.3d", "QGSJet-II.04", "EPOS-LHC"}))
      ->group("Misc.");
  app.add_option("-T,--hadronModelTransitionEnergy",
                 "Transition between high-/low-energy hadronic interaction "
                 "model in GeV")
      ->default_val(std::pow(10, 1.9)) // 79.4 GeV
      ->check(CLI::NonNegativeNumber)
      ->group("Misc.");
  app.add_option("--emthin",
                 "fraction of primary energy at which thinning of EM particles starts")
      ->default_val(1.e-6)
      ->check(CLI::Range(0., 1.))
      ->group("Thinning");
  app.add_option("--max-weight",
                 "maximum weight for thinning of EM particles (0 to select Kobal's "
                 "optimum times 0.5)")
      ->default_val(0)
      ->check(CLI::NonNegativeNumber)
      ->group("Thinning");
  bool multithin = false;
  app.add_flag("--multithin", multithin, "keep thinned particles (with weight=0)")
      ->group("Thinning");
  // Needed input for upgoing showers and non-horizontal observational plan 
  app.add_option("--xpos","x position of injection in TAMBO coordinates")
      ->default_val(0.)
      ->group("Misc.");
  app.add_option("--ypos","y position of injection in TAMBO coordinates")
      ->default_val(0.)
      ->group("Misc.");
  app.add_option("--zpos","z position of injection in TAMBO coordinates")
      ->default_val(0.)
      ->group("Misc.");
  app.add_option("--iteration","obs plane number")
      ->default_val(0.)
      ->group("Misc.");
  app.add_option("--xdir","X component of directional vector for observational level")
      ->default_val(0.)
      ->group("Misc.");
  app.add_option("--ydir","Y component of directional vector for observational level")
      ->default_val(0.)
      ->group("Misc.");
  app.add_option("--zdir","Z component of directional vector for observational level")
      ->default_val(-1.)
      ->group("Misc.");
  app.add_option("--x-intercept","X component of intercept of momentum vector with observational plane")
      ->default_val(0.)
      ->group("Misc.");
  app.add_option("--y-intercept","Y component of intercept of momentum vector with observational plane")
      ->default_val(0.)
      ->group("Misc.");
  app.add_option("--z-intercept","Z component of intercept of momentum vector with observational plane")
      ->default_val(-1.)
      ->group("Misc.");

  // parse the command line options into the variables
  CLI11_PARSE(app, argc, argv);

  if (app.count("--verbosity")) {
    auto const loglevel = app["--verbosity"]->as<std::string>();
    if (loglevel == "warn") {
      logging::set_level(logging::level::warn);
    } else if (loglevel == "info") {
      logging::set_level(logging::level::info);
    } else if (loglevel == "debug") {
      logging::set_level(logging::level::debug);
    } else if (loglevel == "trace") {
#ifndef _C8_DEBUG_
      CORSIKA_LOG_ERROR("trace log level requires a Debug build.");
      return 1;
#endif
      logging::set_level(logging::level::trace);
    }
  }

  // check that we got either PDG or A/Z
  // this can be done with option_groups but the ordering
  // gets all messed up
  if (app.count("--pdg") == 0) {
    if ((app.count("-A") == 0) || (app.count("-Z") == 0)) {
      CORSIKA_LOG_ERROR("If --pdg is not provided, then both -A and -Z are required.");
      return 1;
    }
  }

  // initialize random number sequence(s)
  auto seed = registerRandomStreams(app["--seed"]->as<long>());

  /* === START: SETUP ENVIRONMENT AND ROOT COORDINATE SYSTEM === */
  EnvType env;
  CoordinateSystemPtr const& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};
  Point const surface_{rootCS, 0_m, 0_m, constants::EarthRadius::Mean};
  GeomagneticModel wmm(center, corsika_data("GeoMag/WMM.COF"));

  // build a Linsley US Standard atmosphere into `env`,
  // https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm for Arequipa at 3km 
  create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv>(
      env, AtmosphereId::ColcaValley2022Jan, center, 1.000327, surface_, Medium::AirDry1Atm,
      MagneticFieldVector{rootCS, 22.7266_uT, -2.5322_nT, -4.2859_nT});

  /* === END: SETUP ENVIRONMENT AND ROOT COORDINATE SYSTEM === */

  ofstream atmout("earth.dat");
  for (LengthType h = 0_m; h < 110_km; h += 100_m) {
    Point const ptest{rootCS, 0_m, 0_m, constants::EarthRadius::Mean + h};
    auto rho =
        env.getUniverse()->getContainingNode(ptest)->getModelProperties().getMassDensity(
            ptest);
    atmout << h / 1_m << " " << rho / 1_kg * cube(1_m) << "\n";
  }
  atmout.close();

  /* === START: CONSTRUCT PRIMARY PARTICLE === */

  // parse the primary ID as a PDG or A/Z code
  Code beamCode;

  // check if we want to use a PDG code instead
  if (app.count("--pdg") > 0) {
    beamCode = convert_from_PDG(PDGCode(app["--pdg"]->as<int>()));
  } else {
    // check manually for proton and neutrons
    if ((A == 1) && (Z == 1))
      beamCode = Code::Proton;
    else if ((A == 1) && (Z == 0))
      beamCode = Code::Neutron;
    else
      beamCode = get_nucleus_code(A, Z);
  }
  HEPEnergyType mass = get_mass(beamCode);

  // particle energy
  HEPEnergyType const E0 = 1_GeV * app["--energy"]->as<double>();

  // direction of the shower in (theta, phi) space
  auto const thetaRad = app["--zenith"]->as<double>();
  auto const phiRad = app["--azimuth"]->as<double>();
  auto const iteration = app["--iteration"]->as<int>(); 

  // convert Elab to Plab
  HEPMomentumType P0 = calculate_momentum(E0, mass);

   //all in units of km 
  //z position is relative to origin 0,0,0 which is center of plane 
  //so if z = -0.624 that is -624 m below the observation plane 
  auto const xpos = app["--xpos"]->as<double>() * 1_km; 
  auto const ypos = app["--ypos"]->as<double>() * 1_km; 
  auto const zpos = app["--zpos"]->as<double>() * 1_km;
  auto const xintercept = app["--x-intercept"]->as<double>() * 1_km; 
  auto const yintercept = app["--y-intercept"]->as<double>() * 1_km;
  auto const zintercept = app["--z-intercept"]->as<double>() * 1_km;
  auto const obsHeight = app["--observation-height"]->as<double>();
  auto const cutHeight = app["--cutoff-height"]->as<double>();
    
  auto const particle_xdir = sin(thetaRad) * cos(phiRad);
  auto const particle_ydir = sin(thetaRad) * sin(phiRad);
  auto const particle_zdir = cos(thetaRad); 
 
  std::cout << "Direction Vector is: " << particle_xdir << " " << particle_ydir << " " << particle_zdir << std::endl;

  // convert the momentum to the zenith and azimuth angle of the primary
  auto const [px, py, pz] =
      std::make_tuple(P0 * sin(thetaRad) * cos(phiRad), P0 * sin(thetaRad) * sin(phiRad),
                      P0 * cos(thetaRad));
  auto plab = MomentumVector(rootCS, {px, py, pz});
  /* === END: CONSTRUCT PRIMARY PARTICLE === */

  /* === START: CONSTRUCT GEOMETRY === */

  /* === START: CONSTRUCT GEOMETRY === */
  auto const observationHeight = 1_km * obsHeight + constants::EarthRadius::Mean;
  auto const cutOffHeight = 1_km * cutHeight + constants::EarthRadius::Mean;

 //Point const surfaceCenter{rootCS, 0_m, 0_m, cutOffHeight};
 //auto skyLevel = EnvType::createNode<SeparationPlane>(Plane(surfaceCenter, DirectionVector(rootCS, {0, 0, 1})));

  std::cout << "This is the obs height: " << observationHeight << std::endl; 
  std::cout << "this is what it was originally: " << obsHeight << std::endl; 
  std::cout << "This is the Earth radius: " << constants::EarthRadius::Mean << std::endl;

  // Point const cutOffPoint{rootCS, 0_m, 0_m, cutOffHeight};
  Point const showerCore{rootCS, 0_m, 0_m, observationHeight};
  Point const injectionPos{rootCS,xpos,ypos,observationHeight+zpos};
  Point const interceptPos{rootCS,xintercept,yintercept,observationHeight+zintercept};

  // we make the axis much longer than the inj-core distance since the
  // profile will go beyond the core, depending on zenith angle
  ShowerAxis const showerAxis{injectionPos, (interceptPos - injectionPos) * 5.0, env};
  auto const dX = 10_g / square(1_cm); // Binning of the writers along the shower axis
  uint const nAxisBins = showerAxis.getMaximumX() / dX + 1; // Get maximum number of bins
  /* === END: CONSTRUCT GEOMETRY === */
  
  for (int i = 0; i < 10; ++i) {
  auto const d = i * 1_km;
  auto const test_showeraxis = showerAxis.getX(d); 
  std::cout << "Shower axis test getX(): " << test_showeraxis << std::endl;
  }
  double const emthinfrac = app["--emthin"]->as<double>();
  std::cout << "thinning fraction: " << emthinfrac << std::endl;
  double const maxWeight = std::invoke([&]() {
    if (auto const wm = app["--max-weight"]->as<double>(); wm > 0)
      return wm;
    else
      return 0.5 * emthinfrac * E0 / 1_GeV;
  });
  EMThinning thinning{emthinfrac * E0, maxWeight, !multithin};

  std::stringstream args;
  for (int i = 0; i < argc; ++i) { args << argv[i] << " "; }
  // create the output manager that we then register outputs with
  OutputManager output(app["--filename"]->as<std::string>(), seed, args.str());

  EnergyLossWriter dEdX{showerAxis, dX, nAxisBins};
  output.add("energyloss", dEdX);

  DynamicInteractionProcess<StackType> heModel;

  // have SIBYLL always for PROPOSAL photo-hadronic interactions
  auto sibyll = std::make_shared<corsika::sibyll::Interaction>(env);

 if (auto const modelStr = app["--hadronModel"]->as<std::string>();
      modelStr == "SIBYLL-2.3d") {
    heModel = DynamicInteractionProcess<StackType>{sibyll};
  } else if (modelStr == "QGSJet-II.04") {
    heModel = DynamicInteractionProcess<StackType>{
        std::make_shared<corsika::qgsjetII::Interaction>()};
  } else if (modelStr == "EPOS-LHC") {
    heModel = DynamicInteractionProcess<StackType>{
        std::make_shared<corsika::epos::Interaction>()};
  } else {
    CORSIKA_LOG_CRITICAL("invalid choice \"{}\"; also check argument parser", modelStr);
    return EXIT_FAILURE;
  }

  InteractionCounter heCounted{heModel};

  corsika::pythia8::Decay decayPythia;

  // hadronic photon interactions in resonance region
  corsika::sophia::InteractionModel sophia;

  HEPEnergyType const photoncut = 1_GeV * app["--photoncut"]->as<double>();
  HEPEnergyType const emcut = 1_GeV * app["--emcut"]->as<double>();
  HEPEnergyType const hadcut = 1_GeV * app["--hadcut"]->as<double>();
  HEPEnergyType const mucut = 1_GeV * app["--mucut"]->as<double>();

  std::cout << "Photon cut energy: " << photoncut << std::endl;
  std::cout << "EM cut energy: " << emcut << std::endl;
  std::cout << "Muon cut energy: " << mucut << std::endl;
  std::cout << "Hadron cut energy: " << hadcut << std::endl;

  ParticleCut<SubWriter<decltype(dEdX)>> cut(emcut, photoncut, hadcut, mucut, true, dEdX);

  // tell proposal that we are interested in all energy losses above the 
  set_energy_production_threshold(Code::Electron, emcut);
  set_energy_production_threshold(Code::Positron, emcut);
  set_energy_production_threshold(Code::Photon, photoncut);
  set_energy_production_threshold(Code::MuMinus, mucut);
  set_energy_production_threshold(Code::MuPlus, mucut);

  // energy threshold for high energy hadronic model. Affects LE/HE switch for
  // hadron interactions and the hadronic photon model in proposal
  HEPEnergyType const heHadronModelThreshold =
      1_GeV * app["--hadronModelTransitionEnergy"]->as<double>();

  corsika::proposal::Interaction emCascade(
      env, sophia, sibyll->getHadronInteractionModel(), heHadronModelThreshold);

  // use BetheBlochPDG for hadronic continuous losses, and proposal otherwise
  corsika::proposal::ContinuousProcess<SubWriter<decltype(dEdX)>> emContinuousProposal(
      env, dEdX);
  BetheBlochPDG<SubWriter<decltype(dEdX)>> emContinuousBethe{dEdX};
  struct EMHadronSwitch {
    EMHadronSwitch() = default;
    bool operator()(const Particle& p) const { return is_hadron(p.getPID()); }
  };
auto emContinuous =
      make_select(EMHadronSwitch(), emContinuousBethe, emContinuousProposal);

// LongitudinalWriter profile{showerAxis, nAxisBins, dX};
// output.add("profile", profile);
// LongitudinalProfile<SubWriter<decltype(profile)>> longprof{profile};

std::cout << "With fluka: " << WITH_FLUKA << std::endl; 


#ifdef WITH_FLUKA
  corsika::fluka::Interaction leIntModel{env};
#else
  corsika::urqmd::UrQMD leIntModel{};
#endif
  InteractionCounter leIntCounted{leIntModel};

// corsika::urqmd::UrQMD leIntModel{};
// InteractionCounter leIntCounted{leIntModel};


//orginally was 10000 for stackInspector
  StackInspector<StackType> stackInspect(10000, false, E0);
  //StackActive stack;

  // assemble all processes into an ordered process list
  struct EnergySwitch {
    HEPEnergyType cutE_;
    EnergySwitch(HEPEnergyType cutE)
        : cutE_(cutE) {}
    bool operator()(const Particle& p) const { return (p.getKineticEnergy() < cutE_); }
  };
  auto hadronSequence =
      make_select(EnergySwitch(heHadronModelThreshold), leIntCounted, heCounted);

  // observation plane
  
  auto const xdir = app["--xdir"]->as<double>();
  auto const ydir = app["--ydir"]->as<double>();
  auto const zdir = app["--zdir"]->as<double>();

  DirectionVector const propDir = DirectionVector(rootCS, {particle_xdir,particle_ydir,particle_zdir});

  Plane const obsPlane(showerCore, DirectionVector(rootCS, {xdir, ydir, zdir}));
  DirectionVector const leftVec = DirectionVector(rootCS,{0., particle_zdir/sqrt(pow(particle_ydir,2)+pow(particle_zdir,2)), -particle_ydir/sqrt(pow(particle_ydir,2)+pow(particle_zdir,2))});

//   CORSIKA_LOG_INFO("ShowerCore Center: {}", showerCore);
//   auto& obsPlaneFinal = obsPlanes.emplace_back(
//       obsPlane, DirectionVector(rootCS, {0, -zdir/sqrt(pow(ydir,2)+pow(zdir,2)), ydir/sqrt(pow(ydir,2)+pow(zdir,2))}), true, true);

  // * output module
  //OutputManager output(output_dir);
  // for (int i = 0; i < nPlane; i++) {
  //   output.add(fmt::format("particles_{:}", i), obsPlanes[i]);
  // }
  // // hard coded
  // auto obsPlaneSequence =
  //     make_sequence(obsPlanes[0], obsPlanes[1], obsPlanes[2], obsPlanes[3], obsPlanes[4]);


  // Plane const ceilingPlane(cutOffPoint, DirectionVector(rootCS, {0., 0., -1.})); 
  // ObservationPlane<TrackingType, ParticleWriterParquet> obsPlaneFirst{
  //     obsPlane, DirectionVector(rootCS, {0, -zdir/sqrt(pow(ydir,2)+pow(zdir,2)), ydir/sqrt(pow(ydir,2)+pow(zdir,2))}), true};
  // register ground particle output
  // ObservationPlane<setup::Tracking, ParticleWriterParquet> skyLevel{
  //     ceilingPlane, DirectionVector(rootCS,{1.,0.,0.})};
  //Direction Vector here is the x axis so y axis will be {0,-1,0} (with z being (0,0,-1)
  // output.add("observation_plane_particles", observationLevel);
  // output.add("500m_particles", observationLevel1);
  // output.add("1000m_particles", observationLevel2);
  // output.add("1500m_particles", observationLevel3);
  // output.add("2000m_particles", observationLevel4);
  // output.add("sky", skyLevel);

  // PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(observationLevel);
  // output.add("observation_plane", primaryWriter);
  // PrimaryWriter<TrackingType, ParticleWriterParquet> obs_500m_Writer(observationLevel1);
  // output.add("500m_plane", obs_500m_Writer);
  // PrimaryWriter<TrackingType, ParticleWriterParquet> obs_1000m_Writer(observationLevel2);
  // output.add("1000m_plane", obs_1000m_Writer);
  // PrimaryWriter<TrackingType, ParticleWriterParquet> obs_1500m_Writer(observationLevel3);
  // output.add("1500m_plane", obs_1500m_Writer);
  // PrimaryWriter<TrackingType, ParticleWriterParquet> obs_2000m_Writer(observationLevel4);
  // output.add("2000m_plane", obs_2000m_Writer);


  // DirectionVector const leftVec = DirectionVector(rootCS,{0., -particle_zdir/sqrt(pow(particle_ydir,2)+pow(particle_zdir,2)), particle_ydir/sqrt(pow(particle_ydir,2)+pow(particle_zdir,2))});

  // // * observation plane
  // std::vector<ObservationPlane<tracking_line::Tracking>> obsPlanes;

  ObservationPlane<TrackingType, ParticleWriterParquet> obsPlaneFirst{
       obsPlane, DirectionVector(rootCS, {0, -zdir/sqrt(pow(ydir,2)+pow(zdir,2)), ydir/sqrt(pow(ydir,2)+pow(zdir,2))}), true, true};

  output.add("particles_obs_final", obsPlaneFirst);

  // Point planeCenter{rootCS, xpos+(particle_xdir*iteration*500_m),ypos+(particle_ydir*iteration*500_m), observationHeight+zpos+(particle_zdir*iteration*500_m)};

  // std::cout << "PlaneCenter is: " << planeCenter << std::endl;

  // Plane const obsPlaneFinal(planeCenter, DirectionVector(rootCS, {-particle_xdir, -particle_ydir, -particle_zdir}));

  // ObservationPlane<TrackingType, ParticleWriterParquet> obsPlaneFirst{
  //      obsPlaneFinal, leftVec, true, true};

  // output.add(fmt::format("particles_{:}m", iteration*500), obsPlaneFirst);

  // PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(obsPlaneFinal);
  // output.add("observation_plane", primaryWriter);

  //  // * output module
  // //OutputManager output(output_dir);

  // for (int i = 0; i < nPlane-1; i++) {
  //   PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(obsPlanes[i]);
  //   output.add("observation_plane_{:}m", primaryWriter);

  //   output.add(fmt::format("particles_{:}m", (i + 1) * 500), obsPlanes[i]);
  // }

  // // hard coded
  // auto obsPlaneSequence =
  //     make_sequence(obsPlanes[0], obsPlanes[1], obsPlanes[2], obsPlanes[3], obsPlanes[4]);

  // auto obsPlaneSequence = 
  //   make_sequence(observationLevel1,observationLevel2,observationLevel3,observationLevel4,observationLevel);

  // PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(obsPlanes);
  // output.add("primary", primaryWriter);

  auto const showerCoreX_{showerCore.getCoordinates().getX()};
  auto const showerCoreY_{showerCore.getCoordinates().getY()};
  auto const injectionPosX_{injectionPos.getCoordinates().getX()};
  auto const injectionPosY_{injectionPos.getCoordinates().getY()};
  auto const injectionPosZ_{injectionPos.getCoordinates().getZ()};
  auto const triggerpoint_{Point(rootCS, injectionPosX_, injectionPosY_, injectionPosZ_)};
  std::cout << "Trigger Point is: " << triggerpoint_ << std::endl;
  auto const interceptpoint_{Point(rootCS,xintercept,yintercept,zintercept)};
  std::cout << "Intercept Point is: " << interceptpoint_ << std::endl; 

  // assemble the final process sequence with radio

  //stackInspector was here 

  auto obsPlaneSequence =
      make_sequence(obsPlaneFirst);

  // auto sequence =
  //     make_sequence(stackInspect, hadronSequence, decayPythia, emCascade, emContinuous,
  //             obsPlaneFinal, obsPlane1, thinning, cut);

   auto sequence =
      make_sequence(hadronSequence, decayPythia, emCascade, emContinuous,
              obsPlaneSequence, thinning, cut);



  /* === END: SETUP PROCESS LIST === */

  // create the cascade object using the default stack and tracking
  // implementation

  //tracking_line::Tracking tracking;
  tracking_line::Tracking tracking; 
  StackType stack;

  Cascade EAS(env, tracking, sequence, output, stack);

  // print our primary parameters all in one place
  if (app["--pdg"]->count() > 0) {
    CORSIKA_LOG_INFO("Primary PDG ID: {}", app["--pdg"]->as<int>());
  } else {
    CORSIKA_LOG_INFO("Primary Z/A: {}/{}", Z, A);
  }
  CORSIKA_LOG_INFO("Primary Energy: {}", E0);
  CORSIKA_LOG_INFO("Primary Momentum: {}", P0);
  CORSIKA_LOG_INFO("Primary Direction:  {}", plab.getNorm());
  CORSIKA_LOG_INFO("Point of Injection: {}", injectionPos.getCoordinates());
  CORSIKA_LOG_INFO("Shower Axis Length: {}", (interceptPos - injectionPos).getNorm() * 5.0);
  CORSIKA_LOG_INFO("Intercept Length: {}", (interceptPos - injectionPos).getNorm());

  // trigger the output manager to open the library for writing
  output.startOfLibrary();

  // loop over each shower
  for (int i_shower = 1; i_shower < nevent + 1; i_shower++) {

    CORSIKA_LOG_INFO("Shower {} / {} ", i_shower, nevent);

    // directory for outputs
    string const outdir(app["--filename"]->as<std::string>());
    string const labHist_file = outdir + "/inthist_lab_" + to_string(i_shower) + ".npz";
    string const cMSHist_file = outdir + "/inthist_cms_" + to_string(i_shower) + ".npz";

    // setup particle stack, and add primary particle
    stack.clear();

    // add the desired particle to the stack
    auto const primaryProperties = std::make_tuple(
        beamCode, calculate_kinetic_energy(plab.getNorm(), get_mass(beamCode)),
        plab.normalized(), injectionPos, 0_ns);
    stack.addParticle(primaryProperties);

    // if we want to fix the first location of the shower
    if (force_interaction) {
      CORSIKA_LOG_INFO("Fixing first interaction at injection point.");
      EAS.forceInteraction();
    }

    //primaryWriter.recordPrimary(primaryProperties);

    // run the shower
    EAS.run();

    // HEPEnergyType const Efinal =
    //     dEdX.getEnergyLost() + observationLevel.getEnergyGround();
     HEPEnergyType const Efinal =
          dEdX.getEnergyLost() + obsPlaneFirst.getEnergyGround();

    CORSIKA_LOG_INFO(
        "total energy budget (GeV): {} (dEdX={} ground={}), "
        "relative difference (%): {}",
        Efinal / 1_GeV, dEdX.getEnergyLost() / 1_GeV,
        obsPlaneFirst.getEnergyGround() / 1_GeV, (Efinal / E0 - 1) * 100);

    auto const hists = heCounted.getHistogram() + leIntCounted.getHistogram();

    save_hist(hists.labHist(), labHist_file, true);
    save_hist(hists.CMSHist(), cMSHist_file, true);
  }

  // and finalize the output on disk
  output.endOfLibrary();

  return EXIT_SUCCESS;
}
