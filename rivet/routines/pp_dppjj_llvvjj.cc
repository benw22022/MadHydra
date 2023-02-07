//=================================  Introduction  =================================//<
// -*- C++ -*-
//Mohamed Aly and Ahmed Abdelmotteleb
//Majorana Neutrino Sensitivity at ATLAS
//Rivet analysis code that runs over our data [signal + background] and applies cuts

//=================================  Libraries  =================================//
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <map>
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedMET.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "TFile.h"
#include "TTree.h"

#include <math.h>

//==========================================================================//
//=================================  MAIN  =================================//
//==========================================================================//

namespace Rivet {

  struct ParticleContainer{

    FourMomentum momentum; 
    long   pid;
    bool   from_tau;
    Particles parents;

    ParticleContainer(){}

    ParticleContainer(const Particle& p){
      momentum = p.momentum();
      from_tau = p.fromTau();
      parents = p.parents();
    }

    ParticleContainer(const Jet& j){
      momentum = j.momentum();
      from_tau = j.tauTagged();
    
    }

    ParticleContainer& operator+(const ParticleContainer& other){

      if (this == &other){  return *this;}

      this->momentum = this->momentum + other.momentum;

      return *this;
    }

  };
  

  //divided into constructor + 3 loops: init, analyze and finalize

  //==================================================================================//
  //=================================  Constructer ==================================//
  //=================================================================================//
  class pp_dpp_llvvjj : public Analysis {

    public:

      double calc_deltaR(const double& y1, const double& y2, const double& phi1, const double& phi2){

      // double y1 = p1.momentum.rapidity();
      // double y2 = p2.momentum.rapidity();
      // double phi1 = p1.momentum.phi();
      // double phi2 = p2.momentum.phi();

      return sqrt(pow((y1 - y2), 2) + pow((phi1 - phi2), 2));

      };

      /// Constructor
      pp_dpp_llvvjj() : Analysis("pp_dpp_llvvjj") {   }

      //=================================  Declaring global variables  =================================//

      double sumW = 0;  //sum of weights
      double var_sumW = 0; //variance of sum of weights
      double fid_sumW = 0; //fiducial sum of weights
      double var_fid_sumW = 0; //variance of fiducial sum of weights

      //Event counter
      int event_number = 0;
      int survive_event_number = 0;
      int num_events_before_a_cut=0;
      int num_events_after_a_cut=0;
      //=================================  Declaring global functions  =================================//

      //====================================================================================//
      //=================================  Initialisation  =================================//
      //====================================================================================//

      void init() {
        
        // 0 == no object smearing, 1 == ATLAS Run2 smearing, 2 == CMS Run2 smearing
        SMEAR = 0;
        
        // 0 == no smearing, 1 == smear
        MET_SMEAR = 1;
        JET_SMEAR = 1;
        MUON_SMEAR = 1;

        // Output file
        OutputFile = new TFile("output.root", "RECREATE");
        OutputTree = new TTree("tree", "tree");
        MetaDataTree = new TTree("metadata", "metadata");

        OutputTree->Branch("jet_pt",     &jet_pt,     "jet_pt/F");
        OutputTree->Branch("jet_rap",    &jet_rap,    "jet_rap/F");
        OutputTree->Branch("jet_phi",    &jet_pt,     "jet_phi/F");
        OutputTree->Branch("jet_mass",   &jet_mass,   "jet_mass/F");
        OutputTree->Branch("jet_energy", &jet_energy, "jet_energy/F");
        OutputTree->Branch("jet_px",     &jet_px,     "jet_px/F");
        OutputTree->Branch("jet_py",     &jet_py,     "jet_py/F");
        OutputTree->Branch("jet_pz",     &jet_pz,     "jet_pz/F");
        OutputTree->Branch("jet_pid",    &jet_pid,    "jet_pid/L");

        OutputTree->Branch("lep_pt",     &lep_pt,     "lep_pt/F");
        OutputTree->Branch("lep_rap",    &lep_rap,    "lep_rap/F");
        OutputTree->Branch("lep_phi",    &lep_pt,     "lep_phi/F");
        OutputTree->Branch("lep_mass",   &lep_mass,   "lep_mass/F");
        OutputTree->Branch("lep_energy", &lep_energy, "lep_energy/F");
        OutputTree->Branch("lep_px",     &lep_px,     "lep_px/F");
        OutputTree->Branch("lep_py",     &lep_py,     "lep_py/F");
        OutputTree->Branch("lep_pz",     &lep_pz,     "lep_pz/F");
        OutputTree->Branch("lep_pid",    &lep_pid,    "lep_pid/L");

        OutputTree->Branch("tau_pt",     &tau_pt,     "tau_pt/F");
        OutputTree->Branch("tau_rap",    &tau_rap,    "tau_rap/F");
        OutputTree->Branch("tau_phi",    &tau_pt,     "tau_phi/F");
        OutputTree->Branch("tau_mass",   &tau_mass,   "tau_mass/F");
        OutputTree->Branch("tau_energy", &tau_energy, "tau_energy/F");
        OutputTree->Branch("tau_px",     &tau_px,     "tau_px/F");
        OutputTree->Branch("tau_py",     &tau_py,     "tau_py/F");
        OutputTree->Branch("tau_pz",     &tau_pz,     "tau_pz/F");
        OutputTree->Branch("tau_pid",    &tau_pid,    "tau_pid/L");

        OutputTree->Branch("tau_comp_pt",     &tau_comp_pt,     "tau_comp_pt/F");
        OutputTree->Branch("tau_comp_rap",    &tau_comp_rap,    "tau_comp_rap/F");
        OutputTree->Branch("tau_comp_phi",    &tau_comp_pt,     "tau_comp_phi/F");
        OutputTree->Branch("tau_comp_mass",   &tau_comp_mass,   "tau_comp_mass/F");
        OutputTree->Branch("tau_comp_energy", &tau_comp_energy, "tau_comp_energy/F");
        OutputTree->Branch("tau_comp_px",     &tau_comp_px,     "tau_comp_px/F");
        OutputTree->Branch("tau_comp_py",     &tau_comp_py,     "tau_comp_py/F");
        OutputTree->Branch("tau_comp_pz",     &tau_comp_pz,     "tau_comp_pz/F");
        OutputTree->Branch("tau_comp_pid",    &tau_comp_pid,    "tau_comp_pid/L");

        OutputTree->Branch("MET_pt",     &MET_pt,  "MET_pt/F");
        OutputTree->Branch("MET_rap",    &MET_rap, "MET_rap/F");
        OutputTree->Branch("MET_phi",    &MET_pt,  "MET_phi/F");
        OutputTree->Branch("MET",        &MET,     "MET/F");
        OutputTree->Branch("MET_px",     &MET_px,  "MET_px/F");
        OutputTree->Branch("MET_py",     &MET_py,  "MET_py/F");
        OutputTree->Branch("MET_pz",     &MET_pz,  "MET_pz/F");

        OutputTree->Branch("ATLAS_smeared_MET_pt",     &ATLAS_smeared_MET_pt,  "ATLAS_smeared_MET_pt/F");
        OutputTree->Branch("ATLAS_smeared_MET_rap",    &ATLAS_smeared_MET_rap, "ATLAS_smeared_MET_rap/F");
        OutputTree->Branch("ATLAS_smeared_MET_phi",    &ATLAS_smeared_MET_pt,  "ATLAS_smeared_MET_phi/F");
        OutputTree->Branch("ATLAS_smeared_MET",        &ATLAS_smeared_MET,     "ATLAS_smeared_MET/F");
        OutputTree->Branch("ATLAS_smeared_MET_px",     &ATLAS_smeared_MET_px,  "ATLAS_smeared_MET_px/F");
        OutputTree->Branch("ATLAS_smeared_MET_py",     &ATLAS_smeared_MET_py,  "ATLAS_smeared_MET_py/F");
        OutputTree->Branch("ATLAS_smeared_MET_pz",     &ATLAS_smeared_MET_pz,  "ATLAS_smeared_MET_pz/F");

        MetaDataTree->Branch("total_xs",               &xs_total,              "total_xs_madgraph/F");
        MetaDataTree->Branch("total_xs_error",         &xs_total_err,          "total_xs_error/F");
        MetaDataTree->Branch("sum_of_weights",         &sum_of_weights,         "sum_of_weights/F");
        MetaDataTree->Branch("sum_of_weights_err",     &sum_of_weights_err,     "sum_of_weights_err/F");
        MetaDataTree->Branch("fid_sum_of_weights",     &fid_sum_of_weights,     "fid_sum_of_weights/F");
        MetaDataTree->Branch("fid_sum_of_weights_err", &fid_sum_of_weights_err, "fid_sum_of_weights_err/F");
        MetaDataTree->Branch("fid_xs",                 &xs_fid,                 "fid_xs/F");
        MetaDataTree->Branch("fid_xs_err",             &fid_xs_err,             "fid_xs_err/F");
        MetaDataTree->Branch("nsurvive_events",        &nsurvive_events,        "nsurvive_events/F");
        MetaDataTree->Branch("ninit_events",           &ninit_events,           "ninit_events/F");
        MetaDataTree->Branch("presel_eff",             &presel_eff,             "preself_eff/F");



        //=================================  General Final States Declaration  =================================//
        // define final state
        const FinalState fs(Cuts::abseta < 5);

        // define visible final states (not pT miss)
        VisibleFinalState vfs(Cuts::abseta < 4.5);
        declare(vfs,"vfs");

        // Define usable particles
        declare(UnstableParticles(), "UFS");

        //=================================  MET Declaration  =================================//
        MissingMomentum missing_et(fs);
        declare(missing_et,"missing_et");

        // define MET smeared with ATLAS Run2
        SmearedMET met_smeared_atlas(missing_et, MET_SMEAR_ATLAS_RUN2);
        declare(met_smeared_atlas,"met_smeared_atlas");

        // define MET smeared with CMS Run2
        SmearedMET met_smeared_cms(missing_et, MET_SMEAR_CMS_RUN2);
        declare(met_smeared_cms,"met_smeared_cms");

        //=================================  Photon and Neutrino Decleration  =================================//
        // define photon
        IdentifiedFinalState photon(fs);
        photon.acceptIdPair(PID::PHOTON); 
        declare(photon,"photon");

        // identify the Neutrinos
        IdentifiedFinalState neutrinos(fs);
        neutrinos.acceptNeutrinos();
        declare(neutrinos, "neutrinos");

        //=================================  Muons Decleration  =================================//
        // all Muons
        IdentifiedFinalState muon(fs);
        muon.acceptIdPair(PID::MUON);
        declare(muon, "muon");

        // muons that are visible in a detector 
        VetoedFinalState vfs_muons(vfs);
        vfs_muons.addVetoPair(PID::MUON, Cuts::abseta > 2.5);
        declare(vfs_muons,"vfs_muons");

        // define dressed muons as muons observed with a photon at radius 0.1 at  most from them
        // DressedLeptons dressed_muons(photon, muon, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
        DressedLeptons dressed_muons(photon, muon, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 5*GeV);
        declare(dressed_muons, "dressed_muons");

        SmearedParticles muons_eff_smear_atlas(dressed_muons, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2);
        declare(muons_eff_smear_atlas, "muons_eff_smear_atlas");

        //=================================  Electrons Decleration  =================================//
        // all electrons
        IdentifiedFinalState elec(fs);
        elec.acceptIdPair(PID::ELECTRON);
        declare(elec, "elec");

        // electrons that are visible in a detector       
        VetoedFinalState vfs_elecs(vfs);
        vfs_elecs.addVetoPair(PID::ELECTRON,Cuts::abseta > 2.5);
        declare(vfs_elecs,"vfs_elecs");

        //define dressed electrons as muons observed with a photon at radius 0.1 at  most from them
        DressedLeptons dressed_elecs(photon, elec, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
        declare(dressed_elecs, "dressed_elecs");

        // SmearedParticles elecs_eff_smear_atlas(dressed_elecs, ELECTRON_EFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2);
        // declare(elecs_eff_smear_atlas, "elecs_eff_smear_atlas");

        //=================================  Tau Decleration  =================================//
        

        //=================================  Jets Decleration  =================================//
        
        // use FastJet to select Jets 
        FastJets fj(fs, FastJets::ANTIKT, 0.4);
        declare(fj,"jets");

        // define jets smeared with ATLAS Run2
        SmearedJets jets_smeared_atlas(fj,JET_SMEAR_ATLAS_RUN2);
        declare(jets_smeared_atlas,"jets_smeared_atlas");
      }
      //===============================================================================//
      //=================================  Analysis  =================================//
      //==============================================================================//

      /// Perform the per-event analysis
      //select particles, filter according to cuts, loop over combinations, construct observables, fill histograms.
      // This is where the per-event aspect of the analysis algorithm goes.

      ///Do the event by event analysis here (loop over every event)
      void analyze(const Event& event) {
        
        event_number++;

        if(event_number % 1000 == 0){std::cout << "Proccessing event: " << event_number << std::endl;}

        //================================= VISIBLE FINAL STATE PARTICLES IN EVENT  =================================//
        
        // get visible final state particles: veto events with < 4 vfs particles
        Particles vfs_particles = apply<VisibleFinalState>(event, "vfs").particles();
        if (vfs_particles.size()  < 4) { vetoEvent; }

        // array of leptons after applying second cut (first cut is below-->cand_electrons & cand_muons)
        Particles recon_leptons;

        // array of jets after applying first cut (pT/eta)
        Jets cand_jets;

        // array of jets after applying removing hadronic taus
        Jets nontau_jets;

        // array of jets after applying second cut (lep iso)
        Jets recon_jets;

        //define a 4-momentum array for missing pT
        FourMomentum pTmiss;
        FourMomentum pTmiss_muons;
        FourMomentum pTmiss_elecs;

        //define gap jets
        FourMomentum gapjet(0., 0., 0., 0.);

        //=================== weight and variance ===================//

        //This will calculate the weight for each of the N events that rivet processes

        //weight of each event (Rivet has for some reason removed the weight() method, the 0th element should be the default weight from looking at the output of checkMetaSG.py)
        const double weight = event.weights()[0];

        //summing over all the weights (sum of weights)
        sumW += weight;

        //Sum up the weights^2 --> This happens to be equal to the variance
        var_sumW += std::pow(weight, 2);

        //================================= Building Candidate Leptons  -- Dressed =================================//
        
        //=================== Muons ===================//
        const Particles cand_muons = apply<DressedLeptons>(event, "dressed_muons").particlesByPt(Cuts::pT > 5*GeV);
        const Particles& smeared_muons_atlas = apply<ParticleFinder>(event, "muons_eff_smear_atlas").particlesByPt();       
        
        //=================== Electrons ===================//
        const Particles cand_elecs = apply<DressedLeptons>(event, "dressed_elecs").particlesByPt( Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
        const Particles& smeared_elecs_atlas = apply<ParticleFinder>(event, "elecs_eff_smear_atlas").particlesByPt();
        
        //=================== All Leptons -- Dressed ===================//
        Particles cand_leptons;

        for (unsigned int i{0}; i < cand_muons.size(); i++) { cand_leptons.push_back(cand_muons[i]); }
        for (unsigned int i{0}; i < cand_elecs.size(); i++) { cand_leptons.push_back(cand_elecs[i]); }

        //=================== Building Candidate Neutrinos and Hadronic Neutrino Flag  ===================//

        //=================== Neutrinos ===================//
        Particles neutrinos = applyProjection<IdentifiedFinalState>(event, "neutrinos").particlesByPt();
        
        //=================== Building Candidate Jets ===================//

        // create FastJets array with FastJets projection conditions defined earlier
        const Jets& truth_jets = apply<FastJets>(event, "jets").jetsByPt( (Cuts::abseta < 4.5) && (Cuts::pT > 20*GeV) );
        const Jets& jets_smeared_atlas = apply<JetAlg>(event, "jets_smeared_atlas").jetsByPt( (Cuts::abseta < 4.5) && (Cuts::pT > 20*GeV) );

        // defining candidate Jets based on smearing setting
        if (SMEAR == 0) 
        { 
          for (unsigned int i{0}; i < truth_jets.size(); i++)  { cand_jets.push_back(truth_jets[i]); }
        }
        if (SMEAR == 1 && JET_SMEAR) 
        { 
          for (unsigned int i{0}; i < jets_smeared_atlas.size(); i++)  { cand_jets.push_back(jets_smeared_atlas[i]); }
        }
          
        // need at least two jets
        if(cand_jets.size() < 2){ vetoEvent; }

        //=================== Choosing non-hadronic Leptons ===================//

        // veto leptons from hadrons      
        Particles non_hadronic_leps;
        for (unsigned int i{0}; i < cand_leptons.size(); i++) 
        {
          bool away_from_jet = true;  
          if (cand_leptons[i].fromHadron()) {
              away_from_jet = false;
            }
          if ( away_from_jet ) { non_hadronic_leps.push_back(cand_leptons[i]); }
        }    

        if (non_hadronic_leps.size() < 2) { vetoEvent; }

      //=================== Choosing reconstructed jets ===================//
      
      // use non-hadronic leptons to select jets that are far from prompt leptons 
      for (unsigned int i{0}; i < cand_jets.size(); i++)
      {
        for (unsigned int j{0}; j < non_hadronic_leps.size(); j++){

          auto p1     = non_hadronic_leps[j].momentum();
          auto p2     = cand_jets[i].momentum();
          double y1   = p1.rapidity();
          double y2   = p2.rapidity();
          double phi1 = p1.phi();
          double phi2 = p2.phi();
          double delR = calc_deltaR(y1, y2, phi1, phi2);

          if (delR > 0.3) { recon_jets.push_back(cand_jets[i]); }
        }
      }

      // need at least two jets
      if (recon_jets.size() < 2){ vetoEvent;}


      //=================== work out tau vis ===================//

      Particles non_tau_leps;
      ParticleContainer tau_vis;
      std::vector<std::vector<ParticleContainer>> tau_vis_comp;
      std::vector<ParticleContainer> tau_vis_tmp; // unsorted bag of particles from taus, will use to work out which particle goes with which tau if ditau
      
      // get all particles from tau decays, place rest into non-tau containers
      for(const Jet &j : cand_jets){ 
        if (j.tauTagged()){ tau_vis_tmp.push_back(ParticleContainer(j)); }
        else{nontau_jets.push_back(j);}
      }

      for(const Particle &p : non_hadronic_leps){ 
        if (p.fromTau()){ tau_vis_tmp.push_back(ParticleContainer(p)); }
        else{non_tau_leps.push_back(p);}
      }

      // get set of mother ids
      std::set<Particle> parents;
      for(const ParticleContainer &pc : tau_vis_tmp){ mother_ids.insert(pc.parents[0]); }

      // get all particles orginating from each mother tau
      for(unsigned int i{0}; i<parents.size(); i++){
        std::vector<ParticleContainer> tmp_tau_comp;
        for(const ParticleContainer &pc : tau_vis_tmp){
          if(tau_vis_tmp.parents[0] == parents[i]){tmp_tau_comp.push_back(pc); }
        }
        tau_vis_comp.push_back(ParticleContainer(tmp_tau_comp));
      }

      // sum 4-vecs of tau componants to get combined tau
      for(const ParticleContainer &tau_particles : tau_vis_tmp){
        ParticleContainer comb_tau = tau_particles[0];
        for(int i{1}; i<tau_particles.size(); i++){
          comb_tau += tau_particles[i]
        }
        tau_vis.push_back(tau_comb);
      }

      //=================== Choosing and sorting reconstructed Leptons ===================//
      
      // use non hadronic leptons to select leptons that are far enough from all other leptons
      for (unsigned int i{0}; i < non_tau_leps.size(); i++)
      {
          for (unsigned int j{0}; j < non_tau_leps.size(); j++)
          {

            double y1   = non_tau_leps[j].momentum.rapidity;
            double y2   = non_tau_leps[i].momentum.rapidity;
            double phi1 = non_tau_leps[j].momentum.phi;
            double phi2 = non_tau_leps[i].momentum.phi;
            double delR = calc_deltaR(y1, y2, phi1, phi2);

            if(i != j && delR < 0.2) { recon_leptons.push_back(non_tau_leps[i]); }
          }
        }

      // // Combine reco leps with taus (if any)
      // for(const Particle &p : tau_vis){
      //   recon_leps.push_back(p)
      // }

      //Sort the leptons by pT
      std::sort (recon_leptons.begin(), recon_leptons.end(), [](Particle const &a, Particle const &b) { return a.momentum().pT() > b.momentum().pT(); });

      //=================== Leptons and Jets multiplicity vetoes ===================//

      // need at exactly two leptons
      if (recon_leptons.size() + tau_vis.size() != 2) { vetoEvent; }
      
      // ensure we have at least 2 jets for each event
      if (recon_jets.size() < 2 ) { vetoEvent; }

      //=================== Summing up weight of surviving events ===================//

      //This will calculate the weight for the N_selected events after some cuts. If an event i passes the cuts, its weight will be calculated here
      //and in the weight variable. If it didn't pass, it will only be in the all_weight variable.

      //deduce the weight of the event that made it here
      const double selected_weight = event.weights()[0];

      //sum up the weights
      //fiducial means after applying cuts (our selected events)
      fid_sumW += selected_weight;

      //Calculate the variance
      var_fid_sumW += std::pow(selected_weight, 2);

      //=================== Calculating Missing Transverse Energy (MET) ===================//

      Particle MET_proj = apply<MissingMomentum>(event, "missing_et"); //.vectorEt().mod();
      Particle MET_smeared_atlas_proj = apply<SmearedMET>(event, "met_smeared_atlas"); //.vectorEt().mod();

      const double eTmiss = pTmiss.pT();

      //=================== Calculating the kinematic variables ===================//

      // write 4-momentum to file
      for(const Particle &p : recon_jets){
        jet_pt.push_back(p.pT());
        jet_rap.push_back(p.rapidity());
        jet_phi.push_back(p.phi());
        jet_mass.push_back(p.mass());
        jet_energy.push_back(p.energy());
        jet_px.push_back(p.px() );
        jet_py.push_back(p.py());
        jet_pz.push_back(p.pz());
        jet_pid.push_back(p.pid());
      }

      for(const Particle &p : recon_leps){
        lep_pt.push_back(p.pT());
        lep_rap.push_back(p.rapidity());
        lep_phi.push_back(p.phi());
        lep_mass.push_back(p.mass());
        lep_energy.push_back(p.energy());
        lep_px.push_back(p.px() );
        lep_py.push_back(p.py());
        lep_pz.push_back(p.pz());
        lep_pid.push_back(p.pid());
      }

      for(const Particle &p : tau_vis){
        tau_pt.push_back(p.momentum.pT);
        tau_rap.push_back(p.momentum.rapidity);
        tau_phi.push_back(p.momentum.phi);
        tau_phi.push_back(p.momentum.phi);
        tau_mass.push_back(p.momentum.mass);
        tau_energy.push_back(p.momentum.energy);
        tau_px.push_back(p.momentum.px);
        tau_py.push_back(p.momentum.py);
        tau_pz.push_back(p.momentum.pz);
        tau_pid.push_back(p.pid());
      }

      for(const Particle &t : tau_vis_comp){
        std::vector<float> tmp_tau_comp_pt;
        std::vector<float> tmp_tau_comp_rap;
        std::vector<float> tmp_tau_comp_phi;
        std::vector<float> tmp_tau_comp_mass;
        std::vector<float> tmp_tau_comp_energy;
        std::vector<float> tmp_tau_comp_px;
        std::vector<float> tmp_tau_comp_py;
        std::vector<float> tmp_tau_comp_pz;
        std::vector<long>  tmp_tau_comp_pid;

        for(const Particle &p : t){
          tmp_tau_comp_pt.push_back(p.momentum);
          tmp_tau_comp_rap.push_back(p.momentum.rapidity);
          tmp_tau_comp_phi.push_back(p.momentum.phi);
          tmp_tau_comp_mass.push_back(p.momentum.mass);
          tmp_tau_comp_energy.push_back(p.momentum.energy);
          tmp_tau_comp_px.push_back(p.momentum.px);
          tmp_tau_comp_py.push_back(p.momentum.py);
          tmp_tau_comp_pz.push_back(p.momentum.pz);
          tmp_tau_comp_pid.push_back(p.pid);
        }
        tau_comp_pt.push_back(tmp_tau_comp_pt);
        tau_comp_rap.push_back(tmp_tau_comp_rap);
        tau_comp_phi.push_back(tmp_tau_comp_phi);
        tau_comp_mass.push_back(tmp_tau_comp_mass);
        tau_comp_energy.push_back(tmp_tau_comp_energy);
        tau_comp_px.push_back(tmp_tau_comp_px);
        tau_comp_py.push_back(tmp_tau_comp_py);
        tau_comp_pz.push_back(tmp_tau_comp_pz);
        tau_comp_pid.push_back(tmp_tau_comp_pid);
      }

      for(const Particle &p : neutrinos){
        nu_pt.push_back(p.pT());
        nu_rap.push_back(p.rapidity());
        nu_phi.push_back(p.phi());
        nu_mass.push_back(p.mass());
        nu_energy.push_back(p.energy());
        nu_px.push_back(p.px() );
        nu_py.push_back(p.py());
        nu_pz.push_back(p.pz());
        nu_pid.push_back(p.pid());
      }

      // for(const Particle &p : tau_vis){
      //   tau_pt.push_back(p.pT());
      //   tau_rap.push_back(p.rapidity());
      //   tau_phi.push_back(p.phi());
      //   tau_mass.push_back(p.mass());
      //   tau_energy.push_back(p.energy());
      //   tau_px.push_back(p.px() );
      //   tau_py.push_back(p.py());
      //   tau_pz.push_back(p.pz());
      // }

      MET_pt =  MET_proj.pT();
      MET_rap = MET_proj.rapidity();
      MET_phi = MET_proj.phi();
      MET =     MET_proj.vectorEt().mod();
      MET_px =  MET_proj.px();
      MET_py =  MET_proj.py();
      MET_pz =  MET_proj.pz();

      ATLAS_smeared_MET_pt  = MET_smeared_atlas_proj.pT();
      ATLAS_smeared_MET_rap = MET_smeared_atlas_proj.rapidity();
      ATLAS_smeared_MET_phi = MET_smeared_atlas_proj.phi();
      ATLAS_smeared_MET     = MET_smeared_atlas_proj.vectorEt().mod();
      ATLAS_smeared_MET_px =  MET_smeared_atlas_proj.px();
      ATLAS_smeared_MET_py =  MET_smeared_atlas_proj.py();
      ATLAS_smeared_MET_pz =  MET_smeared_atlas_proj.pz();

      survive_event_number++;
      OutputTree->Fill();
    }

    //================================================================================//
    //=================================  Finalising  =================================//
    //================================================================================//

    /// Normalise histograms, scale/divide histograms etc., after the run
    void finalize() {
      
      //uncertainty on sum of weights
      const double err_sumW=std::pow(var_sumW,0.5);

      //uncertainty on fiducial sum of weights
      const double err_fid_sumW=std::pow(var_fid_sumW,0.5);

      //total cross-section in femtobarns
      const double total_xs = crossSection();

      //fiducial cross-section in femtobarns
      const double fid_xs = fid_sumW*total_xs/sumW;

      //error on cross-section generated from MadGraph (hard-coded in femtobarns)
      const double err_xs_mg = 0;

      //error on fiducial cross-section
      const double err_fid_xs = pow(pow(err_fid_sumW,2)*pow(total_xs/sumW,2)+pow(err_xs_mg,2)*pow(fid_sumW/sumW,2)+pow(err_sumW,2)*pow(fid_sumW*total_xs/pow(sumW,2),2),0.5);
      std::cout << "SMEAR OPTION: " << SMEAR << ", JET_SMEAR: " << JET_SMEAR << ", MUON_SMEAR: " << MUON_SMEAR << std::endl;
      std::cout << total_xs << std::endl;
      std::cout << "******* IMPORTANT INFORMATION *********" << std::endl;
      std::cout << " Number of Events surviving: " << survive_event_number << std::endl;
      std::cout << " Total XS MadGraph:  " << total_xs << std::endl;
      std::cout << " Total XS error: " << err_xs_mg << std::endl;
      std::cout << " Total Sum of Weights: " << sumW << std::endl;
      std::cout << " Error on Total Sum of Weights: " <<  err_sumW << std::endl;
      std::cout << " Selected Events Sum of Weights: " << fid_sumW << std::endl;
      std::cout << " Error on Selected Events Sum of Weights: " << err_fid_sumW << std::endl;
      std::cout << " Fiducial crossection: " << fid_xs << std::endl;
      std::cout << " Error of fiducial crossection: " << err_fid_xs << std::endl;
      std::cout << "************************" << std::endl;
      
      xs_total                =  total_xs;
      xs_total_err            =  err_xs_mg;          
      sum_of_weights          =  sumW;       
      sum_of_weights_err      =  err_sumW;           
      fid_sum_of_weights      =  fid_sumW;
      fid_sum_of_weights_err  =  err_fid_sumW;
      xs_fid                  =  fid_xs;
      fid_xs_err              =  err_fid_xs;   
      nsurvive_events         =  survive_event_number;        
      ninit_events            =  event_number;
      presel_eff              =  survive_event_number / event_number * 100;   

      MetaDataTree->Fill();
      OutputFile->Close();
      delete OutputFile;
      delete MetaDataTree;
      delete OutputTree;


      

    }

  private:

    // Output file
    // TFile OutputFile(TFile.Open("output.root"));
    // std::unique_ptr<TFile> OutputFile(TFile::Open("output.root", "RECREATE"));
    TFile *OutputFile; // = TFile::Open("output.root", "RECREATE");
    TTree *OutputTree; //("tree", "tree");
    TTree *MetaDataTree; //("metadata", "metadata");

  private:
    std::vector<double> jet_pt;
    std::vector<double> jet_rap;
    std::vector<double> jet_phi;
    std::vector<double> jet_mass;
    std::vector<double> jet_energy;
    std::vector<double> jet_px;
    std::vector<double> jet_py;
    std::vector<double> jet_pz;
    std::vector<long>  jet_pid;

    std::vector<double> lep_pt;
    std::vector<double> lep_rap;
    std::vector<double> lep_phi;
    std::vector<double> lep_mass;
    std::vector<double> lep_energy;
    std::vector<double> lep_px;
    std::vector<double> lep_py;
    std::vector<double> lep_pz;
    std::vector<long>  lep_pid;

    std::vector<std::vector<double>> tau_comp_pt;
    std::vector<std::vector<double>> tau_comp_rap;
    std::vector<std::vector<double>> tau_comp_phi;
    std::vector<std::vector<double>> tau_comp_mass;
    std::vector<std::vector<double>> tau_comp_energy;
    std::vector<std::vector<double>> tau_comp_px;
    std::vector<std::vector<double>> tau_comp_py;
    std::vector<std::vector<double>> tau_comp_pz;
    std::vector<std::vector<long>>  tau_comp_pid;

    std::vector<double> tau_pt;
    std::vector<double> tau_rap;
    std::vector<double> tau_phi;
    std::vector<double> tau_mass;
    std::vector<double> tau_energy;
    std::vector<double> tau_pid;
    std::vector<double> tau_px;
    std::vector<double> tau_py;
    std::vector<double> tau_pz;

    std::vector<double> nu_pt;
    std::vector<double> nu_rap;
    std::vector<double> nu_phi;
    std::vector<double> nu_mass;
    std::vector<double> nu_energy;
    std::vector<double> nu_px;
    std::vector<double> nu_py;
    std::vector<double> nu_pz;
    std::vector<long>  nu_pid;

    double MET_pt;
    double MET_rap;
    double MET_phi;
    double MET;
    double MET_px;
    double MET_py;
    double MET_pz;

    double ATLAS_smeared_MET_pt;
    double ATLAS_smeared_MET_rap;
    double ATLAS_smeared_MET_phi;
    double ATLAS_smeared_MET;
    double ATLAS_smeared_MET_px;
    double ATLAS_smeared_MET_py;
    double ATLAS_smeared_MET_pz;

    std::vector<std::vector<double>> met_p4;

    double xs_total;
    double xs_total_err;
    double sum_of_weights;
    double sum_of_weights_err;
    double fid_sum_of_weights;
    double fid_sum_of_weights_err;
    double xs_fid;
    double fid_xs_err;
    double nsurvive_events;
    double ninit_events;
    double presel_eff;

    int SMEAR;
    bool MET_SMEAR;
    bool JET_SMEAR;
    bool MUON_SMEAR;

    /// @name Histograms
    //@{
    std::map<string, Histo1DPtr> _h;
    std::map<string, Histo2DPtr> _h2d;
    std::map<string, Profile1DPtr> _p;
    std::map<string, CounterPtr> _c;
    //@}

  }; //closes the majorana class at the very top

  //=================================  Plugin hook  =================================//
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(pp_dpp_llvvjj);


}  //closes the Rivet namespace
