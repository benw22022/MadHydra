//=================================  Introduction  =================================//<
// -*- C++ -*-
//Mohamed Aly and Ahmed Abdelmotteleb
//Majorana Neutrino Sensitivity at ATLAS
//Rivet analysis code that runs over our data [signal + background] and applies cuts

//=================================  Libraries  =================================//
#include <iostream>
#include <fstream>
#include <ostream>
#include <map>
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedMET.hh"

#include <math.h>

//==========================================================================//
//==========================================================================//
//=================================  MAIN  =================================//
//==========================================================================//
//==========================================================================//
namespace Rivet {

  //divided into constructor + 3 loops: init, analyze and finalize

  //==================================================================================//
  //=================================  Constructer ==================================//
  //=================================================================================//
  class MLRSM_mumujj : public Analysis {
  public:

    /// Constructor
    MLRSM_mumujj() : Analysis("pp_dpp_llvvjj") {   }

    //=================================  Declaring global variables  =================================//

    double sumW = 0;  //sum of weights
    double var_sumW = 0; //variance of sum of weights
    double fid_sumW = 0; //fiducial sum of weights
    double var_fid_sumW = 0; //variance of fiducial sum of weights
    //=================================  Declaring global functions  =================================//

    //=====LEPTONS
    //Lepton centrality
    double calcLC( const Jets& jets, const FourMomentum& lep ){

      //Define important jets
      // .at() checks if we are within bounds of our array number of elements
      // use .at() instead of [] because it will give an error if we exceed elements however [] won't
      FourMomentum J1 = jets.at(0).momentum(); //highest pT jet
      FourMomentum J2 = jets.at(1).momentum(); //second-highest pT jet

      double LCval = -999; //should be negative??
      if (jets.size() < 2 ) {
        //std::cout << "ERROR: should not allow calculation of Lepton Centrality on jets with less than multiplicity two! " << std::endl;
        return LCval;
      }

      else {  //centrality equation
        double LCval = fabs( ( lep.rapidity() - ( (J1.rapidity() + J2.rapidity() ) /2.0  ) ) / ( J1.rapidity() - J2.rapidity() ) );
        return LCval;
      }

    }

    //boolean that checks if your lepton passes the centrality cut (ATLAS detector)
    bool PassesLeptonCentralityCut( const Jets& jets, const Particle& lep ){
      if ( jets.size() < 2 ) return false;

      if ( calcLC(jets,lep) < 0.4 ) {
        return true;
      }

      else {
        return false;
      }
    }

    //=====JETS
    //Calculate Jet Centrality
    double calcJC( const Jets& jets ){
      //Define important jets
      FourMomentum J1 = jets.at(0).momentum();  //highest pT jet
      FourMomentum J2 = jets.at(1).momentum();  //second-highest pT jet
      double JCval = 999; //should be positive??
      if ( jets.size() < 3 ) {
        //std::cout << "ERROR: should not allow calculation of Jet Centrality on jets with less than multiplicity three! " << std::endl;
        return JCval;
      }

      else {
        FourMomentum J = jets.at(2).momentum();  //take the third jet (jets are always ordered in descending pT order)
        //calculating the centrality
        JCval = fabs(( J.rapidity() - ( (J1.rapidity() + J2.rapidity() ) /2.0  ) ) / ( J1.rapidity() - J2.rapidity() ));
        return JCval;
      }
    }

    //boolean that checks if your lepton (Jet? Ben 14/10/2020) passes the centrality cut (ATLAS detector)
    bool PassesJetCentralityCut( const Jets& jets ){
      if ( jets.size() < 3 ) return true;

      double JCval = calcJC(jets);
      if (JCval < 0.4) { 
        return false;
      }

      else {
        return true;
      }
    }

    //function to find jets between two jet rapidities --> returns true or false
    bool _isBetween(const Jet j, const Jet bj1, const Jet bj2) {
      double y_j = j.rapidity(); //probe jet
      double y_bj1 = bj1.rapidity(); //boundary jet 1
      double y_bj2 = bj2.rapidity(); //boundary jet 2

      double y_min = std::min(y_bj1, y_bj2); //boundary jet with lower rapidity
      double y_max = std::max(y_bj1, y_bj2); //boundary jet with higher rapidity

      if (y_j > y_min && y_j < y_max) return true; //check if our probe lies between the boundary jets
      else return false;
    }


    //function to find the number of gap jets
    int _getNumGapJets(const Jets jets, FourMomentum& thirdJet) {
      if (jets.size() <= 2) return 0;

      // The vector of jets is already sorted by pT. So the boundary jets will be the first two.
      const Jet bj1 = jets.at(0);
      const Jet bj2 = jets.at(1);

      int n_gap_jets = 0;
      // Start loop at the 3rd hardest pT jet
      for (size_t i = 2; i < jets.size(); i++) {
        const Jet j = jets.at(i);
        // If this jet is between the boundary jets and is hard enough, increment counter
        if (_isBetween(j, bj1, bj2)) {
          if (n_gap_jets == 0) thirdJet = j.momentum();
          n_gap_jets++;
        }
      }
      return n_gap_jets;
    }

    //Event counter
    int event_number = 0;
    int survive_event_number = 0;
    int num_events_before_a_cut=0;
    int num_events_after_a_cut=0;
    };

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

      SmearedParticles elecs_eff_smear_atlas(dressed_elecs, ELECTRON_EFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2);
      declare(elecs_eff_smear_atlas, "elecs_eff_smear_atlas");

    //=================================  Tau Decleration  =================================//
    

    //=================================  Jets Decleration  =================================//
    
    // use FastJet to select Jets 
    FastJets fj(fs, FastJets::ANTIKT, 0.4);
    declare(fj,"jets");

    // define jets smeared with ATLAS Run2
    SmearedJets jets_smeared_atlas(fj,JET_SMEAR_ATLAS_RUN2);
    declare(jets_smeared_atlas,"jets_smeared_atlas");

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
        if (deltaR(non_hadronic_leps[j].momentum(), cand_jets[i].momentum()) > 0.3) { recon_jets.push_back(cand_jets[i]); }
      }
    }

    // need at least two jets
    if (recon_jets.size() < 2){ vetoEvent;}


    //=================== work out tau vis ===================//

    Particles non_tau_leps;
    Particles tau_vis;
    std::vector<Particles> tau_vis_comp;
    Particles tau_vis_tmp // unsorted bag of particles from taus, will use to work out which particle goes with which tau if ditau
    
    // get all particles from tau decays, place rest into non-tau containers
    for(const Particle &p : cand_jets){ 
      if (p.fromTau()){ tau_vis_tmp.append(p); }
      else(nontau_jets.append(p);)
    }

    for(const Particle &p : non_hadronic_leps){ 
      if (p.fromTau()){ tau_vis_tmp.append(p); }
      else(nontau_leps.append(p);)
    }

    // get set of mother ids
    std::set<long> mother_ids;
    for(const Particle &p : tau_vis_tmp){ mother_ids.insert(p.mother()); }

    // get all particles orginating from each mother tau
    for(int i{0}; i<mother_ids.size(); i++){
      Particles tmp_tau_comp;
      for(const Particle &p : tau_vis_tmp){
        if(tau_vis_tmp.mother() == mother_ids[i]){tmp_tau_comp.append(p); }
      }
      tau_vis_comp.append(tmp_tau);
    }

    // sum 4-vecs of tau componants to get combined tau
    for(const Particles &tau_particles : tau_vis_tmp){
      Particle comb_tau = tau_particles[0].momentum();
      for(int i{1}; i<tau_particles.size(); i++){
        comb_tau += tau_particles[i].momentum()
      }
      tau_vis.append(tau_comb);
    }

    //=================== Choosing and sorting reconstructed Leptons ===================//
    
    // use non hadronic leptons to select leptons that are far enough from all other leptons
    for (unsigned int i{0}; i < non_tau_leps.size(); i++)
    {
        for (unsigned int j{0}; j < non_tau_leps.size(); j++)
        {
          if(i != j) && deltaR(non_tau_leps[i].momentum(), non_tau_leps[j].momentum()) < 0.2) { recon_leptons.push_back(non_tau_leps[i]); }
        }
      }

    // // Combine reco leps with taus (if any)
    // for(const Particle &p : tau_vis){
    //   recon_leps.append(p)
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

    const double MET = apply<MissingMomentum>(event, "missing_et").vectorEt().mod();
    const double MET_smeared_atlas = apply<SmearedMET>(event, "met_smeared_atlas").vectorEt().mod();

    const double eTmiss = pTmiss.pT();

    //=================== Calculating the kinematic variables ===================//

    // write 4-momentum to file
    for(const Particle &p : recon_jets){
      jet_pt.append(p.pT());
      jet_rap.append(p.rapidity());
      jet_phi.append(p.phi());
      jet_mass.append(p.mass());
      jet_energy.append(p.energy());
      jet_px.append(p.px() );
      jet_py.append(p.py());
      jet_pz.append(p.pz());
      jet_pid.append(p.pid());
    }

    for(const Particle &p : recon_leps){
      lep_pt.append(p.pT());
      lep_rap.append(p.rapidity());
      lep_phi.append(p.phi());
      lep_mass.append(p.mass());
      lep_energy.append(p.energy());
      lep_px.append(p.px() );
      lep_py.append(p.py());
      lep_pz.append(p.pz());
      lep_pid.append(p.pid());
    }

    for(const Particle &p : tau_vis){
      tau_pt.append(p.pT());
      tau_rap.append(p.rapidity());
      tau_phi.append(p.phi());
      tau_mass.append(p.mass());
      tau_energy.append(p.energy());
      tau_px.append(p.px() );
      tau_py.append(p.py());
      tau_pz.append(p.pz());
      tau_pid.append(p.pid());
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
        tmp_tau_comp_pt.append(p.pT());
        tmp_tau_comp_rap.append(p.rapidity());
        tmp_tau_comp_phi.append(p.phi());
        tmp_tau_comp_mass.append(p.mass());
        tmp_tau_comp_energy.append(p.energy());
        tmp_tau_comp_px.append(p.px() );
        tmp_tau_comp_py.append(p.py());
        tmp_tau_comp_pz.append(p.pz());
        tmp_tau_comp_pid.append(p.pid());
      }
      tau_comp_pt.append(tmp_tau_comp_pt)
      tau_comp_rap.append(tmp_tau_comp_rap)
      tau_comp_phi.append(tmp_tau_comp_phi)
      tau_comp_mass.append(tmp_tau_comp_mass)
      tau_comp_energy.append(tmp_tau_comp_energy)
      tau_comp_px.append(tmp_tau_comp_px)
      tau_comp_py.append(tmp_tau_comp_py)
      tau_comp_pz.append(tmp_tau_comp_pz)
      tau_comp_pid.append(tmp_tau_comp_pid)
    }

    for(const Particle &p : neutrinos){
      nu_pt.append(p.pT());
      nu_rap.append(p.rapidity());
      nu_phi.append(p.phi());
      nu_mass.append(p.mass());
      nu_energy.append(p.energy());
      nu_px.append(p.px() );
      nu_py.append(p.py());
      nu_pz.append(p.pz());
      nu_pid.append(p.pid());
    }

    for(const Particle &p : tau_vis){
      tau_pt.append(p.pT());
      tau_rap.append(p.rapidity());
      tau_phi.append(p.phi());
      tau_mass.append(p.mass());
      tau_energy.append(p.energy());
      tau_px.append(p.px() );
      tau_py.append(p.py());
      tau_pz.append(p.pz());
    }

    const double MET = apply<MissingMomentum>(event, "missing_et");
    const double ATLAS_smeared_MET = apply<SmearedMET>(event, "met_smeared_atlas");

    MET_pt = MET.pT();
    MET_rap = MET.rapidity();
    MET_phi = MET.phi();
    MET = MET.vectorEt().mod();
    MET_px = MET.px();
    MET_py = MET.py();
    MET_pz = MET.pz();

    ATLAS_smeared_MET_pt =  ATLAS_smeared_MET.pT();
    ATLAS_smeared_MET_rap = ATLAS_smeared_MET.rapidity();
    ATLAS_smeared_MET_phi = ATLAS_smeared_MET.phi();
    ATLAS_smeared_MET =     ATLAS_smeared_MET.vectorEt().mod();
    ATLAS_smeared_MET_px =  ATLAS_smeared_MET.px();
    ATLAS_smeared_MET_py =  ATLAS_smeared_MET.py();
    ATLAS_smeared_MET_pz =  ATLAS_smeared_MET.pz();

    survive_event_number++;
    OutputTree.Fill();
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
  fid_xs                  =  fid_xs;
  fid_xs_err              =  err_fid_xs;   
  nsurvive_events         =  survive_event_number;        
  ninit_events            =  event_number;
  presel_eff              =  survive_event_number / event_number * 100;   

  MetaDataTree.Fill()

}

private:

  // Output file
  TFile OutputFile(TFile::open("output.root"))

  TTree OutputTree("tree", "tree")
  TTree MetaDataTree("metadata", "metadata")

  std::vector<float> jet_pt;
  std::vector<float> jet_rap;
  std::vector<float> jet_phi;
  std::vector<float> jet_mass;
  std::vector<float> jet_energy;
  std::vector<float> jet_px;
  std::vector<float> jet_py;
  std::vector<float> jet_pz;
  std::vector<long>  jet_pid;

  std::vector<float> lep_pt;
  std::vector<float> lep_rap;
  std::vector<float> lep_phi;
  std::vector<float> lep_mass;
  std::vector<float> lep_energy;
  std::vector<float> lep_px;
  std::vector<float> lep_py;
  std::vector<float> lep_pz;
  std::vector<long>  lep_pid;

  std::vector<std::vector<float>> tau_comp_pt;
  std::vector<std::vector<float>> tau_comp_rap;
  std::vector<std::vector<float>> tau_comp_phi;
  std::vector<std::vector<float>> tau_comp_mass;
  std::vector<std::vector<float>> tau_comp_energy;
  std::vector<std::vector<float>> tau_comp_px;
  std::vector<std::vector<float>> tau_comp_py;
  std::vector<std::vector<float>> tau_comp_pz;
  std::vector<std::vector<long>>  tau_comp_pid;

  std::vector<std::vector<float>> tau_pt;
  std::vector<std::vector<float>> tau_rap;
  std::vector<std::vector<float>> tau_phi;
  std::vector<std::vector<float>> tau_mass;
  std::vector<std::vector<float>> tau_energy;
  std::vector<std::vector<float>> tau_px;
  std::vector<std::vector<float>> tau_py;
  std::vector<std::vector<float>> tau_pz;

  std::vector<float> nu_pt;
  std::vector<float> nu_rap;
  std::vector<float> nu_phi;
  std::vector<float> nu_mass;
  std::vector<float> nu_energy;
  std::vector<float> nu_px;
  std::vector<float> nu_py;
  std::vector<float> nu_pz;
  std::vector<long>  nu_pid;

  std::vector<float> MET_pt;
  std::vector<float> MET_rap;
  std::vector<float> MET_phi;
  std::vector<float> MET;
  std::vector<float> MET_px;
  std::vector<float> MET_py;
  std::vector<float> MET_pz;

  std::vector<float> ATLAS_smeared_MET_pt;
  std::vector<float> ATLAS_smeared_MET_rap;
  std::vector<float> ATLAS_smeared_MET_phi;
  std::vector<float> ATLAS_smeared_MET;
  std::vector<float> ATLAS_smeared_MET_px;
  std::vector<float> ATLAS_smeared_MET_py;
  std::vector<float> ATLAS_smeared_MET_pz;

  std::vector<std::vector<float>> met_p4;

  OutputTree.Branch("jet_pt",     &jet_pt,     "jet_pt/F");
  OutputTree.Branch("jet_rap",    &jet_rap,    "jet_rap/F");
  OutputTree.Branch("jet_phi",    &jet_pt,     "jet_phi/F");
  OutputTree.Branch("jet_mass",   &jet_mass,   "jet_mass/F");
  OutputTree.Branch("jet_energy", &jet_energy, "jet_energy/F");
  OutputTree.Branch("jet_px",     &jet_px,     "jet_px/F");
  OutputTree.Branch("jet_py",     &jet_py,     "jet_py/F");
  OutputTree.Branch("jet_pz",     &jet_pz,     "jet_pz/F");
  OutputTree.Branch("jet_pid",    &jet_pid,    "jet_pid/L");

  OutputTree.Branch("lep_pt",     &lep_pt,     "lep_pt/F");
  OutputTree.Branch("lep_rap",    &lep_rap,    "lep_rap/F");
  OutputTree.Branch("lep_phi",    &lep_pt,     "lep_phi/F");
  OutputTree.Branch("lep_mass",   &lep_mass,   "lep_mass/F");
  OutputTree.Branch("lep_energy", &lep_energy, "lep_energy/F");
  OutputTree.Branch("lep_px",     &lep_px,     "lep_px/F");
  OutputTree.Branch("lep_py",     &lep_py,     "lep_py/F");
  OutputTree.Branch("lep_pz",     &lep_pz,     "lep_pz/F");
  OutputTree.Branch("lep_pid",    &lep_pid,    "lep_pid/L");

  OutputTree.Branch("tau_pt",     &tau_pt,     "tau_pt/F");
  OutputTree.Branch("tau_rap",    &tau_rap,    "tau_rap/F");
  OutputTree.Branch("tau_phi",    &tau_pt,     "tau_phi/F");
  OutputTree.Branch("tau_mass",   &tau_mass,   "tau_mass/F");
  OutputTree.Branch("tau_energy", &tau_energy, "tau_energy/F");
  OutputTree.Branch("tau_px",     &tau_px,     "tau_px/F");
  OutputTree.Branch("tau_py",     &tau_py,     "tau_py/F");
  OutputTree.Branch("tau_pz",     &tau_pz,     "tau_pz/F");
  OutputTree.Branch("tau_pid",    &tau_pid,    "tau_pid/L");

  OutputTree.Branch("tau_comp_pt",     &tau_comp_pt,     "tau_comp_pt/F");
  OutputTree.Branch("tau_comp_rap",    &tau_comp_rap,    "tau_comp_rap/F");
  OutputTree.Branch("tau_comp_phi",    &tau_comp_pt,     "tau_comp_phi/F");
  OutputTree.Branch("tau_comp_mass",   &tau_comp_mass,   "tau_comp_mass/F");
  OutputTree.Branch("tau_comp_energy", &tau_comp_energy, "tau_comp_energy/F");
  OutputTree.Branch("tau_comp_px",     &tau_comp_px,     "tau_comp_px/F");
  OutputTree.Branch("tau_comp_py",     &tau_comp_py,     "tau_comp_py/F");
  OutputTree.Branch("tau_comp_pz",     &tau_comp_pz,     "tau_comp_pz/F");
  OutputTree.Branch("tau_comp_pid",    &tau_comp_pid,    "tau_comp_pid/L");

  OutputTree.Branch("MET_pt",     &MET_pt,  "MET_pt/F");
  OutputTree.Branch("MET_rap",    &MET_rap, "MET_rap/F");
  OutputTree.Branch("MET_phi",    &MET_pt,  "MET_phi/F");
  OutputTree.Branch("MET",        &MET,     "MET/F");
  OutputTree.Branch("MET_px",     &MET_px,  "MET_px/F");
  OutputTree.Branch("MET_py",     &MET_py,  "MET_py/F");
  OutputTree.Branch("MET_pz",     &MET_pz,  "MET_pz/F");

  OutputTree.Branch("ATLAS_smeared_MET_pt",     &ATLAS_smeared_MET_pt,  "ATLAS_smeared_MET_pt/F");
  OutputTree.Branch("ATLAS_smeared_MET_rap",    &ATLAS_smeared_MET_rap, "ATLAS_smeared_MET_rap/F");
  OutputTree.Branch("ATLAS_smeared_MET_phi",    &ATLAS_smeared_MET_pt,  "ATLAS_smeared_MET_phi/F");
  OutputTree.Branch("ATLAS_smeared_MET",        &ATLAS_smeared_MET,     "ATLAS_smeared_MET/F");
  OutputTree.Branch("ATLAS_smeared_MET_px",     &ATLAS_smeared_MET_px,  "ATLAS_smeared_MET_px/F");
  OutputTree.Branch("ATLAS_smeared_MET_py",     &ATLAS_smeared_MET_py,  "ATLAS_smeared_MET_py/F");
  OutputTree.Branch("ATLAS_smeared_MET_pz",     &ATLAS_smeared_MET_pz,  "ATLAS_smeared_MET_pz/F");

  MetaDataTree.Branch("total_xs",               &xs_total,              "total_xs_madgraph/F");
  MetaDataTree.Branch("total_xs_error",         &xs_total_err,          "total_xs_error/F");
  MetaDataTree.Branch("sum_of_weights",         &sum_of_weights,         "sum_of_weights/F");
  MetaDataTree.Branch("sum_of_weights_err",     &sum_of_weights_err,     "sum_of_weights_err/F");
  MetaDataTree.Branch("fid_sum_of_weights",     &fid_sum_of_weights,     "fid_sum_of_weights/F");
  MetaDataTree.Branch("fid_sum_of_weights_err", &fid_sum_of_weights_err, "fid_sum_of_weights_err/F");
  MetaDataTree.Branch("fid_xs",                 &fid_xs,                 "fid_xs/F");
  MetaDataTree.Branch("fid_xs_err",             &fid_xs_err,             "fid_xs_err/F");
  MetaDataTree.Branch("nsurvive_events",        &nsurvive_events,        "nsurvive_events/F");
  MetaDataTree.Branch("ninit_events",           &ninit_events,           "ninit_events/F");
  MetaDataTree.Branch("presel_eff",             &presel_eff,             "preself_eff/F");

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
