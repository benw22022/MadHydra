//=================================  Introduction  =================================//
// -*- C++ -*-
// Writes an NTuple of lepton, jet, tau, photon and neutrino kinematics
// Requires at least two jets and two leptons with very loose selection requirements
//
//=================================  Libraries  =================================//
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "TFile.h"
#include "TTree.h"
#include <math.h>


  // Helpful CutFlow struct
  struct CutFlow{

    std::map<std::string, int> m_cutflow;

    void BookCut(const std::string& cut){
      std::pair<std::string, int> entry{cut, 0};
      m_cutflow.insert(entry);
    }

    void Count(const std::string& cut, const int& count=1){
      if (m_cutflow.find(cut) == m_cutflow.end()){
        BookCut(cut);
      }
      m_cutflow[cut] += count;
    }

    void Print() const {
      std::cout << "*************************" << std::endl;
      std::cout << "Cut : count" << std::endl;
      for (auto const& [key, val] : m_cutflow){
        std::cout << key << " : " << val << std::endl;
      }
      std::cout << "*************************" << std::endl;
    }

  };

//==========================================================================//
//=================================  MAIN  =================================//
//==========================================================================//

//divided into constructor + 3 loops: init, analyze and finalize

namespace Rivet {
  

  //==================================================================================//
  //=================================  Constructer ==================================//
  //=================================================================================//
  class minimal_2l2j : public Analysis {

    public:

      double calc_deltaR(const double& y1, const double& y2, const double& phi1, const double& phi2){

      return sqrt(pow((y1 - y2), 2) + pow((phi1 - phi2), 2));

      };

      /// Constructor
      minimal_2l2j() : Analysis("minimal_2l2j") {   }

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

      // Cutflow
      CutFlow cutflow;

      // Output file
      TFile *m_OutputFile;
      TTree *m_OutputTree;
      TTree *m_MetaDataTree;
      bool m_has_finalized{false};

      std::vector<double> b_jet_pt;
      std::vector<double> b_jet_rap;
      std::vector<double> b_jet_phi;
      std::vector<double> b_jet_mass;
      std::vector<double> b_jet_energy;
      std::vector<double> b_jet_px;
      std::vector<double> b_jet_py;
      std::vector<double> b_jet_pz;

      std::vector<double> b_lep_pt;
      std::vector<double> b_lep_rap;
      std::vector<double> b_lep_phi;
      std::vector<double> b_lep_mass;
      std::vector<double> b_lep_energy;
      std::vector<double> b_lep_px;
      std::vector<double> b_lep_py;
      std::vector<double> b_lep_pz;
      std::vector<long>   b_lep_pid;

      std::vector<double> b_tau_pt;
      std::vector<double> b_tau_rap;
      std::vector<double> b_tau_phi;
      std::vector<double> b_tau_mass;
      std::vector<double> b_tau_energy;
      std::vector<double> b_tau_pid;
      std::vector<double> b_tau_px;
      std::vector<double> b_tau_py;
      std::vector<double> b_tau_pz;

      std::vector<double> b_tau_daughter_pt;
      std::vector<double> b_tau_daughter_rap;
      std::vector<double> b_tau_daughter_phi;
      std::vector<double> b_tau_daughter_mass;
      std::vector<double> b_tau_daughter_energy;
      std::vector<double> b_tau_daughter_px;
      std::vector<double> b_tau_daughter_py;
      std::vector<double> b_tau_daughter_pz;
      std::vector<double> b_tau_daughter_pid;

      std::vector<double> b_neutrino_pt;
      std::vector<double> b_neutrino_rap;
      std::vector<double> b_neutrino_phi;
      std::vector<double> b_neutrino_mass;
      std::vector<double> b_neutrino_energy;
      std::vector<double> b_neutrino_px;
      std::vector<double> b_neutrino_py;
      std::vector<double> b_neutrino_pz;
      std::vector<long>   b_neutrino_pid;

      std::vector<double> b_photon_pt;
      std::vector<double> b_photon_rap;
      std::vector<double> b_photon_phi;
      std::vector<double> b_photon_mass;
      std::vector<double> b_photon_energy;
      std::vector<double> b_photon_px;
      std::vector<double> b_photon_py;
      std::vector<double> b_photon_pz;
      std::vector<long>   b_photon_pid;

      double b_HTlep;
      double b_HTjets;
      double b_eTmiss;
      double b_meff;

      double b_weight;

      double b_total_xs;
      double b_err_xs_mg;
      double b_sumW;
      double b_err_sumW;
      double b_fid_sumW;
      double b_err_fid_sumW;
      double b_fid_xs;
      double b_err_fid_xs;
      double b_survive_event_number;
      double b_event_number;
      double b_presel_eff;
      bool _use_fiducial_lepton_efficiency;

      //=================================  Declaring global functions  =================================//

      //====================================================================================//
      //=================================  Initialisation  =================================//
      //====================================================================================//

      void init() {
        
        

        // Output file
        m_OutputFile = new TFile("output.root", "RECREATE");
        m_OutputTree = new TTree("tree", "tree");
        m_MetaDataTree = new TTree("metadata", "metadata");
        m_OutputFile->cd();

        m_OutputTree->Branch("jet_pt",     &b_jet_pt    );
        m_OutputTree->Branch("jet_rap",    &b_jet_rap   );
        m_OutputTree->Branch("jet_phi",    &b_jet_phi    );
        m_OutputTree->Branch("jet_mass",   &b_jet_mass  );
        m_OutputTree->Branch("jet_energy", &b_jet_energy );
        m_OutputTree->Branch("jet_px",     &b_jet_px    );
        m_OutputTree->Branch("jet_py",     &b_jet_py    );
        m_OutputTree->Branch("jet_pz",     &b_jet_pz    );

        m_OutputTree->Branch("lep_pt",     &b_lep_pt    );
        m_OutputTree->Branch("lep_rap",    &b_lep_rap   );
        m_OutputTree->Branch("lep_phi",    &b_lep_phi    );
        m_OutputTree->Branch("lep_mass",   &b_lep_mass  );
        m_OutputTree->Branch("lep_energy", &b_lep_energy);
        m_OutputTree->Branch("lep_px",     &b_lep_px    );
        m_OutputTree->Branch("lep_py",     &b_lep_py    );
        m_OutputTree->Branch("lep_pz",     &b_lep_pz    );
        m_OutputTree->Branch("lep_pid",    &b_lep_pid   );

        m_OutputTree->Branch("tau_pt",     &b_tau_pt   );
        m_OutputTree->Branch("tau_rap",    &b_tau_rap   );
        m_OutputTree->Branch("tau_phi",    &b_tau_phi    );
        m_OutputTree->Branch("tau_mass",   &b_tau_mass  );
        m_OutputTree->Branch("tau_energy", &b_tau_energy);
        m_OutputTree->Branch("tau_px",     &b_tau_px    );
        m_OutputTree->Branch("tau_py",     &b_tau_py    );
        m_OutputTree->Branch("tau_pz",     &b_tau_pz    );
        m_OutputTree->Branch("tau_pid",    &b_tau_pid   );

        m_OutputTree->Branch("tau_daughter_pt",     &b_tau_daughter_pt    );
        m_OutputTree->Branch("tau_daughter_rap",    &b_tau_daughter_rap   );
        m_OutputTree->Branch("tau_daughter_phi",    &b_tau_daughter_phi    );
        m_OutputTree->Branch("tau_daughter_mass",   &b_tau_daughter_mass  );
        m_OutputTree->Branch("tau_daughter_energy", &b_tau_daughter_energy);
        m_OutputTree->Branch("tau_daughter_px",     &b_tau_daughter_px    );
        m_OutputTree->Branch("tau_daughter_py",     &b_tau_daughter_py    );
        m_OutputTree->Branch("tau_daughter_pz",     &b_tau_daughter_pz    );
        m_OutputTree->Branch("tau_daughter_pid",    &b_tau_daughter_pid   );

        m_OutputTree->Branch("photon_pt",     &b_photon_pt   );
        m_OutputTree->Branch("photon_rap",    &b_photon_rap   );
        m_OutputTree->Branch("photon_phi",    &b_photon_phi    );
        m_OutputTree->Branch("photon_mass",   &b_photon_mass  );
        m_OutputTree->Branch("photon_energy", &b_photon_energy);
        m_OutputTree->Branch("photon_px",     &b_photon_px    );
        m_OutputTree->Branch("photon_py",     &b_photon_py    );
        m_OutputTree->Branch("photon_pz",     &b_photon_pz    );
        m_OutputTree->Branch("photon_pid",    &b_photon_pid   );

        m_OutputTree->Branch("neutrino_pt",     &b_neutrino_pt   );
        m_OutputTree->Branch("neutrino_rap",    &b_neutrino_rap   );
        m_OutputTree->Branch("neutrino_phi",    &b_neutrino_phi    );
        m_OutputTree->Branch("neutrino_mass",   &b_neutrino_mass  );
        m_OutputTree->Branch("neutrino_energy", &b_neutrino_energy);
        m_OutputTree->Branch("neutrino_px",     &b_neutrino_px    );
        m_OutputTree->Branch("neutrino_py",     &b_neutrino_py    );
        m_OutputTree->Branch("neutrino_pz",     &b_neutrino_pz    );
        m_OutputTree->Branch("neutrino_pid",    &b_neutrino_pid   );

        m_OutputTree->Branch("eTmiss", &b_eTmiss,  "eTmiss/D");
        m_OutputTree->Branch("HTlep",  &b_HTlep,   "HTlep/D");
        m_OutputTree->Branch("HTjets", &b_HTjets,  "HTjets/D");
        m_OutputTree->Branch("meff",   &b_meff,    "meff/D");

        m_OutputTree->Branch("weight",   &b_weight,    "weight/D");
        
        m_MetaDataTree->Branch("total_xs",               &b_total_xs,             "total_xs_madgraph/D");
        m_MetaDataTree->Branch("total_xs_error",         &b_err_xs_mg,            "total_xs_error/D");
        m_MetaDataTree->Branch("sum_of_weights",         &b_sumW,                 "sum_of_weights/D");
        m_MetaDataTree->Branch("sum_of_weights_err",     &b_err_sumW,             "sum_of_weights_err/D");
        m_MetaDataTree->Branch("fid_sum_of_weights",     &b_fid_sumW,             "fid_sum_of_weights/D");
        m_MetaDataTree->Branch("fid_sum_of_weights_err", &b_err_fid_sumW,         "fid_sum_of_weights_err/D");
        m_MetaDataTree->Branch("fid_xs",                 &b_fid_xs,               "fid_xs/D");
        m_MetaDataTree->Branch("fid_xs_err",             &b_err_fid_xs,           "fid_xs_err/D");
        m_MetaDataTree->Branch("nsurvive_events",        &b_survive_event_number, "nsurvive_events/D");
        m_MetaDataTree->Branch("ninit_events",           &b_event_number,          "ninit_events/D");
        m_MetaDataTree->Branch("presel_eff",             &b_presel_eff,            "preself_eff/D");
        m_MetaDataTree->Branch("cutflow",                &cutflow,                 "cutflow");

        // Final state including all charged and neutral particles
        const FinalState fs((Cuts::etaIn(-5.0, 5.0) && Cuts::pT >=  1*GeV));
        declare(fs, "FS");

        // Final state including all charged particles
        declare(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 1*GeV), "CFS");

        // Final state including all visible particles (to calculate MET, Jets etc.)
        declare(VisibleFinalState(Cuts::abseta < 5.0), "VFS");

        // Final state including all AntiKt 04 Jets
        VetoedFinalState vfs;
        vfs.addVetoPairId(PID::MUON); // this removes Muons/Taus from visible final state - prevents inclusion in jets?
        declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

        // Final state including all unstable particles (including taus)
        declare(UnstableParticles(Cuts::abseta < 5.0 && Cuts::pT > 5*GeV), "UFS");

        // Final state including all electrons
        IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 5*GeV);
        elecs.acceptIdPair(PID::ELECTRON);
        declare(elecs, "elecs");

        // Final state including all muons
        IdentifiedFinalState muons(Cuts::abseta < 2.5 && Cuts::pT > 5*GeV);
        muons.acceptIdPair(PID::MUON);
        declare(muons, "muons");

        //=================================  Photon and Neutrino Decleration  =================================//
        // define photon
        IdentifiedFinalState photon(fs);
        photon.acceptIdPair(PID::PHOTON); 
        declare(photon,"photons");

        // identify the Neutrinos
        IdentifiedFinalState neutrinos(fs);
        neutrinos.acceptNeutrinos();
        declare(neutrinos, "neutrinos");
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
        cutflow.Count("Initial");

        if(event_number % 1000 == 0){std::cout << "Proccessing event: " << event_number << std::endl;}

        //=================== weight and variance ===================//

        // This will calculate the weight for each of the N events that rivet processes

        // weight of each event (Rivet has for some reason removed the weight() method, the 0th element should be the default weight from looking at the output of checkMetaSG.py)
        const double weight = event.weights()[0];

        // summing over all the weights (sum of weights)
        sumW += weight;

        // Sum up the weights^2 --> This happens to be equal to the variance
        var_sumW += std::pow(weight, 2);

        //================================= Building Candidate Particles =================================//
        
        // Jets (all anti-kt R=0.4 jets with pT > 20 GeV and eta < 4.9)
        Jets jet_candidates;
        Jets tau_jet_candidates;
        for (const Jet& jet : apply<FastJets>(event, "AntiKtJets04").jetsByPt(20*GeV)) {
          if (jet.abseta() < 4.9) {
            if(!jet.tauTagged()) jet_candidates.push_back(jet); 
            else tau_jet_candidates.push_back(jet);
          }
          
        }

        // Veto event if we have < 2 jets
        cutflow.BookCut("< 2 jets");
        if (jet_candidates.size() < 2){
          cutflow.Count("< 2 jets");
          vetoEvent;
        }

        Particles recon_leptons;
        Particles tau_daughter_leptons;
        const Particles charged_tracks    = apply<ChargedFinalState>(event, "CFS").particles();
        const Particles visible_particles = apply<VisibleFinalState>(event, "VFS").particles();
        
        // Muons
        Particles muon_candidates;
        Particles tau_muon_candidates;
        for (const Particle& mu : apply<IdentifiedFinalState>(event, "muons").particlesByPt()) {
          if (!mu.fromTau()){ muon_candidates.push_back(mu); recon_leptons.push_back(mu);}
          else{  tau_daughter_leptons.push_back(mu); }
        }


        // Electrons
        Particles electron_candidates;
        Particles tau_electron_candidates;
        for (const Particle& e : apply<IdentifiedFinalState>(event, "elecs").particlesByPt()) {
            if (!e.fromTau()){ electron_candidates.push_back(e); recon_leptons.push_back(e);}
            else{ tau_daughter_leptons.push_back(e); }
        }


        // Taus
        /// @todo This could benefit from a tau finder projection
        Particles tau_candidates;
        for (const Particle& tau : apply<UnstableParticles>(event, "UFS").particlesByPt()) {

          // Only pick taus out of all unstable particles
          if (tau.abspid() != PID::TAU) continue;

          // Check that tau has decayed into daughter particles
          /// @todo Huh? Unstable taus with no decay vtx? Can use Particle.isStable()? But why in this situation?
          if (tau.genParticle()->end_vertex() == 0) continue;

          // Calculate visible tau pT from pT of tau neutrino in tau decay for pT and |eta| cuts
          FourMomentum daughter_tau_neutrino_momentum = get_tau_neutrino_mom(tau);
          Particle tau_vis = tau;
          tau_vis.setMomentum(tau.momentum()-daughter_tau_neutrino_momentum);
          tau_candidates.push_back(tau_vis);
        }


        // Photons
        Particles photon_candidates;
        for (const Particle& a : apply<IdentifiedFinalState>(event, "photons").particlesByPt()) {
          photon_candidates.push_back(a);
        }

        // neutrinos
        Particles neutrino_candidates;
        for (const Particle& nu : apply<IdentifiedFinalState>(event, "neutrinos").particlesByPt()) {
          neutrino_candidates.push_back(nu);
        }


        // ETmiss
        Particles vfs_particles = apply<VisibleFinalState>(event, "VFS").particles();
        FourMomentum pTmiss;
        for (const Particle& p : vfs_particles) pTmiss -= p.momentum();
        double eTmiss = pTmiss.pT()/GeV;

        // Veto event if we have < 2 leptons
        if (recon_leptons.size() + tau_candidates.size() < 2){
          cutflow.Count("< 2 leptons");
          vetoEvent;
        }

        // Calculate HTlep, fill lepton pT histograms & store chosen combination of 3 leptons
        double HTlep = 0.;
        Particles chosen_leptons;
        
        for(const Particle &p : recon_leptons){
            
          b_lep_pt.push_back(p.pT());
          b_lep_rap.push_back(p.rapidity());
          b_lep_phi.push_back(p.phi());
          b_lep_mass.push_back(p.mass());
          b_lep_energy.push_back(p.energy());
          b_lep_px.push_back(p.px() );
          b_lep_py.push_back(p.py());
          b_lep_pz.push_back(p.pz());
          b_lep_pid.push_back(p.pid());
        }

        for(const Particle &p : tau_candidates){
          b_tau_pt.push_back(p.pT());
          b_tau_rap.push_back(p.rapidity());
          b_tau_phi.push_back(p.phi());
          b_tau_mass.push_back(p.mass());
          b_tau_energy.push_back(p.energy());
          b_tau_px.push_back(p.px() );
          b_tau_py.push_back(p.py());
          b_tau_pz.push_back(p.pz());
          b_tau_pid.push_back(p.pid());
        }

        for(const Jet &p : jet_candidates){
          b_jet_pt.push_back(p.pT());
          b_jet_rap.push_back(p.rapidity());
          b_jet_phi.push_back(p.phi());
          b_jet_mass.push_back(p.mass());
          b_jet_energy.push_back(p.energy());
          b_jet_px.push_back(p.px() );
          b_jet_py.push_back(p.py());
          b_jet_pz.push_back(p.pz());
        }

        for(const Particle &p : photon_candidates){
            b_photon_pt.push_back(p.pT());
            b_photon_rap.push_back(p.rapidity());
            b_photon_phi.push_back(p.phi());
            b_photon_mass.push_back(p.mass());
            b_photon_energy.push_back(p.energy());
            b_photon_px.push_back(p.px() );
            b_photon_py.push_back(p.py());
            b_photon_pz.push_back(p.pz());
            b_photon_pid.push_back(p.pid());
        }
        
        for(const Particle &p : neutrino_candidates){
          b_neutrino_pt.push_back(p.pT());
          b_neutrino_rap.push_back(p.rapidity());
          b_neutrino_phi.push_back(p.phi());
          b_neutrino_mass.push_back(p.mass());
          b_neutrino_energy.push_back(p.energy());
          b_neutrino_px.push_back(p.px() );
          b_neutrino_py.push_back(p.py());
          b_neutrino_pz.push_back(p.pz());
          b_neutrino_pid.push_back(p.pid());
        }

        for(const Particle &p : tau_daughter_leptons){
          b_tau_daughter_pt.push_back(p.pT());
          b_tau_daughter_rap.push_back(p.rapidity());
          b_tau_daughter_phi.push_back(p.phi());
          b_tau_daughter_mass.push_back(p.mass());
          b_tau_daughter_energy.push_back(p.energy());
          b_tau_daughter_px.push_back(p.px() );
          b_tau_daughter_py.push_back(p.py());
          b_tau_daughter_pz.push_back(p.pz());
          b_tau_daughter_pid.push_back(p.pid());
        }

        for(const Jet &p : tau_jet_candidates){
          b_tau_daughter_pt.push_back(p.pT());
          b_tau_daughter_rap.push_back(p.rapidity());
          b_tau_daughter_phi.push_back(p.phi());
          b_tau_daughter_mass.push_back(p.mass());
          b_tau_daughter_energy.push_back(p.energy());
          b_tau_daughter_px.push_back(p.px() );
          b_tau_daughter_py.push_back(p.py());
          b_tau_daughter_pz.push_back(p.pz());
          b_tau_daughter_pid.push_back(-999);
        }

      // Calculate HTjets
      double HTjets = 0.;
      for ( const Jet & jet : jet_candidates)
        HTjets += jet.perp()/GeV;


      // Calculate meff
      double meff = eTmiss + HTjets;
      for ( const Particle & e  : electron_candidates )  {
        meff += e.perp()/GeV;
      }
      for ( const Particle & mu : muon_candidates)  {
        meff += mu.perp()/GeV;
      }
      for ( const Particle & tau : tau_candidates)  {
        meff += tau.perp()/GeV;
      }

      b_HTlep = HTlep;
      b_HTjets = HTjets;
      b_eTmiss = eTmiss;
      b_meff = meff;

      b_weight= event.weights()[0];

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

      survive_event_number++;
      m_OutputTree->Fill();

      // Clear vectors
      b_jet_pt.clear();
      b_jet_rap.clear();
      b_jet_phi.clear();
      b_jet_mass.clear();
      b_jet_energy.clear();
      b_jet_px.clear();
      b_jet_py.clear();
      b_jet_pz.clear();
      
      b_lep_pt.clear();
      b_lep_rap.clear();
      b_lep_phi.clear();
      b_lep_mass.clear();
      b_lep_energy.clear();
      b_lep_px.clear();
      b_lep_py.clear();
      b_lep_pz.clear();
      b_lep_pid.clear();
      
      b_tau_pt.clear();
      b_tau_rap.clear();
      b_tau_phi.clear();
      b_tau_mass.clear();
      b_tau_energy.clear();
      b_tau_pid.clear();
      b_tau_px.clear();
      b_tau_py.clear();
      b_tau_pz.clear();

      b_tau_daughter_pt.clear();
      b_tau_daughter_rap.clear();
      b_tau_daughter_phi.clear();
      b_tau_daughter_mass.clear();
      b_tau_daughter_energy.clear();
      b_tau_daughter_pid.clear();
      b_tau_daughter_px.clear();
      b_tau_daughter_py.clear();
      b_tau_daughter_pz.clear();
      
      b_photon_pt.clear();
      b_photon_rap.clear();
      b_photon_phi.clear();
      b_photon_mass.clear();
      b_photon_energy.clear();
      b_photon_px.clear();
      b_photon_py.clear();
      b_photon_pz.clear();
      b_photon_pid.clear();
      
      b_neutrino_pt.clear();
      b_neutrino_rap.clear();
      b_neutrino_phi.clear();
      b_neutrino_mass.clear();
      b_neutrino_energy.clear();
      b_neutrino_px.clear();
      b_neutrino_py.clear();
      b_neutrino_pz.clear();
      b_neutrino_pid.clear();

      // Count #events that survived selection
      cutflow.Count("final");

    }

    //================================================================================//
    //=================================  Finalising  =================================//
    //================================================================================//

    /// Normalise histograms, scale/divide histograms etc., after the run
    void finalize() {

      // Finalize gets called twice for some reason - without this flag output file gets deleted twice -> seg fault      
      if (!m_has_finalized){

        m_has_finalized = true;
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
        
        b_total_xs                =  total_xs;
        b_err_xs_mg               =  err_xs_mg;
        b_sumW                    =  sumW;
        b_err_sumW                =  err_sumW;
        b_fid_sumW                =  fid_sumW;
        b_err_fid_sumW            =  err_fid_sumW;
        b_fid_xs                  =  fid_xs;
        b_err_fid_xs              =  err_fid_xs;
        b_survive_event_number    =  survive_event_number;
        b_event_number            =  event_number;
        b_presel_eff              =  ((double)survive_event_number / (double)event_number) * 100;

        cutflow.Print();
        
        m_MetaDataTree->Fill();
        m_MetaDataTree->Write();
        m_OutputTree->Write();
        
        delete m_OutputFile;
      }
    }

  FourMomentum get_tau_neutrino_mom(const Particle& p)  {
      assert(p.abspid() == PID::TAU);
      ConstGenVertexPtr dv = p.genParticle()->end_vertex();
      assert(dv != nullptr);
      for(ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){
        if (abs(pp->pdg_id()) == PID::NU_TAU) return FourMomentum(pp->momentum());
      }
      return FourMomentum();
    }

  void get_prong_number(ConstGenParticlePtr p, unsigned int& nprong, bool& lep_decaying_tau) {
      assert(p != nullptr);
      //const int tau_barcode = p->barcode();
      ConstGenVertexPtr dv = p->end_vertex();
      assert(dv != nullptr);
      for(ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){
        // If they have status 1 and are charged they will produce a track and the prong number is +1
        if (pp->status() == 1 )  {
          const int id = pp->pdg_id();
          if (Rivet::PID::charge(id) != 0 ) ++nprong;
          // Check if tau decays leptonically
          // @todo Can a tau decay include a tau in its decay daughters?!
          if ((abs(id) == PID::ELECTRON || abs(id) == PID::MUON || abs(id) == PID::TAU) && abs(p->pdg_id()) == PID::TAU) lep_decaying_tau = true;
        }
        // If the status of the daughter particle is 2 it is unstable and the further decays are checked
        else if (pp->status() == 2 )  {
          get_prong_number(pp, nprong, lep_decaying_tau);
        }
      }
    }

  }; //closes the majorana class at the very top

  //=================================  Plugin hook  =================================//
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(minimal_2l2j);


}  //closes the Rivet namespace
