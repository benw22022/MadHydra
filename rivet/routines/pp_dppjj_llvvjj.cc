//=================================  Introduction  =================================//<
// -*- C++ -*-
//Mohamed Aly and Ahmed Abdelmotteleb
//Majorana Neutrino Sensitivity at ATLAS
//Rivet analysis code that runs over our data [signal + background] and applies cuts

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

//==========================================================================//
//=================================  MAIN  =================================//
//==========================================================================//

namespace Rivet {
  

  //divided into constructor + 3 loops: init, analyze and finalize

  //==================================================================================//
  //=================================  Constructer ==================================//
  //=================================================================================//
  class pp_dpp_llvvjj : public Analysis {

    public:

      double calc_deltaR(const double& y1, const double& y2, const double& phi1, const double& phi2){

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

      std::vector<double> b_nu_pt;
      std::vector<double> b_nu_rap;
      std::vector<double> b_nu_phi;
      std::vector<double> b_nu_mass;
      std::vector<double> b_nu_energy;
      std::vector<double> b_nu_px;
      std::vector<double> b_nu_py;
      std::vector<double> b_nu_pz;
      std::vector<long>   b_nu_pid;

      double b_HTlep;
      double b_HTjets;
      double b_eTmiss;
      double b_meff;

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

      int SMEAR;
      bool MET_SMEAR;
      bool JET_SMEAR;
      bool MUON_SMEAR;

      bool _use_fiducial_lepton_efficiency;

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
        m_OutputFile = new TFile("output.root", "RECREATE");
        m_OutputTree = new TTree("tree", "tree");
        m_MetaDataTree = new TTree("metadata", "metadata");
        m_OutputFile->cd();

        m_OutputTree->Branch("jet_pt",     &b_jet_pt    );
        m_OutputTree->Branch("jet_rap",    &b_jet_rap   );
        m_OutputTree->Branch("jet_phi",    &b_jet_pt    );
        m_OutputTree->Branch("jet_mass",   &b_jet_mass  );
        m_OutputTree->Branch("jet_energy", &b_jet_energy );
        m_OutputTree->Branch("jet_px",     &b_jet_px    );
        m_OutputTree->Branch("jet_py",     &b_jet_py    );
        m_OutputTree->Branch("jet_pz",     &b_jet_pz    );

        m_OutputTree->Branch("lep_pt",     &b_lep_pt    );
        m_OutputTree->Branch("lep_rap",    &b_lep_rap   );
        m_OutputTree->Branch("lep_phi",    &b_lep_pt    );
        m_OutputTree->Branch("lep_mass",   &b_lep_mass  );
        m_OutputTree->Branch("lep_energy", &b_lep_energy);
        m_OutputTree->Branch("lep_px",     &b_lep_px    );
        m_OutputTree->Branch("lep_py",     &b_lep_py    );
        m_OutputTree->Branch("lep_pz",     &b_lep_pz    );
        m_OutputTree->Branch("lep_pid",    &b_lep_pid   );

        m_OutputTree->Branch("tau_pt",     &b_tau_pt    );
        m_OutputTree->Branch("tau_rap",    &b_tau_rap   );
        m_OutputTree->Branch("tau_phi",    &b_tau_pt    );
        m_OutputTree->Branch("tau_mass",   &b_tau_mass  );
        m_OutputTree->Branch("tau_energy", &b_tau_energy);
        m_OutputTree->Branch("tau_px",     &b_tau_px    );
        m_OutputTree->Branch("tau_py",     &b_tau_py    );
        m_OutputTree->Branch("tau_pz",     &b_tau_pz    );
        m_OutputTree->Branch("tau_pid",    &b_tau_pid   );

        m_OutputTree->Branch("eTmiss",     &b_eTmiss,  "eTmiss/D");
        m_OutputTree->Branch("HTlep",      &b_HTlep,   "HTlep/D");
        m_OutputTree->Branch("HTjets",     &b_HTjets,  "HTjets/D");
        m_OutputTree->Branch("meff",       &b_meff,    "meff/D");
        
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



        // To calculate the acceptance without having the fiducial lepton efficiencies included, this part can be turned off
        _use_fiducial_lepton_efficiency = true;

        // Random numbers for simulation of ATLAS detector reconstruction efficiency
        srand(160385);

        // Final state including all charged and neutral particles
        const FinalState fs((Cuts::etaIn(-5.0, 5.0) && Cuts::pT >=  1*GeV));
        declare(fs, "FS");

        // Final state including all charged particles
        declare(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 1*GeV), "CFS");

        // Final state including all visible particles (to calculate MET, Jets etc.)
        declare(VisibleFinalState(Cuts::abseta < 5.0), "VFS");

        // Final state including all AntiKt 04 Jets
        VetoedFinalState vfs;
        vfs.addVetoPairId(PID::MUON);
        declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

        // Final state including all unstable particles (including taus)
        declare(UnstableParticles(Cuts::abseta < 5.0 && Cuts::pT > 5*GeV), "UFS");

        // Final state including all electrons
        IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 10*GeV);
        elecs.acceptIdPair(PID::ELECTRON);
        declare(elecs, "elecs");

        // Final state including all muons
        IdentifiedFinalState muons(Cuts::abseta < 2.5 && Cuts::pT > 10*GeV);
        muons.acceptIdPair(PID::MUON);
        declare(muons, "muons");

        //=================================  Photon and Neutrino Decleration  =================================//
        // define photon
        IdentifiedFinalState photon(fs);
        photon.acceptIdPair(PID::PHOTON); 
        declare(photon,"photon");

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

        if(event_number % 1000 == 0){std::cout << "Proccessing event: " << event_number << std::endl;}


        
        //=================== weight and variance ===================//

        //This will calculate the weight for each of the N events that rivet processes

        //weight of each event (Rivet has for some reason removed the weight() method, the 0th element should be the default weight from looking at the output of checkMetaSG.py)
        const double weight = event.weights()[0];

        //summing over all the weights (sum of weights)
        sumW += weight;

        //Sum up the weights^2 --> This happens to be equal to the variance
        var_sumW += std::pow(weight, 2);

        //================================= Building Candidate Leptons  -- Dressed =================================//
        // Muons
        Particles muon_candidates;
        const Particles charged_tracks    = apply<ChargedFinalState>(event, "CFS").particles();
        const Particles visible_particles = apply<VisibleFinalState>(event, "VFS").particles();
        for (const Particle& mu : apply<IdentifiedFinalState>(event, "muons").particlesByPt()) {
          // Calculate pTCone30 variable (pT of all tracks within dR<0.3 - pT of muon itself)
          double pTinCone = -mu.pT();
          for (const Particle& track : charged_tracks) {
            if (deltaR(mu.momentum(), track.momentum()) < 0.3)
              pTinCone += track.pT();
          }

          // Calculate eTCone30 variable (pT of all visible particles within dR<0.3)
          double eTinCone = 0.;
          for (const Particle& visible_particle : visible_particles) {
            if (visible_particle.abspid() != PID::MUON && inRange(deltaR(mu.momentum(), visible_particle.momentum()), 0.1, 0.3))
              eTinCone += visible_particle.pT();
          }

          // Apply reconstruction efficiency and simulate reco
          int muon_id = 13;
          if ( mu.hasAncestor(15) || mu.hasAncestor(-15)) muon_id = 14;
          const double eff = (_use_fiducial_lepton_efficiency) ? apply_reco_eff(muon_id, mu) : 1.0;
          const bool keep_muon = rand()/static_cast<double>(RAND_MAX) <= eff;

          // Keep muon if pTCone30/pT < 0.15 and eTCone30/pT < 0.2 and reconstructed
          if (keep_muon && pTinCone/mu.pT() <= 0.15 && eTinCone/mu.pT() < 0.2)
            muon_candidates.push_back(mu);
        }


        // Electrons
        Particles electron_candidates;
        for (const Particle& e : apply<IdentifiedFinalState>(event, "elecs").particlesByPt()) {
          // Neglect electrons in crack regions
          if (inRange(e.abseta(), 1.37, 1.52)) continue;

          // Calculate pTCone30 variable (pT of all tracks within dR<0.3 - pT of electron itself)
          double pTinCone = -e.pT();
          for (const Particle& track : charged_tracks) {
            if (deltaR(e.momentum(), track.momentum()) < 0.3) pTinCone += track.pT();
          }

          // Calculate eTCone30 variable (pT of all visible particles (except muons) within dR<0.3)
          double eTinCone = 0.;
          for (const Particle& visible_particle : visible_particles) {
            if (visible_particle.abspid() != PID::MUON && inRange(deltaR(e.momentum(), visible_particle.momentum()), 0.1, 0.3))
              eTinCone += visible_particle.pT();
          }

          // Apply reconstruction efficiency and simulate reco
          int elec_id = 11;
          if (e.hasAncestor(15) || e.hasAncestor(-15)) elec_id = 12;
          const double eff = (_use_fiducial_lepton_efficiency) ? apply_reco_eff(elec_id, e) : 1.0;
          const bool keep_elec = rand()/static_cast<double>(RAND_MAX) <= eff;

          // Keep electron if pTCone30/pT < 0.13 and eTCone30/pT < 0.2 and reconstructed
          if (keep_elec && pTinCone/e.pT() <= 0.13 && eTinCone/e.pT() < 0.2)
            electron_candidates.push_back(e);
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
          // keep only taus in certain eta region and above 15 GeV of visible tau pT
          if ( tau_vis.pT() <= 15.0*GeV || tau_vis.abseta() > 2.5) continue;

          // Get prong number (number of tracks) in tau decay and check if tau decays leptonically
          unsigned int nprong = 0;
          bool lep_decaying_tau = false;
          get_prong_number(tau.genParticle(), nprong, lep_decaying_tau);

          // Apply reconstruction efficiency
          int tau_id = 15;
          if (nprong == 1) tau_id = 15;
          else if (nprong == 3) tau_id = 16;

          // Get fiducial lepton efficiency simulate reco efficiency
          const double eff = (_use_fiducial_lepton_efficiency) ? apply_reco_eff(tau_id, tau_vis) : 1.0;
          const bool keep_tau = rand()/static_cast<double>(RAND_MAX) <= eff;

          // Keep tau if nprong = 1, it decays hadronically, and it's reconstructed by the detector
          if ( !lep_decaying_tau && nprong == 1 && keep_tau) tau_candidates.push_back(tau_vis);
        }


        // Jets (all anti-kt R=0.4 jets with pT > 25 GeV and eta < 4.9)
        Jets jet_candidates;
        for (const Jet& jet : apply<FastJets>(event, "AntiKtJets04").jetsByPt(25*GeV)) {
          if (jet.abseta() < 4.9) jet_candidates.push_back(jet);
        }


        // ETmiss
        Particles vfs_particles = apply<VisibleFinalState>(event, "VFS").particles();
        FourMomentum pTmiss;
        for (const Particle& p : vfs_particles) pTmiss -= p.momentum();
        double eTmiss = pTmiss.pT()/GeV;


        //------------------
        // Overlap removal

        // electron - electron
        Particles electron_candidates_2;
        for (size_t ie = 0; ie < electron_candidates.size(); ++ie) {
          const Particle & e = electron_candidates[ie];
          bool away = true;
          // If electron pair within dR < 0.1: remove electron with lower pT
          for (size_t ie2=0; ie2 < electron_candidates_2.size(); ++ie2) {
            if ( deltaR( e.momentum(), electron_candidates_2[ie2].momentum()) < 0.1 ) {
              away = false;
              break;
            }
          }
          // If isolated keep it
          if ( away )
            electron_candidates_2.push_back( e );
        }
        // jet - electron
        Jets recon_jets;
        for (const Jet& jet : jet_candidates) {
          bool away = true;
          // if jet within dR < 0.2 of electron: remove jet
          for (const Particle& e : electron_candidates_2) {
            if (deltaR(e.momentum(), jet.momentum()) < 0.2) {
              away = false;
              break;
            }
          }
          // jet - tau
          if (away)  {
            // If jet within dR < 0.2 of tau: remove jet
            for (const Particle& tau : tau_candidates) {
              if (deltaR(tau.momentum(), jet.momentum()) < 0.2) {
                away = false;
                break;
              }
            }
          }
          // If isolated keep it
          if ( away )
            recon_jets.push_back( jet );
        }


        // electron - jet
        Particles recon_leptons, recon_e;
        for (size_t ie = 0; ie < electron_candidates_2.size(); ++ie) {
          const Particle& e = electron_candidates_2[ie];
          // If electron within 0.2 < dR < 0.4 from any jets: remove electron
          bool away = true;
          for (const Jet& jet : recon_jets) {
            if (deltaR(e.momentum(), jet.momentum()) < 0.4) {
              away = false;
              break;
            }
          }
          // electron - muon
          // if electron within dR < 0.1 of a muon: remove electron
          if (away) {
            for (const Particle& mu : muon_candidates) {
              if (deltaR(mu.momentum(), e.momentum()) < 0.1) {
                away = false;
                break;
              }
            }
          }
          // If isolated keep it
          if (away)  {
            recon_e += e;
            recon_leptons += e;
          }
        }


        // tau - electron
        Particles recon_tau;
        for ( const Particle& tau : tau_candidates ) {
          bool away = true;
          // If tau within dR < 0.2 of an electron: remove tau
          for ( const Particle& e : recon_e ) {
            if (deltaR( tau.momentum(), e.momentum()) < 0.2) {
              away = false;
              break;
            }
          }
          // tau - muon
          // If tau within dR < 0.2 of a muon: remove tau
          if (away)  {
            for (const Particle& mu : muon_candidates) {
              if (deltaR(tau.momentum(), mu.momentum()) < 0.2) {
                away = false;
                break;
              }
            }
          }
          // If isolated keep it
          if (away) recon_tau.push_back( tau );
        }

        // Muon - jet isolation
        Particles recon_mu, trigger_mu;
        // If muon within dR < 0.4 of a jet, remove muon
        for (const Particle& mu : muon_candidates) {
          bool away = true;
          for (const Jet& jet : recon_jets) {
            if ( deltaR( mu.momentum(), jet.momentum()) < 0.4 ) {
              away = false;
              break;
            }
          }
          if (away) {
            recon_mu.push_back( mu );
            recon_leptons.push_back( mu );
            if (mu.abseta() < 2.4) trigger_mu.push_back( mu );
          }
        }

        // End overlap removal
        //------------------


        // Jet cleaning
        if (rand()/static_cast<double>(RAND_MAX) <= 0.42) {
          for (const Jet& jet : recon_jets) {
            const double eta = jet.rapidity();
            const double phi = jet.azimuthalAngle(MINUSPI_PLUSPI);
            if (jet.pT() > 25*GeV && inRange(eta, -0.1, 1.5) && inRange(phi, -0.9, -0.5)) vetoEvent;
          }
        }


        // Post-isolation event cuts
        // Require at least 1 charged tracks in event
        if (charged_tracks.size() < 1) vetoEvent;

        // And at least one e/mu passing trigger
        if (!( !recon_e   .empty() && recon_e[0]   .pT() > 25*GeV)  &&
            !( !trigger_mu.empty() && trigger_mu[0].pT() > 25*GeV) ) {
          MSG_DEBUG("Hardest lepton fails trigger");
          vetoEvent;
        }

        // And only accept events with at least 2 lepton
        if (recon_leptons.size() < 2) vetoEvent;

        // Sort leptons by decreasing pT
        sortByPt(recon_leptons);
        sortByPt(recon_tau);

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

        for(const Particle &p : recon_tau){
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

      for(const Jet &p : recon_jets){
        // std::cout << "Jet pT  = " << p.pT() << std::endl;
        
        // float* p_pT = new float(*p.pT());

        // b_jet_pt.push_back(*p_pT);
        b_jet_pt.push_back(p.pT());
        b_jet_rap.push_back(p.rapidity());
        b_jet_phi.push_back(p.phi());
        b_jet_mass.push_back(p.mass());
        b_jet_energy.push_back(p.energy());
        b_jet_px.push_back(p.px() );
        b_jet_py.push_back(p.py());
        b_jet_pz.push_back(p.pz());

        // delete p_pT;
      }

      for(unsigned int i{0}; i<b_jet_pt.size(); i++){
        std::cout << "b_jet_pt[" << i << "] = " << b_jet_pt[i] << std::endl;
      }

      // for(const Particle &p : neutrinos){
      //   b_nu_pt.push_back(p.pT());
      //   b_nu_rap.push_back(p.rapidity());
      //   b_nu_phi.push_back(p.phi());
      //   b_nu_mass.push_back(p.mass());
      //   b_nu_energy.push_back(p.energy());
      //   b_nu_px.push_back(p.px() );
      //   b_nu_py.push_back(p.py());
      //   b_nu_pz.push_back(p.pz());
      //   b_nu_pid.push_back(p.pid());
      // }

      // Calculate HTjets
      double HTjets = 0.;
      for ( const Jet & jet : recon_jets )
        HTjets += jet.perp()/GeV;


      // Calculate meff
      double meff = eTmiss + HTjets;
      for ( const Particle & e  : recon_e  )  {
        meff += e.perp()/GeV;
      }
      for ( const Particle & mu : recon_mu )  {
        meff += mu.perp()/GeV;
      }
      for ( const Particle & tau : recon_tau )  {
        meff += tau.perp()/GeV;
      }

      b_HTlep = HTlep;
      b_HTjets = HTjets;
      b_eTmiss = eTmiss;
      b_meff = meff;

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
      
      b_nu_pt.clear();
      b_nu_rap.clear();
      b_nu_phi.clear();
      b_nu_mass.clear();
      b_nu_energy.clear();
      b_nu_px.clear();
      b_nu_py.clear();
      b_nu_pz.clear();
      b_nu_pid.clear();

      // m_OutputTree->Print();
      std::cout << "Filled output tree" << std::endl;
    }

    //================================================================================//
    //=================================  Finalising  =================================//
    //================================================================================//

    /// Normalise histograms, scale/divide histograms etc., after the run
    void finalize() {
      
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
        b_presel_eff              =  survive_event_number / event_number * 100;   

      
        std::cout << "Filling metadata tree" << std::endl;
        m_MetaDataTree->Fill();
        std::cout << "Filled metadata tree" << std::endl;
        
        std::cout << "Writing MetaDataTree" << std::endl;
        m_MetaDataTree->Write();
        std::cout << "Written MetaDataTree" << std::endl;
        
        std::cout << "Writing OutputTree" << std::endl;
        m_OutputTree->Write();
        std::cout << "Written OutputTree" << std::endl;

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


    /// Function giving fiducial lepton efficiency
    double apply_reco_eff(int flavor, const Particle& p) {
      float pt = p.pT()/GeV;
      float eta = p.eta();

      double eff = 0.;
      //double err = 0.;

      if (flavor == 11) { // weight prompt electron -- now including data/MC ID SF in eff.
        //float rho = 0.820;
        float p0 = 7.34;  float p1 = 0.8977;
        //float ep0= 0.5 ;  float ep1= 0.0087;
        eff = p1 - p0/pt;

        //double err0 = ep0/pt; // d(eff)/dp0
        //double err1 = ep1;    // d(eff)/dp1
        //err = sqrt(err0*err0 + err1*err1 - 2*rho*err0*err1);

        double avgrate = 0.6867;
        float wz_ele_eta[] = {0.588717,0.603674,0.666135,0.747493,0.762202,0.675051,0.751606,0.745569,0.665333,0.610432,0.592693,};
        //float ewz_ele_eta[] ={0.00292902,0.002476,0.00241209,0.00182319,0.00194339,0.00299785,0.00197339,0.00182004,0.00241793,0.00245997,0.00290394,};
        int ibin = 3;

        if (eta >= -2.5 && eta < -2.0) ibin = 0;
        if (eta >= -2.0 && eta < -1.5) ibin = 1;
        if (eta >= -1.5 && eta < -1.0) ibin = 2;
        if (eta >= -1.0 && eta < -0.5) ibin = 3;
        if (eta >= -0.5 && eta < -0.1) ibin = 4;
        if (eta >= -0.1 && eta <  0.1) ibin = 5;
        if (eta >=  0.1 && eta <  0.5) ibin = 6;
        if (eta >=  0.5 && eta <  1.0) ibin = 7;
        if (eta >=  1.0 && eta <  1.5) ibin = 8;
        if (eta >=  1.5 && eta <  2.0) ibin = 9;
        if (eta >=  2.0 && eta <  2.5) ibin = 10;

        double eff_eta = wz_ele_eta[ibin];
        //double err_eta = ewz_ele_eta[ibin];

        eff = (eff*eff_eta)/avgrate;
      }

      if (flavor == 12)  { // weight electron from tau
        //float rho = 0.884;
        float p0 = 6.799;  float p1 = 0.842;
        //float ep0= 0.664;  float ep1= 0.016;
        eff = p1 - p0/pt;

        //double err0 = ep0/pt; // d(eff)/dp0
        //double err1 = ep1;    // d(eff)/dp1
        //err = sqrt(err0*err0 + err1*err1 - 2*rho*err0*err1);

        double avgrate = 0.5319;
        float wz_elet_eta[] = {0.468945,0.465953,0.489545,0.58709,0.59669,0.515829,0.59284,0.575828,0.498181,0.463536,0.481738,};
        //float ewz_elet_eta[] ={0.00933795,0.00780868,0.00792679,0.00642083,0.00692652,0.0101568,0.00698452,0.00643524,0.0080002,0.00776238,0.0094699,};
        int ibin = 3;

        if (eta >= -2.5 && eta < -2.0) ibin = 0;
        if (eta >= -2.0 && eta < -1.5) ibin = 1;
        if (eta >= -1.5 && eta < -1.0) ibin = 2;
        if (eta >= -1.0 && eta < -0.5) ibin = 3;
        if (eta >= -0.5 && eta < -0.1) ibin = 4;
        if (eta >= -0.1 && eta <  0.1) ibin = 5;
        if (eta >=  0.1 && eta <  0.5) ibin = 6;
        if (eta >=  0.5 && eta <  1.0) ibin = 7;
        if (eta >=  1.0 && eta <  1.5) ibin = 8;
        if (eta >=  1.5 && eta <  2.0) ibin = 9;
        if (eta >=  2.0 && eta <  2.5) ibin = 10;

        double eff_eta = wz_elet_eta[ibin];
        //double err_eta = ewz_elet_eta[ibin];

        eff = (eff*eff_eta)/avgrate;

      }

      if (flavor == 13)  {// weight prompt muon

        //if eta>0.1
        float p0 = -18.21;  float p1 = 14.83;  float p2 = 0.9312;
        //float ep0= 5.06;    float ep1= 1.9;    float ep2=0.00069;

        if ( fabs(eta) < 0.1)  {
          p0  = 7.459; p1 = 2.615; p2  = 0.5138;
          //ep0 = 10.4; ep1 = 4.934; ep2 = 0.0034;
        }

        double arg = ( pt-p0 )/( 2.*p1 ) ;
        eff = 0.5 * p2 * (1.+erf(arg));
        //err = 0.1*eff;
      }

      if (flavor == 14)  {// weight muon from tau

        if (fabs(eta) < 0.1) {
          float p0 = -1.756;  float p1 = 12.38;  float p2 = 0.4441;
          //float ep0= 10.39;   float ep1= 7.9;  float ep2=0.022;
          double arg = ( pt-p0 )/( 2.*p1 ) ;
          eff = 0.5 * p2 * (1.+erf(arg));
          //err = 0.1*eff;
        }
        else {
          float p0 = 2.102;  float p1 = 0.8293;
          //float ep0= 0.271;  float ep1= 0.0083;
          eff = p1 - p0/pt;
          //double err0 = ep0/pt; // d(eff)/dp0
          //double err1 = ep1;    // d(eff)/dp1
          //err = sqrt(err0*err0 + err1*err1 - 2*rho*err0*err1);
        }
      }

      if (flavor == 15)  {// weight hadronic tau 1p

        float wz_tau1p[] = {0.0249278,0.146978,0.225049,0.229212,0.21519,0.206152,0.201559,0.197917,0.209249,0.228336,0.193548,};
        //float ewz_tau1p[] ={0.00178577,0.00425252,0.00535052,0.00592126,0.00484684,0.00612941,0.00792099,0.0083006,0.0138307,0.015568,0.0501751,};
        int ibin = 0;
        if (pt > 15)  ibin = 1;
        if (pt > 20)  ibin = 2;
        if (pt > 25)  ibin = 3;
        if (pt > 30)  ibin = 4;
        if (pt > 40)  ibin = 5;
        if (pt > 50)  ibin = 6;
        if (pt > 60)  ibin = 7;
        if (pt > 80)  ibin = 8;
        if (pt > 100) ibin = 9;
        if (pt > 200) ibin = 10;

        eff = wz_tau1p[ibin];
        //err = ewz_tau1p[ibin];


        double avgrate = 0.1718;
        float wz_tau1p_eta[] = {0.162132,0.176393,0.139619,0.178813,0.185144,0.210027,0.203937,0.178688,0.137034,0.164216,0.163713,};
        //float ewz_tau1p_eta[] ={0.00706705,0.00617989,0.00506798,0.00525172,0.00581865,0.00865675,0.00599245,0.00529877,0.00506368,0.00617025,0.00726219,};

        ibin = 3;
        if (eta >= -2.5 && eta < -2.0) ibin = 0;
        if (eta >= -2.0 && eta < -1.5) ibin = 1;
        if (eta >= -1.5 && eta < -1.0) ibin = 2;
        if (eta >= -1.0 && eta < -0.5) ibin = 3;
        if (eta >= -0.5 && eta < -0.1) ibin = 4;
        if (eta >= -0.1 && eta <  0.1) ibin = 5;
        if (eta >=  0.1 && eta <  0.5) ibin = 6;
        if (eta >=  0.5 && eta <  1.0) ibin = 7;
        if (eta >=  1.0 && eta <  1.5) ibin = 8;
        if (eta >=  1.5 && eta <  2.0) ibin = 9;
        if (eta >=  2.0 && eta <  2.5) ibin = 10;

        double eff_eta = wz_tau1p_eta[ibin];
        //double err_eta = ewz_tau1p_eta[ibin];

        eff = (eff*eff_eta)/avgrate;
      }

      if (flavor == 16)  { //weight hadronic tau 3p

        float wz_tau3p[] = {0.000587199,0.00247181,0.0013031,0.00280112,};
        //float ewz_tau3p[] ={0.000415091,0.000617187,0.000582385,0.00197792,};

        int ibin = 0;
        if (pt > 15) ibin = 1;
        if (pt > 20) ibin = 2;
        if (pt > 40) ibin = 3;
        if (pt > 80) ibin = 4;

        eff = wz_tau3p[ibin];
        //err = ewz_tau3p[ibin];
      }

      return eff;
    }

  }; //closes the majorana class at the very top

  //=================================  Plugin hook  =================================//
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(pp_dpp_llvvjj);


}  //closes the Rivet namespace
