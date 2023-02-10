//=================================  Introduction  =================================//<
// -*- C++ -*-
// Test routine - will produce dummy file with two trees
// If test fails - seg fault, malloc error, corruption etc... This is due to an error
// in setting up enviroment


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
  class TEST : public Analysis {

    public:

      /// Constructor
      TEST() : Analysis("TEST") {  
        std::cout << "TEST Constructor called" << std::endl;
      }

      //=================================  Declaring global variables  =================================//

      //Event counter
      int event_number = 0;
      int events_survive = 0;
      bool has_finalized = false;
      bool has_init = false;

      // Output file
      
      TFile *m_OutputFile;
      TTree *m_OutputTree;
      TTree *m_MetaDataTree;

      //=================================  Declaring global functions  =================================//

      //====================================================================================//
      //=================================  Initialisation  =================================//
      //====================================================================================//

      void init() {
        
        std::cout << "Entering init" << std::endl;

        if (!has_init){
          // Output file
          has_init = true;
          m_OutputFile = new TFile("output.root", "RECREATE");
          m_OutputTree = new TTree("tree", "tree");
          m_MetaDataTree = new TTree("metadata", "metadata");
          m_OutputTree->Branch("counter",           &event_number,   "event_number/I");
          m_MetaDataTree->Branch("events_survive",  &events_survive, "events_survive/I");
        }
      }
      //===============================================================================//
      //=================================  Analysis  =================================//
      //==============================================================================//

      /// Perform the per-event analysis
      //select particles, filter according to cuts, loop over combinations, construct observables, fill histograms.
      // This is where the per-event aspect of the analysis algorithm goes.

      ///Do the event by event analysis here (loop over every event)
      void analyze(const Event& event) {
        
        std::cout << "Entering Analyze" << std::endl;

        event_number++;

        if(event_number % 1000 == 0){std::cout << "Proccessing event: " << event_number << std::endl;}

        m_OutputFile->cd();
        m_OutputTree->Fill();
        
        events_survive++;

        std::cout << "Filled output tree" << std::endl;
    }

    //================================================================================//
    //=================================  Finalising  =================================//
    //================================================================================//

    /// Normalise histograms, scale/divide histograms etc., after the run
    void finalize() {
      
      std::cout << "Entering Finalize" << std::endl;
      
      if (!has_finalized){
        has_finalized = true;
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
    
  }; //closes the routine class at the very top

  //=================================  Plugin hook  =================================//
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TEST);


}  //closes the Rivet namespace
