#ifndef computeAngles_h
#define computeAngles_h

/*
 *  MELA - cf. http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/sbologne/MELAproject/
 *
 *  $Date: 2012/10/04 13:51:21 $
 *  $Revision: 1.5 $
 */

#include "TLorentzVector.h"

namespace mela {
  /// Leptons have to be massless for ME calculations.
  /// Remove lepton mass with constraint on mV if the flag is set to true
  extern bool forbidMassiveLeptons;
  void applyLeptonMassCorrection(bool flag=false);
  void constrainedRemoveLeptonMass(TLorentzVector& p1, TLorentzVector& p2);

  /// Compute decay angles from the lepton four-vectors and pdgIds.  
  /// Theta1 is the angle corresponding to Z1.
  /// Z1_lept1 and  Z1_lept2 are supposed to come from the same Z.
  /// Leptons are re-ordered internally according to a standard convention:
  /// lept1 = negative-charged lepton (for OS pairs).
  void computeAngles(TLorentzVector Z1_lept1, int Z1_lept1Id,
		     TLorentzVector Z1_lept2, int Z1_lept2Id,
		     TLorentzVector Z2_lept1, int Z2_lept1Id,
		     TLorentzVector Z2_lept2, int Z2_lept2Id,
		     float& costhetastar, 
		     float& costheta1, 
		     float& costheta2, 
		     float& Phi, 
		     float& Phi1);

  void computeAnglesCS(TLorentzVector Z1_lept1, int Z1_lept1Id,
			 TLorentzVector Z1_lept2, int Z1_lept2Id,
			 TLorentzVector Z2_lept1, int Z2_lept1Id,
			 TLorentzVector Z2_lept2, int Z2_lept2Id,
			 float pbeam,  
			 float& costhetastar, 
			 float& costheta1, 
			 float& costheta2, 
			 float& Phi, 
			 float& Phi1);

  /// Jets
  void computeJetMassless(TLorentzVector massiveJet, TLorentzVector& masslessJet); // Input massive -> output massless
  void computeFakeJet(TLorentzVector realJet, TLorentzVector others, TLorentzVector& fakeJet); // Input massive + higgs -> output massless fake jet
}
#endif
