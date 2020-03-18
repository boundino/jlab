#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"//calculate trajectory distance

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include "RecoVertex/TrimmedKalmanVertexFinder/interface/KalmanTrimmedVertexFinder.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"//proper covariance error calculation
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "Muonset/Muonset/interface/format.h"

class Muonset : public edm::EDAnalyzer
{
public:
  explicit Muonset(const edm::ParameterSet&);
  ~Muonset();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  edm::ParameterSet theConfig;
  edm::EDGetTokenT< std::vector<pat::Muon> > muonLabel_;
  edm::EDGetTokenT< reco::BeamSpot > bsLabel_;
  edm::EDGetTokenT< reco::VertexCollection > pvLabel_;

  edm::Service<TFileService> fs;
  TTree *root;
  MuonInfoBranches    MuonInfo;

};

void Muonset::beginJob()
{
  root = fs->make<TTree>("root","root");
  MuonInfo.regTree(root);
}

Muonset::Muonset(const edm::ParameterSet& iConfig):theConfig(iConfig)
{
  muonLabel_      = consumes< std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("MuonLabel"));
  bsLabel_        = consumes< reco::BeamSpot >(iConfig.getParameter<edm::InputTag>("BSLabel"));
  pvLabel_        = consumes< reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("PVLabel"));
}

Muonset::~Muonset()
{
}

void Muonset::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonLabel_, muons);
  std::vector<pat::Muon> input_muons = *muons;

  memset(&MuonInfo    ,0x00,sizeof(MuonInfo)  );

  MuonInfo.RunNo   = iEvent.id().run();
  MuonInfo.EvtNo   = iEvent.id().event();
  MuonInfo.LumiNo  = iEvent.luminosityBlock();

  Vertex thePrimaryV;
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(bsLabel_, beamSpotHandle);
  edm::Handle<reco::VertexCollection> VertexHandle;
  iEvent.getByToken(pvLabel_, VertexHandle);
  double PVBS_Pt_Max = -100.;
  if( VertexHandle.isValid() && !VertexHandle.failedToGet() && VertexHandle->size() > 0) 
    {
      if(int(VertexHandle->end()-VertexHandle->begin())==1)
        {
          thePrimaryV = *(VertexHandle->begin());
        }
      else
        {
          for(std::vector<reco::Vertex>::const_iterator it_vtx = VertexHandle->begin();it_vtx != VertexHandle->end(); it_vtx++ ) 
            {
              double Pt_Sum = 0;
              for(reco::Vertex::trackRef_iterator it = it_vtx->tracks_begin(); it != it_vtx->tracks_end(); it++) { Pt_Sum += (*it)->pt(); }
              if(Pt_Sum >= PVBS_Pt_Max) { thePrimaryV = *it_vtx; }
            }
        }
    }
  else
    { 
      thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }
  math::XYZPoint RefVtx = thePrimaryV.position();

  MuonInfo.size   = 0;
  for(std::vector<pat::Muon>::const_iterator mu_it=input_muons.begin();
      mu_it != input_muons.end() ; mu_it++)
    {
      if(mu_it->pt() < 4) continue; //
      if(!(mu_it->isGlobalMuon() && mu_it->isTrackerMuon() &&
           mu_it->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
           && mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 
           && mu_it->innerTrack()->quality(reco::TrackBase::highPurity)
           && fabs(mu_it->innerTrack()->dxy(thePrimaryV.position())) < 0.3
           && fabs(mu_it->innerTrack()->dz(thePrimaryV.position())) < 20.)) continue;

      MuonInfo.charge         [MuonInfo.size] = mu_it->charge();
      MuonInfo.pt             [MuonInfo.size] = mu_it->pt();
      MuonInfo.eta            [MuonInfo.size] = mu_it->eta();
      MuonInfo.phi            [MuonInfo.size] = mu_it->phi();
      // MuonInfo.isTrackerMuon  [MuonInfo.size] = mu_it->isTrackerMuon();
      // MuonInfo.isGlobalMuon   [MuonInfo.size] = mu_it->isGlobalMuon();
      if(mu_it->innerTrack().isNonnull())
        {
          MuonInfo.i_nStripLayer           [MuonInfo.size] = mu_it->innerTrack()->hitPattern().stripLayersWithMeasurement();
          MuonInfo.i_nPixelLayer           [MuonInfo.size] = mu_it->innerTrack()->hitPattern().pixelLayersWithMeasurement();
          MuonInfo.dzPV                    [MuonInfo.size] = mu_it->track()->dz(RefVtx); // ==mu_it->innerTrack()->dxy(thePrimaryV.position());
          MuonInfo.dxyPV                   [MuonInfo.size] = mu_it->track()->dxy(RefVtx); // ==mu_it->innerTrack()->dz(thePrimaryV.position());
        }
      MuonInfo.size++;
    }
  std::cout<<"found \e[33;1m"<<MuonInfo.size<<"\e[0m muons"<<std::endl;
  root->Fill();

}

// ------------ method called once each job just after ending the event loop  ------------{{{
void Muonset::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void Muonset::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void Muonset::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void Muonset::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void Muonset::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Muonset::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}//}}}

//define this as a plug-in
DEFINE_FWK_MODULE(Muonset);
