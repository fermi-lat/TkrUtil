// $Header$

// Include files
// Gaudi system includes

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "TkrUtil/TkrFailureModeSvc.h"

// Define the class here instead of in a header file: not needed anywhere but here!
//------------------------------------------------------------------------------
/** 
A simple algorithm.
*/

class test_TkrUtil : public Algorithm {
public:
    test_TkrUtil(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private: 

  //! number of times called

  int m_count; 

  /// pointer to failure mode service
  ITkrFailureModeSvc* m_FailSvc;

};

//------------------------------------------------------------------------

// necessary to define a Factory for this algorithm
// expect that the xxx_load.cxx file contains a call     
//     DLL_DECL_ALGORITHM( test_TkrUtil );

static const AlgFactory<test_TkrUtil>  Factory;
const IAlgFactory& test_TkrUtilFactory = Factory;

//------------------------------------------------------------------------
//! ctor

test_TkrUtil::test_TkrUtil(const std::string& name, ISvcLocator* pSvcLocator)
:Algorithm(name, pSvcLocator)
,m_count(0)
{

}



//------------------------------------------------------------------------

//! set parameters and attach to various perhaps useful services.

StatusCode test_TkrUtil::initialize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;

    sc = service("TkrFailureModeSvc", m_FailSvc);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to find TkrFailureMode service" << endreq;
        return sc;
    }

    return sc;
}

//------------------------------------------------------------------------

//! process an event
StatusCode test_TkrUtil::execute()
{
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << " " << endreq;

    const int tower1 = 10;
    const int tower2 = 11;

    const int tower3 = 6;

    const int tower4 = 5;
    const int layer1 = 3;
    const int view1 = 0;
    const int layer2 = 4;
    const int view2 = 1;

    if (m_FailSvc == 0) return StatusCode::FAILURE;

    if (m_FailSvc->isFailed(tower1)) {
      log << MSG::INFO << "removed tower " << tower1 << endreq;
    } else {
      log << MSG::ERROR << "failed to remove tower " << tower1<< endreq;
      sc = StatusCode::FAILURE;
    }

    if (m_FailSvc->isFailed(tower2)) {
      log << MSG::INFO << "removed tower " << tower2 << endreq;
    } else {
      log << MSG::ERROR << "failed to remove tower " << tower2<< endreq;
      sc = StatusCode::FAILURE;
    }
    if (!m_FailSvc->isFailed(tower3)) {
      log << MSG::INFO << "correctly left tower " << tower3 << endreq;
    } else {
      log << MSG::ERROR << "erroneously removed tower " << tower3 << endreq;
      sc = StatusCode::FAILURE;
    }

    if (m_FailSvc->isFailed(tower4, layer1, view1)) {
      log << MSG::INFO << "removed tower " << tower4 << " layer " << layer1 << " view " << view1 << endreq;
    } else {
      log << MSG::ERROR << "failed to remove tower " << tower4 << " layer " << layer1 << " view " << view1 << endreq;
      sc = StatusCode::FAILURE;
    }

    if (m_FailSvc->isFailed(tower4, layer2, view2)) {
      log << MSG::INFO << "removed tower " << tower4 << " layer " << layer2 << " view " << view2 << endreq;
    } else {
      log << MSG::ERROR << "failed to remove tower " << tower4 << " layer " << layer2 << " view " << view2 << endreq;
      sc = StatusCode::FAILURE;
    }

    if (!m_FailSvc->isFailed(tower3, layer2, view2)) {
      log << MSG::INFO << "correctly left tower " << tower3<< " layer " << layer2 << " view " << view2 << endreq;
    } else {
      log << MSG::ERROR << "erroneously removed tower " << tower4 << " layer " << layer2 << " view " << view2 << endreq;
      sc = StatusCode::FAILURE;
    }

    log << MSG::INFO << " " << endreq;

    return sc;
}



//------------------------------------------------------------------------

//! clean up, summarize
StatusCode test_TkrUtil::finalize(){
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize after " << m_count << " calls." << endreq;

    return sc;
}







