import FWCore.ParameterSet.Config as cms
from CalibPPS.ESProducers.ctppsCompositeESSource_cfi import ctppsCompositeESSource as _esComp
from CalibPPS.ESProducers.ppsAssociationCuts_non_DB_cff import use_single_infinite_iov_entry, p2021
from CalibPPS.ESProducers.ppsAssociationCuts_non_DB_cff import ppsAssociationCutsESSource as _esAssCuts
from Geometry.VeryForwardGeometry.commons_cff import cloneGeometry
from SimPPS.DirectSimProducer.profiles_2021_cff import profile_2021_default
from SimPPS.DirectSimProducer.simPPS2017_cfi import rpIds

ppsAssociationCutsESSource = _esAssCuts.clone()
use_single_infinite_iov_entry(ppsAssociationCutsESSource, p2021)
XMLIdealGeometryESSource_CTPPS, _ctppsGeometryESModule = cloneGeometry('Geometry.VeryForwardGeometry.geometryRPFromDD_2021_cfi')
# not cloning the ctppsGeometryESModule, as it is replaced by the composite ES source

ctppsCompositeESSource = _esComp.clone(
    generateEveryNEvents = 100,
    periods = [profile_2021_default],
    # geometry (using 2017 here is OK)
    compactViewTag = _ctppsGeometryESModule.compactViewTag,
    isRun2 = _ctppsGeometryESModule.isRun2
)
