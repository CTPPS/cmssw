import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from RecoPPS.Local.ctppsDiamondLocalReconstruction_cff import *

process = cms.Process("TIMINGSTUDY")
options = VarParsing ('analysis')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

options.register ('outputFileName',
		  'run_output.root',
		  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "output file name")

options.register ('calibInput',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Calibration input file for this iteration")

options.register ('planesConfig',
				  'planes.json',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "planes config file")

options.register ('geometryFile',
				  'Geometry.VeryForwardGeometry.geometryRPFromDD_2018_cfi',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Geometry input file")

options.register ('validOOT',
				  -1,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "valid OOT slice")				  
				  
options.parseArguments()
				  

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    *(
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/240000/1ECC0154-9688-6248-ADF5-7B9180CB7CCF.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/240000/21FEFBBE-0E7B-2A4D-A207-917AB49701A2.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/240000/6EE9461B-FC71-7144-94F8-A3AD82581D20.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/240000/C97D30B3-21DB-6447-9596-3220587A02B5.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/240000/D428E596-E313-2848-B326-5C4374865723.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/240000/DBD2A4DF-B43F-FA49-BD25-286E234FCA58.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/240000/DC146602-04E9-A94D-BD8A-8CE20ADA4996.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/240000/E1655977-576F-974F-B197-8898546015F4.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/240000/F57C9DFD-920B-994A-85AD-08983D85A3DA.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/2410000/080968F9-5C0E-BE4D-8D8B-D20AAB71C291.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/2410000/7D2AF654-F36D-894F-9922-EEEC29AEF1BA.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/2410000/FCC39027-5D9F-6942-9D94-9EDB65D41B77.root",
    "/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/002EF773-026C-2C49-981E-602DC3F43DF5.root",
    "/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/01BDECC6-EDA1-A141-94C7-2BC88B97A586.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/027CA2FB-8559-FF47-9026-E3A716F7569C.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/038AE40E-E498-0E46-AE33-E396FCA00136.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/04E49D12-050E-2243-AEE5-77EA9CEED5F7.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/054E7E4C-A75A-0E45-A4C7-07F328D9C0BB.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/07AC3C61-0EF7-A14F-AB4C-84A61E888B6C.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/08319D13-3375-2B44-BB09-43C91825EC7E.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/0BFE7589-6C55-D649-8DE4-54D18625BB18.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/0D2142AD-9C6D-2C48-B9AE-55185ED382FA.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/0FB29033-CEC1-F541-98AD-33365C270C26.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/12628C6B-81F4-134C-944E-22C0E536F4F4.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/13332FC7-A67A-EF40-9719-4D77B07C4FBC.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/1365ED78-1316-4041-AFB6-EB8DCFB4FE4F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/15A7D2D0-0A82-7749-B930-C10D065CE793.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/164B0F40-3659-CE4C-BFA2-6835AC225C3F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/1710300F-216D-BC44-8A1F-3B79D53ABE09.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/17C26920-7618-B74A-A261-007A1914376E.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/17F11341-84E8-0F4D-803A-7B0FC40D9623.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/1B4E20D4-587D-7942-B5C6-5737E62BDDDC.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/1BF800E8-9712-E349-9404-4DC4B35839F8.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/1C1D5637-360E-7F4A-AF2B-0C6E6BE03CA6.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/1D4E5A98-DD40-C140-A870-8D35611F46F4.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/1DBFE360-7E73-BA40-98BE-3851C984C19F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/1E9A0DED-6E9B-6F4D-999E-C84AAA0E4B45.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/1EA6651C-DA91-D446-9B21-2349DA0A00BF.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/21AD55BA-C4A5-F84D-8052-56FAB0339BC2.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/220227D8-1966-B24F-85D0-50D92049F881.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/222EA925-28A9-5B4D-BB94-33FE97983456.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/25013A84-C838-7F40-881F-BD3F82C79CF9.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/276328D7-5E7A-EF45-A8C6-3511B7358FE7.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/28A34A0B-2269-DD43-8FDC-3BEEE414D7C0.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/2936BDF5-0538-284E-87BD-9C733117640A.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/2A07C993-ADEA-5842-BD61-F270C46AF7E3.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/2A132A03-80DB-1D46-B875-FB7F8053220F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/2C8C5457-2AF6-D744-BEC0-3840C660BDE1.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/2CA805EF-FFE0-AE40-8A1F-84702A003F05.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/2D0E6263-426C-C84B-A087-2C169CD0738F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/30D8053C-877F-A74C-BA89-65C850639765.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/32C4F351-95D6-BD45-9228-2A5850DCDEEF.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/33B13A47-BA71-9B42-8362-321BA8CBC19D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/36492643-087E-6F47-8879-68A1A2387223.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/36BA6EBA-00FC-C743-AD4F-1FBEEA6F2FB4.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/38364BFE-7F84-CC4C-8545-0656BD19DBB9.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/38E4DDC6-686D-424E-B839-0C999A6ABE1D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/39FAD98A-8337-C044-B76B-68E0A1530BC4.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/3A0865E1-03FA-A441-A6EC-A91D8DDE9BED.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/3E094F5C-223A-4644-B618-350A30153F9E.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/3E98AAB2-FDFA-DD47-929B-F2A66BD66AE6.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/3F7367D9-8D26-D24D-BDCC-BA8F7504BA53.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/406D4CBF-BF41-DF4B-AB75-68D59D52EC32.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/42DBE057-ED41-4742-9AFA-DDEF8B5F5A86.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/43E5D3CA-6A8C-E044-AF25-EECF2755CC0A.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/4501C237-E353-2047-8F0E-EFB03CFE8C6A.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/45E5E0AE-30C1-3143-AE89-6521CC1618B8.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/4664D4C3-7B36-7A43-8475-78ADD908CB10.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/47049D58-A17C-CD45-A986-B8B237D6CA24.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/470D6BF5-7E1A-6541-BD08-132BC82BD978.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/47578029-D405-6B42-8F8D-B709323BAD15.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/47764251-509B-BD40-81BA-8E5BA70158CB.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/4899A728-08E3-5A44-A386-46AA46BD0B21.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/49B8E525-D0F8-3549-A0F6-45A328D7C9B5.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/4AAD90D1-5F97-9E42-8E65-F91EBC637C03.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/4AF3C43A-487E-B646-90D5-641964AF8965.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/4B5AA617-BB2B-104F-861E-D1321ED8F7F5.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/4BC3C2A1-AEF4-D343-AE37-8876C2C2928A.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/51CE0CD1-6597-F548-BD91-DCC1391C2023.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/5456BC8A-5E20-ED4B-B1C7-B97FC897202D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/5752A569-7537-D64C-B3B4-5D0370EE0120.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/59409B26-E339-CF47-B787-25F2BEC344A0.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/5CED1947-044B-6740-868B-65C4FF84B145.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/6533150E-5BA7-6E48-84A9-AC6AF7F8A7C3.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/68BA8E6C-8563-8248-93C3-B820A5ACB6E5.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/69B7D5EB-12C3-5C42-B38E-1722F6E06DB7.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/6D29D075-1873-BE42-A0B6-D8C1D86D1BC3.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/6D4ED335-A5F9-784D-92A4-6A9A315D4BAB.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/71A39B2D-2CA8-C24B-925B-D24B0DB79679.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/723674DE-5DC0-0B42-AFFF-83683A7D7985.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/729BCAE1-8005-5843-9463-9FE821FB8ED7.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/72EED3E8-2935-8846-BC0D-56087E6555C9.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/74D0935F-94EE-594A-9503-AE2B1343A8AB.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/7592CE6B-95EB-424D-ABDE-3B17CA8ED415.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/77027D93-4242-384E-AE90-3FC90BBFCB90.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/77D9A1C6-9C9D-2B4F-A806-244ED956B379.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/7AF0F37B-498A-854E-8006-D55978A3CFEE.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/7BDF6B5E-720F-FF4B-91F6-244ACB92FB6B.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/7C68F877-FC78-7040-ABF5-56D747098B4D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/7CC3687C-4A49-2A43-9639-27CE10A62552.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/7E8B4627-1B0B-844D-A8D4-A754FF005DF4.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/80038494-37AD-8643-89CC-9D9AE345CD1F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/82A27FD0-5985-8842-BE03-37D32D985517.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/849B01BA-7A24-8545-881D-F29C621606AE.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/84DE1CCD-571D-094A-A0FB-ED30CB83F23E.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/86CE0C6C-8EC3-434B-B528-D48B76CC4D4D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/89D1A307-C9E6-144A-956A-7F38B2085C69.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/92071577-519B-CC44-B266-3237A52093A5.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/93D853EC-5D44-C14A-BA6E-A52C1B7ED615.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/95F494ED-2199-524B-9BA0-6C267ABCF81D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/9738533D-1DB4-E140-9AB6-F287B1AA4757.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/97470E76-19E4-7449-AB42-4208BAE1D3E5.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/97996C9D-2645-8441-A55B-97E9BD23056D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/97F04906-C522-044E-A39A-99EB534A6564.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/9C1EF435-0F21-AF4B-94E7-2308C82AF846.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/9C3A1404-06D0-F647-85DB-480B3BEE6997.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/9CF44F16-BB98-2945-A01B-8FCFF93C5DE0.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/9CF45680-C6C3-7245-89B9-C40794149DDA.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/A3901447-064E-6746-B03B-8B8D9562F389.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/A41E98DF-AD86-EA45-AE4B-6EF7C509DB34.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/A819A39C-7016-F74A-94CC-482C4DB57D98.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/AC3E70C9-445C-6945-A6DE-759F429F452B.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/ACED4D7D-C2A9-2240-942D-878ED5622D89.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B129608F-C4D3-7C4E-9A29-E76713770924.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B1CC3074-9EE0-7742-AD70-81A2A0D261A6.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B24E3050-D9AB-2E4D-9CC0-EF17BE64CF03.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B3C75DE7-2DCB-9B46-AF7D-4BAC93DC2F48.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B48619A8-DC75-F649-9475-924E7ED58054.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B4B9F261-1DA2-1241-B3B6-05CE5A1BBE55.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B4C7E07D-FE3D-B44A-9F2C-F03BF5CF7FD4.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B4F1CD28-A6B8-5F48-A0DA-3B09FEE3AFF1.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B51390E7-5CCC-0348-B37F-A4FC3AC98106.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B57C7FBE-E51B-0948-B58E-06AB7B1F496C.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B5EBB61C-C221-B14C-95C6-DEECC85FE1E3.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/B9AE0A53-93D8-E547-9D63-BA02D547E003.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/BD1737FD-214F-B14E-8B71-5B1B5C6AD63F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/C0181147-6841-054B-B5FC-99AF6457F96F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/C19927F6-A6DE-2942-B11A-9DC4D9C4CE09.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/C4B43B24-47F9-CA47-83DE-8C32CA9758F6.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/C537EA44-703B-8B40-B4A0-0C2391614DBD.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/C61B2E9B-FDE1-CD4E-B390-D2EF8CFB1AB0.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/C79AC11A-D1DB-6E47-B49C-367DB3807C34.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/C866A534-8D3F-AB42-AB22-787E6A74FD9D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/CBF9E541-4071-8A40-9DD8-E3F9E56D8B21.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/CC9DC205-0C3A-CA4B-BB00-50CD5C7EC9BA.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/CEF04E1B-C335-0744-851A-D768B9BFEF75.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/D02DF83C-BD9C-5542-8D9F-B8CF2FF627C2.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/D432ACB7-515E-EF41-B3D1-BA3797D64484.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/D635D5F6-798F-C045-A6A7-DA0FDFF72708.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/D6C1D7F2-6930-6244-B66D-A5CAAAC16651.root",
    #
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/D6D9BACA-9BEB-014B-9B2F-EACC8695E1FD.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/D7A0F357-9CB3-0142-817F-45C80A9C67FD.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/D7F25FA3-05E1-A240-87BB-44CCE41707C3.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/DA193752-B3D8-FE4E-BD1D-55D26F74F1CC.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/DA302C5C-DF22-074F-9D24-33B745C70675.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/DAE4267D-BE43-A44B-A10B-1F277DEC7015.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/DCA0D214-B345-C34C-9CCB-C5E63BCE4E26.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/DCBBDD3E-922D-BA4C-BFF3-B268DAB2B0AA.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/DD253128-B0B9-864C-8938-E7A11C0BD9D1.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/E0F4A3F2-9828-A14C-90BE-2DAB7ABBAF18.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/E0F5D595-4DA9-C94E-A2AD-8CC4F41112D6.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/E1022087-07CE-864F-B4C6-A796523CEFD3.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/E5506605-BA91-4F41-94C8-95EB18FF60F2.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/E73A1C8E-EB43-CA4F-9A10-8D0DC8F2735F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/E75B9FB0-B81B-B746-8539-7D71ED1618BE.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/E804237E-3D38-4F4A-955A-43BE5C347262.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/EBCA9A50-E8E5-3F41-871B-9AE7DFBE1899.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/EEE10597-1BBC-2449-96AD-ADB7C6E542AE.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/F0109E6E-052B-FC4C-8F46-6B66C0DF71B1.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/F091C6E3-BFDB-1E40-AB74-859735D7D379.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/F18B495D-E863-7040-9576-3B59E497E1A2.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/F2A56CF5-73EA-B843-9562-D51CE269309E.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/F39D20DD-6A74-924E-9931-822D29DA605A.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/F4E7E036-353B-AD4F-BE25-7ECBCF58986F.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/F63B6E04-ED28-BF4D-9832-1774BF8D4C2B.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/F77A1AF4-A437-B746-BBC5-A65D3283746D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/FC639B31-C4EC-4249-8DF3-F24D204C93AD.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/FE16D5D5-06C0-454E-8AC2-FF6CDF84BE2C.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/FE61A0D8-CEDC-2142-81AA-06301F452792.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280000/FF7E23DA-5D9A-3C40-9391-8E8CDC88CCA2.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280001/541255AD-EA7F-1B42-9E03-DA7311B9078D.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280001/A3AEAC55-28C9-1845-BB5D-B4D67070DB59.root",
    #"/store/data/Run2018D/ZeroBias/AOD/12Nov2019_UL2018_rsb-v1/280001/DA5E8AD6-43A1-C14C-8AD3-92C04117E3EE.root"
    )  
))

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v26', '')

if options.calibInput == '':
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string('PPSTimingCalibrationRcd'),
                    tag = cms.string('PPSDiamondTimingCalibration_v1'),
                connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                )
    )
else:
    process.ppsTimingCalibrationESSource = cms.ESSource('PPSTimingCalibrationESSource',
        calibrationFile = cms.string(options.calibInput),
        subDetector = cms.uint32(2),
        appendToDataLabel = cms.string('')
    )

# rechits production
process.load(options.geometryFile)
process.load('RecoPPS.Local.ctppsDiamondLocalReconstruction_cff')

process.diamondTimingWorker = DQMEDAnalyzer("DiamondTimingWorker",
    tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
    tagRecHit = cms.InputTag("ctppsDiamondRecHits"),
    tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracks"),
    tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracks"),
    tagValidOOT = cms.int32(-1),
    planesConfig = cms.string(options.planesConfig),
    Ntracks_Lcuts = cms.vint32([-1,1,-1,1]),
    Ntracks_Ucuts = cms.vint32([-1,6,-1,6]),
)

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string(options.outputFile)
#)

process.DQMStore = cms.Service("DQMStore")
process.dqmOutput = cms.OutputModule("DQMRootOutputModule", fileName=cms.untracked.string(options.outputFileName))

#print(vars(process))

process.path = cms.Path(
   process.diamondTimingWorker
)

process.end_path = cms.EndPath(process.dqmOutput)

process.schedule = cms.Schedule(
    process.path,
    process.end_path
)
