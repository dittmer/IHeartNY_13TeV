from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'PowhegPythia8_herwigpp2_fullTruth'
config.General.workArea = 'test'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab.py', 'iheartny_topxs_fwlite.py', 'JECs', 'pileup_reweight_mu.root']
config.JobType.outputFiles = ['test_iheartNY.root']
config.JobType.scriptExe = 'execute_iheartNY_ttbar_fullTruth.sh'

config.section_("Data")
config.Data.inputDataset = '/TT_TuneEE5C_13TeV-powheg-herwigpp/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_ext2-v1_v2p4-bfea8033ab2d179bbb8e0faf6e2dc0cf/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/eyandel'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
