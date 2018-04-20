from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'PowhegHerwigPP_PL'
config.General.workArea = 'test'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab.py', 'iheartny_topxs_fwlite.py', 'JECs', 'pileup_reweight_mu.root']
config.JobType.outputFiles = ['test_iheartNY.root']
config.JobType.scriptExe = 'execute_iheartNY_ttbar_fullTruth.sh'

config.section_("Data")
config.Data.inputDataset = ' /TT_TuneEE5C_13TeV-powheg-herwigpp/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v1_v2p4-bfea8033ab2d179bbb8e0faf6e2dc0cf/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 1000
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/skinnari/13TeV_full2016'
# This string is used to construct the output dataset name

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.whitelist = ["T2_HU_Budapest"]
