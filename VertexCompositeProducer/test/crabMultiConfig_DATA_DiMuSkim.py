from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ppRefSkim2016_DiMuContBoth_cfg.py'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/5TeV/ReReco/Cert_306546-306826_5TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt'
config.Data.runRange = '306546-306826'
config.Data.publication = True
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*','T2_US_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

dataMap = {
            "SingleMuon": { "PD": "/SingleMuon/Run2017G-17Nov2017-v1/AOD", "Units": 7, "Memory": 2400, "RunTime": 620 },
            "DoubleMuon": { "PD": "/DoubleMuon/Run2017G-17Nov2017-v1/AOD", "Units": 7, "Memory": 2500, "RunTime": 820 }
          }

## Submit the muon PDs
for key, val in dataMap.items():
    config.General.requestName = 'VertexCompositeSkim_'+key+'_Run2017G_DiMuMassMin2_20190416'
    config.Data.inputDataset = val["PD"]
    config.Data.unitsPerJob = val["Units"]
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/RiceHIN/ppRef2017/SKIM/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    #config.Data.outLFNDirBase = '/store/user/%s/RiceHIN/ppRef2017/SKIM/%s' % (getUsernameFromSiteDB(), config.General.requestName)
    submit(config)