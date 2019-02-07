if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    from CRABClient.UserUtilities import config, getUsernameFromSiteDB
    config = config()

    config.General.workArea = 'VertexCompositeAna_mtd'
    config.General.transferOutputs = True
    config.General.transferLogs = False
    config.JobType.pluginName = 'Analysis'
    config.JobType.maxMemoryMB = 3000
    config.JobType.maxJobRuntimeMin = 2750
    config.Data.unitsPerJob = 1
    config.Data.splitting = 'FileBased'
    config.Data.outLFNDirBase = '/eos/cms/store/group/phys_heavyions/yousen/promptd0'
    config.Data.publication = True
    #config.Data.inputDBS = 'phys03'
    config.Site.storageSite = 'T2_US_MIT'
#    config.Site.storageSite = 'T3_US_Rice'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = 'PromptD0_mc_mtd_v1'
    config.JobType.psetName = '../test/PromptD0_mc_mtd_cfg.py'
    config.Data.outputDatasetTag = 'promptd0_mc_mtd_Skim_v1'
    config.Data.inputDataset = '/D0_PiK_prompt_5p02TeV_TuneCP5_MTD/anstahll-D0_PiK_prompt_5p02TeV_TuneCP5_MTD_RECO_20190126-8c6abfe98cf25191ce912f4177fc40d8/USER'
    submit(config)
