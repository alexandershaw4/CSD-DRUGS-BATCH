

Config.Files  = '/cubric/scratch/waans/RAWVE_13_30/Datasets.txt'; % input files [array]
Config.MODELS = 'MODELS.m'; % models function 
Config.DESIGN = 'Design.m';

Config.mPath = '/cubric/scratch/waans/RAWVE_13_30/CSD-DRUGS-BATCH/';
Config.dPath = '/cubric/scratch/waans/RAWVE_13_30/CSD-DRUGS-BATCH/';


Config.TrialCode  = 1;
Config.pst        = [0 2000];
Config.FreqWin    = [13 30];
Config.Han        = 1;
Config.Downsample = 0;
Config.Filter     = [1 30];


Config.Anal    = 'CSD';
Config.Type    = 'LFP';
Config.Spatial = 'LFP';
Config.Neural  = 'cmm_nmda';

Config.Files = ReadverifyDatasets(Config.Files);

NEW_CSD_NMDA_DCM_AS(Config);
