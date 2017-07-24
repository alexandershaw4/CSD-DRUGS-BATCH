




Config.Files  = '/cubric/scratch/waans/RAWVE_13_30/Datasets.txt'; % input files [array]
Config.MODELS = '/cubric/scratch/waans/RAWVE_13_30/CSD-DRUGS-BATCH/MODELS.m'; % models function 
Config.DESIGN = '/cubric/scratch/waans/RAWVE_13_30/CSD-DRUGS-BATCH/Design.m';

Config.TrialCode  = 1;
Config.pst        = [0 1200];
Config.FreqWin    = [13 30];
Config.Han        = 1;
Config.Downsample = 0;
Config.Filter     = [1 30];

Config.Anal    = 'CSD';
Config.Type    = 'LFP';
Config.Spatial = 'LFP';
Config.Neural  = 'cmm_nmda';

try handles.TrialCode ; catch ; handles.TrialCode  = 1;        end
try handles.pst;        catch ; handles.pst        = [0 200];  end
try handles.FreqWin;    catch ; handles.FreqWin    = [1 100];  end
try handles.Han;        catch ; handles.Han        = 1;        end
try handles.Downsample; catch ; handles.Downsample = 0;        end
try handles.Filter;     catch ; handles.Filter     = [];       end

try handles.Anal;       catch ; handles.Anal       = 'CSD';    end
try handles.Type;       catch ; handles.Type       = 'LFP';    end
try handles.Spatial;    catch ; handles.Spatial    = 'LFP';    end
try handles.Neural;     catch ; handles.Neural     = 'cmm_nmda'; end
