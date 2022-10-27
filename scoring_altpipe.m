function scoring_altpipe(wrkdir,savefile_ssa,savefile_sng)

%
%Function to score data using the alternative replication pipeline for the
% #EEGManyLabs replication of the Hajcak et al. (2003) study of the
% relationship between ERN and anxiety.
%
%
%scoring_altpipe(wrkdir,savefile_ssa,savefile_sng)
%
% This pilot script was developed using data from USF's psychophysiology
% lab. EEG was recorded using MagStim EGI hardware. Functionality for other
% formats will be built out as data are collected by other labs.
%
% The 'prepoc.m' script should have already been run on raw EEG data. Then
% data should have been processed using the instructions through the ERP
% PCA Toolkit using the instructions in erp_pca_toolkit_processing.pdf.
%
%Required Inputs
% dirwork - directory where processed .set/.fdt files are located
% savefile_ssa - where to save .csv file with subject average scores
% savefile_sng - where to save .csv file with single trial scores
%
%Output
% No variables are outputted to the Matlab workspace
% Two .csv files will be save in the locations specified in savefile_ssa
%  and savefile_sng
%
%Required software and plugins (and tested versions)
% MATLAB (R2021a) is required
%

%History
% by Peter Clayson (10/20/22)
% peter.clayson@gmail.com
%
%
%
%
%

%get names of all .set files in directory
rawfilesloc = dir(fullfile(wrkdir,'*.ept'));
nSub = length(rawfilesloc);
files = struct2cell(rawfilesloc)';
files = files(:,1);

%specify where ERN should be scored
ERPinfo.mff.fz = 'E11';
ERPinfo.mff.cz = 'E129';
ERPinfo.mff.pz = 'E62';

%specify time window when ERN should be scored
ERPinfo.ern_wind = [0 100];

%specify event types
ERPinfo.events = {'cor' 'err'};

for subj = 1:nSub
    
    %load ep dataset
    EPdata = load(fullfile(wrkdir,char(files(subj))),'-mat');
    
    %score single-trial and subject averages
    [individ_ssa, individ_sng] = ern_score(files(subj), EPdata.EPdata, ERPinfo);
    
    %store information
    if ~exist('master_ssa','var')
        master_ssa = individ_ssa;
        master_sng = individ_sng;
    else
        master_ssa = vertcat(master_ssa,individ_ssa); %#ok<AGROW>
        master_sng = vertcat(master_sng,individ_sng); %#ok<AGROW>
    end
end


%save .csv files
writetable(master_ssa,savefile_ssa);
writetable(master_sng,savefile_sng);

end


function [individ_ssa, individ_sng] = ern_score(fname, EPdata, ERPinfo)

%remove extension for saving filename
savename = fname{1}(1:end-4);

%find location of fz, cz, and pz
chan_fz = find(strcmp(EPdata.chanNames,ERPinfo.mff.fz));
chan_cz = find(strcmp(EPdata.chanNames,ERPinfo.mff.cz));
chan_pz = find(strcmp(EPdata.chanNames,ERPinfo.mff.pz));

%cycle through trials
for ii = 1:length(EPdata.cellNames)
    
    %make sure the trial is "good"
    if EPdata.analysis.badTrials(ii) == 0
        
        %make sure only a correct or error trial (just in case)
        if any(strcmpi(EPdata.cellNames(ii),ERPinfo.events))
            
            %score ERN at Fz, Cz, and Pz
            ern_fz = mean(EPdata.data(chan_fz,...
                knnsearch(EPdata.timeNames,ERPinfo.ern_wind(1)):...
                knnsearch(EPdata.timeNames,ERPinfo.ern_wind(2)),...
                ii));
            
            ern_cz = mean(EPdata.data(chan_cz,...
                knnsearch(EPdata.timeNames,ERPinfo.ern_wind(1)):...
                knnsearch(EPdata.timeNames,ERPinfo.ern_wind(2)),...
                ii));
            
            ern_pz = mean(EPdata.data(chan_pz,...
                knnsearch(EPdata.timeNames,ERPinfo.ern_wind(1)):...
                knnsearch(EPdata.timeNames,ERPinfo.ern_wind(2)),...
                ii));
            
            %store single-trial information
            if ~exist('individ_sng','var')
                
                individ_sng = table;
                individ_sng.subjid = cellstr(savename);
                individ_sng.event = cellstr(EPdata.cellNames(ii));
                individ_sng.ern_fz = ern_fz;
                individ_sng.ern_cz = ern_cz;
                individ_sng.ern_pz = ern_pz;
                
            else
                
                row = table;
                row.subjid = cellstr(savename);
                row.event = cellstr(EPdata.cellNames(ii));
                row.ern_fz = ern_fz;
                row.ern_cz = ern_cz;
                row.ern_pz = ern_pz;
                
                individ_sng = vertcat(individ_sng,row); %#ok<AGROW>
                
            end
        end
    end
end

%because a time-window mean amplitude was used, the average of single-trial
%scores is equal to the time-window mean amplitde of the average waveform
%therefore, no need to separately average and score data

%pull event information and index correct or error trials
events_cor = strcmp(EPdata.cellNames,'cor');
events_err = strcmp(EPdata.cellNames,'err');



%store information for subject averaged data
individ_ssa = table;
individ_ssa.subjid = cellstr(savename);
individ_ssa.n_cor = sum(events_cor);
individ_ssa.n_err = sum(events_err);
individ_ssa.crn_fz = mean(individ_sng.ern_fz(strcmp(individ_sng.event,'cor')));
individ_ssa.crn_cz = mean(individ_sng.ern_cz(strcmp(individ_sng.event,'cor')));
individ_ssa.crn_pz = mean(individ_sng.ern_pz(strcmp(individ_sng.event,'cor')));
individ_ssa.ern_fz = mean(individ_sng.ern_fz(strcmp(individ_sng.event,'err')));
individ_ssa.ern_cz = mean(individ_sng.ern_cz(strcmp(individ_sng.event,'err')));
individ_ssa.ern_pz = mean(individ_sng.ern_pz(strcmp(individ_sng.event,'err')));

end
