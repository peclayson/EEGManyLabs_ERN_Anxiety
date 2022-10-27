function scoring_orgpipe(wrkdir,savefile_ssa,savefile_sng)
%
%Function to score data using the direct replication pipeline for the
% #EEGManyLabs replication of the Hajcak et al. (2003) study of the
% relationship between ERN and anxiety.
%
%scoring_orgpipe(wrkdir,savefile_ssa,savefile_sng)
%
% This pilot script was developed using data from USF's psychophysiology
% lab. EEG was recorded using MagStim EGI hardware. Functionality for other
% formats will be built out as data are collected by other labs.
%
% The 'prepoc.m' script should have already been run on raw EEG data.
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
% MATLAB (R2021a) and EEGLab (v2022.1) are required
%

%History
% by Peter Clayson (10/20/22)
% peter.clayson@gmail.com
%
%
%
%
%

%check whether EEGLab is contained in the Matlab path
fprintf('\nEnsuring EEGLab is found in the Matlab path\n');

%Check for EEGLab
if exist('eeglab.m','file') ~= 2
    
    dlg = {'Warning: EEGlab is not found. EEGlab may not be installed';...
        'or EEGLab may not be located in the MATLAB path.';...
        'This script requires EEGLab to run.'};
    
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
    return;
else
    fprintf('EEGlab found\n');
end

%get names of all .set files in directory
rawfilesloc = dir(fullfile(wrkdir,'*.set'));
nSub = length(rawfilesloc);
files = struct2cell(rawfilesloc)';
files = files(:,1);

%specify where ERN should be scored
ERPinfo.mff.fz = 'E11';
ERPinfo.mff.cz = 'E129';
ERPinfo.mff.pz = 'E62';

%specify time window when ERN should be scored
ERPinfo.ern_wind = [0 150];

%specify event types
ERPinfo.events = {'cor' 'err'};

%loop through all subjects
for subj = 1:nSub
    
    %just keeping track of the participant being processed
    fprintf('Working on participant %d of %d\n', subj, length(files));
    
    %load EEG file
    EEG = pop_loadset('filename',files(subj),...
        'filepath',wrkdir);
    
    %score single-trial and subject averages
    [individ_ssa, individ_sng] = ern_score(files(subj), EEG, ERPinfo);
    
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


function [individ_ssa, individ_sng] = ern_score(fname, EEG, ERPinfo)

%remove extension for saving filename
savename = fname{1}(1:end-4);

%find location of fz, cz, and pz
chan_fz = find(strcmp(extractfield(EEG.chanlocs,'labels'),ERPinfo.mff.fz));
chan_cz = find(strcmp(extractfield(EEG.chanlocs,'labels'),ERPinfo.mff.cz));
chan_pz = find(strcmp(extractfield(EEG.chanlocs,'labels'),ERPinfo.mff.pz));

%cycle through all events
%data are already proessed using the 'prepoc.m' script and should only
%contain correct (cor) and error (err) trials.
for ii = 1:size(EEG.event,2)
    
    %make sure only a correct or error trial (just in case)
    if any(strcmpi(EEG.event(ii).type,ERPinfo.events))
        
        %score ERN at Fz, Cz, and Pz
        ern_fz = min(EEG.data(chan_fz,...
            knnsearch(EEG.times',ERPinfo.ern_wind(1)):...
            knnsearch(EEG.times',ERPinfo.ern_wind(2)),...
            ii));
        
        ern_cz = min(EEG.data(chan_cz,...
            knnsearch(EEG.times',ERPinfo.ern_wind(1)):...
            knnsearch(EEG.times',ERPinfo.ern_wind(2)),...
            ii));
        
        ern_pz = min(EEG.data(chan_pz,...
            knnsearch(EEG.times',ERPinfo.ern_wind(1)):...
            knnsearch(EEG.times',ERPinfo.ern_wind(2)),...
            ii));
        
        %store single-trial information
        if ~exist('individ_sng','var')
            
            individ_sng = table;
            individ_sng.subjid = cellstr(savename);
            individ_sng.event = cellstr(EEG.event(ii).type);
            individ_sng.ern_fz = ern_fz;
            individ_sng.ern_cz = ern_cz;
            individ_sng.ern_pz = ern_pz;
            
        else
            
            row = table;
            row.subjid = cellstr(savename);
            row.event = cellstr(EEG.event(ii).type);
            row.ern_fz = ern_fz;
            row.ern_cz = ern_cz;
            row.ern_pz = ern_pz;
            
            individ_sng = vertcat(individ_sng,row); %#ok<AGROW>
            
        end
    end
    
end

%pull event information and index correct or error trials
events_all = extractfield(EEG.event,'type');
events_cor = strcmp(events_all,'cor');
events_err = strcmp(events_all,'err');

%make subject averages for crn and ern
ssa_crn = mean(EEG.data(:,:,events_cor),3);
ssa_ern = mean(EEG.data(:,:,events_err),3);

%score data
crn_fz = min(ssa_crn(chan_fz,...
    knnsearch(EEG.times',ERPinfo.ern_wind(1)):...
    knnsearch(EEG.times',ERPinfo.ern_wind(2))));

crn_cz = min(ssa_crn(chan_cz,...
    knnsearch(EEG.times',ERPinfo.ern_wind(1)):...
    knnsearch(EEG.times',ERPinfo.ern_wind(2))));

crn_pz = min(ssa_crn(chan_pz,...
    knnsearch(EEG.times',ERPinfo.ern_wind(1)):...
    knnsearch(EEG.times',ERPinfo.ern_wind(2))));


ern_fz = min(ssa_ern(chan_fz,...
    knnsearch(EEG.times',ERPinfo.ern_wind(1)):...
    knnsearch(EEG.times',ERPinfo.ern_wind(2))));

ern_cz = min(ssa_ern(chan_cz,...
    knnsearch(EEG.times',ERPinfo.ern_wind(1)):...
    knnsearch(EEG.times',ERPinfo.ern_wind(2))));

ern_pz = min(ssa_ern(chan_pz,...
    knnsearch(EEG.times',ERPinfo.ern_wind(1)):...
    knnsearch(EEG.times',ERPinfo.ern_wind(2))));

%store information
individ_ssa = table;
individ_ssa.subjid = cellstr(savename);
individ_ssa.n_cor = sum(events_cor);
individ_ssa.n_err = sum(events_err);
individ_ssa.crn_fz = crn_fz;
individ_ssa.crn_cz = crn_cz;
individ_ssa.crn_pz = crn_pz;
individ_ssa.ern_fz = ern_fz;
individ_ssa.ern_cz = ern_cz;
individ_ssa.ern_pz = ern_pz;

end
