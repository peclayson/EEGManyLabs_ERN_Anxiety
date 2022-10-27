function preproc(dirwork, dirsave, dirsave_alt, useParallel)
%
%Function to import data for the #EEGManyLabs replication of the Hajcak et
% al. (2003) study of the relationship between ERN and anxiety.
%
%preproc(dirwork, dirsave, dirsave_alt, useParallel)
%
% This pilot script was developed using data from USF's psychophysiology
% lab. EEG was recorded using MagStim EGI hardware. Functionality for other
% formats will be built out as data are collected by other labs.
%
% Location: Hardware (extension)
% USF: MagStim EGI (.mff)
%Required Inputs
% dirwork - directory where raw files are located
% dirsave - directory where files should be saved
% dirsave_alt - directory where files for EP Toolkit should be saved
%
%Optional Input
% useParallel - 1 = run files in parallel; 0 = run files serially (default)
%
%Output
% No variables are outputted to the Matlab workspace
% An eeglab dataset will be saved in the dirsave directory
%
%Required software and plugins (and tested versions)
% MATLAB (R2021a), EEGLab (v2022.1), and ERPLab (V9.00) are required for
%  data importing and processing. Additionally, the following plugins are
%  required for importing various formats.
%    MFFMatlabIO (v4.0) for .mff files
% ERP PCA toolkit (v2.97) will be used in the alternate data-processing 
%  pipeline.
%
%
% gratton_emcp.m from
%  https://github.com/kylemath/MathewsonMatlabTools/blob/master/EEG_analysis/gratton_emcp.m
%  The script will be used for ocular artifact correction. The script
%  contains detailed comments and is included in the GitHub repo for the
%  current project. It is used with permission from the author, Bill Gehring.

%History
% by Peter Clayson (10/6/22)
% peter.clayson@gmail.com
%
%
%
%
%


%check whether EEGLab, ERPLab, and EP Toolkit are contained in the Matlab path
fprintf('\nEnsuring EEGLab, ERPLab, and the ERP PCA Toolkit are found in the Matlab path\n');

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

%Check for ERPLab
if exist('eegplugin_erplab.m','file') ~= 2
    
    dlg = {'Warning: ERPLab is not found. ERPLab may not be installed';...
        'or ERPLab may not be located in the MATLAB path.';...
        'This script requires ERPLab to run.'};
    
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
    return;
else
    fprintf('ERPlab found\n');
end

%Check for ERP PCA Toolkit
if exist('ep.m','file') ~= 2
    
    dlg = {'Warning: EP Toolkit is not found. EP Toolkit may not be installed';...
        'or EP Toolkit may not be located in the MATLAB path.';...
        'This script requires the ERP PCA Toolkit to run.'};
    
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
    return;
else
    fprintf('ERP PCA Toolkit found\n');
end


%If user does not indicate whether parallel processing should be used,
%assume that the user does NOT want to process data using all processors in
%parallel
if nargin < 4
    useParallel = 0;
end

%If user indicates that files should be processed in parallel, ensure that
%the parallel processing toolbox is intstalled. If the toolbox is not
%installed, fun files serially.
if useParallel == 1
    if exist('parpool.m','file') ~= 2
        
        dlg = {'Warning: Parallel toolbox not installed';...
            'Data will be run serially rather than in parallel.'};
        
        for ii = 1:length(dlg)
            fprintf('%s\n',dlg{ii});
        end
        useParallel = 0;
        fprintf('\n\n');
    else
        fprintf('Parallel Toolbox installed\n');
    end
end



%Set a variable to determine whether any files are actually found. It will
%remain 0 until appropriate files are located.
filesfound = 0;

%Check for file in dirwork and filter based on extensions. If an extension
%is found ensure that the appropriate plugin, if necessary, is installed.

%Pull .mff file names from dirwork.
mfflist=dir(fullfile(dirwork,'*.mff'));
if ~isempty(mfflist)
    if filesfound == 1
        filenames=[filenames,{mfflist.name}];
    elseif filesfound == 0
        filenames={mfflist.name};
        filesfound = 1;
    end
    
    %Check for MFFMatlabIO plugin
    if exist('eegplugin_mffmatlabio.m','file') ~= 2
        
        plg_dlg = {'Warning: MFFMatlabIO plugin is not found.';...
            'This script requires the MFFMatlabIO plugin to load .mff files.'};
        
        for ii = 1:length(plg_dlg)
            fprintf('%s\n',plg_dlg{ii});
        end
        fprintf('\n\n');
        return;
    else
        fprintf('MFFMatlabIO plugin found\n');
    end
    
else
    dlg = {'Warning: There were no .mff files found in the specified directory.';...
        'Please specify the correct directory with .mff files to process.';...
        'Ignore this message if no .mff files should be processed.'};
    
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
end





%Print the date and time that processing began
fprintf(...
    '\n******\nProcessing began on %s at %s \n******\n\n',...
    date, datestr(now, 'HH:MM:SS'));


if filesfound == 0
    dlg = {'Warning: There were no raw EEG files found in the';...
        'specified directory (dirwork).';...
        'Please specify the correct directory with EEG files to process.'};
    for ii = 1:length(dlg)
        fprintf('%s\n',dlg{ii});
    end
    fprintf('\n\n');
    return;
end


%start a timer
t1 = tic;

%if the user would like to run data in parallel
if useParallel == 1
    
    %use all but one physical core for data processing
    nCores = feature('numCores') - 1;
    
    %create a parallel pool
    parpool(nCores);
    
    %Loop through all subjects in parallel
    parfor ii = 1:length(filenames)
        
        %indicate which participant is being processed and the progress
        fprintf('\n******\nProcessing participant %s\n******\n',filenames{ii});
        
        %process the data file
        prepare_eeg(filenames{ii}, dirwork, dirsave, dirsave_alt);
        
    end
    
    %close the pool of workers
    delete(gcp);
    
    %otherwise run data serially
elseif useParallel == 0
    
    %Loop through all subjects in parallel
    for ii = 1:length(filenames)
        
        %indicate which participant is being processed and the progress
        fprintf('\n******\nProcessing participant %s\n******\n',filenames{ii});
        
        %process the data file
        prepare_eeg(filenames{ii}, dirwork, dirsave, dirsave_alt);
        
        %just keeping track of the participant being processed
        fprintf('Finished participant %d of %d\n', ii, length(filenames));
        
    end
    
end


%stop timer
t2 = toc(t1);

%incdicate how long it took to process the files
tMinutes = round(t2/60);
disp(['It took ' num2str(tMinutes) ' minutes to process these data']);

end

function prepare_eeg(subject, dirwork, dirsave, dirsave_alt)

%parse the name of the file
[~,savename,eeg_ext] = fileparts(subject);

%Load eeg file, then rereference to link mastoids, and if needed load
% channel locations
switch eeg_ext
    case '.mff'
        %load .mff file
        EEG = pop_mffimport(fullfile(dirwork,subject), {'code'});
        
        %locate erp pca toolkit to make a path to the location of the channel
        % location file
        eplabdir = fileparts(which('ep.m'));
        
        %attach channel locations
        EEG = pop_chanedit(EEG, 'lookup',...
            fullfile(eplabdir,...
            'electrodes',...
            'old EGI Hydrocel',...
            'GSN-Hydrocel-129.ced'));
        
        %relevant channels for analysis later
        %E129 (Cz) is not included here, because it is the reference channel
        %It will be interpolated later
        %E11 = Fz
        %E62 = Pz
        imp_chans = {'E11','E62'};
           
end

%identify relevant channel indices for later artifact rejection
chan_ind = [];

for jj = 1:length(imp_chans)
    chan_ind(end+1) = find(strcmp({EEG.chanlocs.labels},...
        imp_chans{jj}));
end


%Resample to 200Hz
EEG = pop_resample(EEG, 200);



%The following code will re-mark the file based on specified bins.
%The code will overwrite the EEG.event structure and leave the
%EEG.urevent structure intact.
switch eeg_ext
    case '.mff'
        
        %the mff format is a bit troublesome
        %information is stored in the EEG.event.code = TRSP event
        %that information pertains to the prevous trial
        %the code below searches for responses, looks at the TRSP,
        %determines whether the response was a correct or error trial, and
        %whether it was a response within 200-800ms
        
        %mffkey_eval = 1 Correct
        %mffkey_eval = 0 Error
        %mffkey_rtim = Response times
        
        %get the number of events to cycle through
        [~,numb] = size(EEG.event);
        
        %serially go through all trials to recode based on response times
        for ii = 1:numb
            
            %find participant response
            if strcmp(EEG.event(ii).code,'resp')
                
                %Make sure the event info was recorded, otherwise
                %script will crash if task crashed (and the last trial
                %event information was not recorded in e-prime)
                if (ii+1 <= numb)
                    
                    %correct trial (mffkey_eval = 1)
                    if strcmp(EEG.event(ii+1).mffkey_eval,'1')
                        
                        %check RT is >200 and <800ms
                        if (str2double(EEG.event(ii+1).mffkey_rtim) > 200 &&...
                                str2double(EEG.event(ii+1).mffkey_rtim) < 800)
                            
                            %Code the trial type as correct
                            EEG.event(ii).type = '1'; %correct
                            
                        end
                        
                        %error trial (mffkey_eval = 0)
                    elseif strcmp(EEG.event(ii+1).mffkey_eval,'0')
                        
                        %check RT is >200 and <800ms
                        if (str2double(EEG.event(ii+1).mffkey_rtim) > 200 &&...
                                str2double(EEG.event(ii+1).mffkey_rtim) < 800)
                            
                            %Code the trial type as error
                            EEG.event(ii).type = '0'; %error
                            
                        end
                        
                    end
                end
            end %if strcmp(EEG.event(ii).code,'resp')
        end %for ii = 1:numb
        
        extra_fields = {'description' 'label'...
            'mffkey_age_' 'mffkey_exp_' 'mffkey_hand'...
            'mffkey_sex_' 'mffkey_subj' 'mffkey_cel' 'mffkey_obs'...
            'mffkey_spc' 'mffkey_pos' 'mffkey_argu' 'mffkey_rsp'...
            'mffkey_eval' 'mffkey_rtim' 'mffkey_trl' 'mffkey_Corr'...
            'mffkey_Accu'};
        
        EEG.event = rmfield(EEG.event,extra_fields);
        
        
end


%Data will go through two processing pipelines. One will mirror the
%original study, and the other will follow an alternative pipeline. This
%alternative dataset will be used later in the script.
EEG_alt = EEG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%Original Data Processing Pipepline%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%per Hajcak et al. (2003)%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Use ERPLab for filtering
%IIR Butterworth
%4th order (24 dB/oct)
%.05 to 35 Hz half-amplitude cutoffs
EEG = pop_basicfilter(EEG, 1:EEG.nbchan,...
    'Cutoff', [0.05 35],...
    'Design', 'butter',...
    'Filter', 'bandpass',...
    'Order', 4,...
    'Boundary','');

%epoch -100 to 800ms
EEG = pop_epoch(EEG, {'0' '1'},...
    [-.1 .8], 'newname', 'epoched', 'epochinfo', 'yes');

%remove unused epoch markers
EEG = pop_selectevent(EEG,... 
    'type', [0 1],...
    'deleteevents','on',...
    'deleteepochs','on',...
    'invertepochs','off');

%baseline adjust before emcp
EEG = pop_rmbase(EEG, [-100 0]);

%emcp using Gratton/Miller approach for regression-based correction
switch eeg_ext
    case '.mff'
        EEG = gratton_emcp(EEG,...
            {'1', '0'},... %events to consider
            {'E25 E8','E127 E126'},... %vertical EOG channels
            {'E128','E125'}); %horizontal EOG channels
end

%baseline adjust after emcp
EEG = pop_rmbase(EEG, [-100 0]);


%reject epochs according to these criteria...
% out of A/D range
% flat signal longer than 25 ms (i.e., 5 samples since data are
%  downsampled to 200Hz)

% only channel indices in chan_ind will be used for artifact rejection

%start with all good trials (good = 1, bad = 0)
ntrials = size(EEG.epoch,2);
good_ind = zeros(1,ntrials);

for trl = 1:ntrials
    
    %good trial unless found bad
    good_trl = 1;
    
    for chn = 1:length(chan_ind)
        
        %First: out of A/D range (any NaN in data?)
        %determine whether data contain any NaN
        if any(isnan(EEG.data(chan_ind(chn),:,trl)))
            good_trl = 0;
        end
        
        %Second: flat signal longer than 25 ms (5 samples)
        consec_diffs = diff(EEG.data(chan_ind(chn),:,trl),5);
        if any(consec_diffs == 0)
            good_trl = 0;
        end
        
    end
    
    good_ind(trl) = good_trl;
    
end

%remove any bad epochs based on above criteria
if any(good_ind == 0)
    
    %remove bad epochs
    EEG = pop_selectevent(EEG,...
        'epoch', find(good_ind == 0),...
        'deleteevents','off',...
        'deleteepochs','on',...
        'invertepochs','on');
    
end



%For EGI data, channel 129 (Cz) needs to be interpolated due to its
% use as the online reference channel
%All data need to be re-referenced to algabraically linked near mastoid
% channels
switch eeg_ext
    case '.mff'
        EEG = pop_interp(EEG, 129, 'spherical');
        
             %re-reference to linked mastoids
        EEG = pop_eegchanoperator(EEG, {...
            'nch1 = ch1 - ((ch57 + ch100)/2) Label E1',...
            'nch2 = ch2 - ((ch57 + ch100)/2) Label E2',...
            'nch3 = ch3 - ((ch57 + ch100)/2) Label E3',...
            'nch4 = ch4 - ((ch57 + ch100)/2) Label E4',...
            'nch5 = ch5 - ((ch57 + ch100)/2) Label E5',...
            'nch6 = ch6 - ((ch57 + ch100)/2) Label E6',...
            'nch7 = ch7 - ((ch57 + ch100)/2) Label E7',...
            'nch8 = ch8 - ((ch57 + ch100)/2) Label E8',...
            'nch9 = ch9 - ((ch57 + ch100)/2) Label E9',...
            'nch10 = ch10 - ((ch57 + ch100)/2) Label E10',...
            'nch11 = ch11 - ((ch57 + ch100)/2) Label E11',...
            'nch12 = ch12 - ((ch57 + ch100)/2) Label E12',...
            'nch13 = ch13 - ((ch57 + ch100)/2) Label E13',...
            'nch14 = ch14 - ((ch57 + ch100)/2) Label E14',...
            'nch15 = ch15 - ((ch57 + ch100)/2) Label E15',...
            'nch16 = ch16 - ((ch57 + ch100)/2) Label E16',...
            'nch17 = ch17 - ((ch57 + ch100)/2) Label E17',...
            'nch18 = ch18 - ((ch57 + ch100)/2) Label E18',...
            'nch19 = ch19 - ((ch57 + ch100)/2) Label E19',...
            'nch20 = ch20 - ((ch57 + ch100)/2) Label E20',...
            'nch21 = ch21 - ((ch57 + ch100)/2) Label E21',...
            'nch22 = ch22 - ((ch57 + ch100)/2) Label E22',...
            'nch23 = ch23 - ((ch57 + ch100)/2) Label E23',...
            'nch24 = ch24 - ((ch57 + ch100)/2) Label E24',...
            'nch25 = ch25 - ((ch57 + ch100)/2) Label E25',...
            'nch26 = ch26 - ((ch57 + ch100)/2) Label E26',...
            'nch27 = ch27 - ((ch57 + ch100)/2) Label E27',...
            'nch28 = ch28 - ((ch57 + ch100)/2) Label E28',...
            'nch29 = ch29 - ((ch57 + ch100)/2) Label E29',...
            'nch30 = ch30 - ((ch57 + ch100)/2) Label E30',...
            'nch31 = ch31 - ((ch57 + ch100)/2) Label E31',...
            'nch32 = ch32 - ((ch57 + ch100)/2) Label E32',...
            'nch33 = ch33 - ((ch57 + ch100)/2) Label E33',...
            'nch34 = ch34 - ((ch57 + ch100)/2) Label E34',...
            'nch35 = ch35 - ((ch57 + ch100)/2) Label E35',...
            'nch36 = ch36 - ((ch57 + ch100)/2) Label E36',...
            'nch37 = ch37 - ((ch57 + ch100)/2) Label E37',...
            'nch38 = ch38 - ((ch57 + ch100)/2) Label E38',...
            'nch39 = ch39 - ((ch57 + ch100)/2) Label E39',...
            'nch40 = ch40 - ((ch57 + ch100)/2) Label E40',...
            'nch41 = ch41 - ((ch57 + ch100)/2) Label E41',...
            'nch42 = ch42 - ((ch57 + ch100)/2) Label E42',...
            'nch43 = ch43 - ((ch57 + ch100)/2) Label E43',...
            'nch44 = ch44 - ((ch57 + ch100)/2) Label E44',...
            'nch45 = ch45 - ((ch57 + ch100)/2) Label E45',...
            'nch46 = ch46 - ((ch57 + ch100)/2) Label E46',...
            'nch47 = ch47 - ((ch57 + ch100)/2) Label E47',...
            'nch48 = ch48 - ((ch57 + ch100)/2) Label E48',...
            'nch49 = ch49 - ((ch57 + ch100)/2) Label E49',...
            'nch50 = ch50 - ((ch57 + ch100)/2) Label E50',...
            'nch51 = ch51 - ((ch57 + ch100)/2) Label E51',...
            'nch52 = ch52 - ((ch57 + ch100)/2) Label E52',...
            'nch53 = ch53 - ((ch57 + ch100)/2) Label E53',...
            'nch54 = ch54 - ((ch57 + ch100)/2) Label E54',...
            'nch55 = ch55 - ((ch57 + ch100)/2) Label E55',...
            'nch56 = ch56 - ((ch57 + ch100)/2) Label E56',...
            'nch57 = ch57 - ((ch57 + ch100)/2) Label E57',...
            'nch58 = ch58 - ((ch57 + ch100)/2) Label E58',...
            'nch59 = ch59 - ((ch57 + ch100)/2) Label E59',...
            'nch60 = ch60 - ((ch57 + ch100)/2) Label E60',...
            'nch61 = ch61 - ((ch57 + ch100)/2) Label E61',...
            'nch62 = ch62 - ((ch57 + ch100)/2) Label E62',...
            'nch63 = ch63 - ((ch57 + ch100)/2) Label E63',...
            'nch64 = ch64 - ((ch57 + ch100)/2) Label E64',...
            'nch65 = ch65 - ((ch57 + ch100)/2) Label E65',...
            'nch66 = ch66 - ((ch57 + ch100)/2) Label E66',...
            'nch67 = ch67 - ((ch57 + ch100)/2) Label E67',...
            'nch68 = ch68 - ((ch57 + ch100)/2) Label E68',...
            'nch69 = ch69 - ((ch57 + ch100)/2) Label E69',...
            'nch70 = ch70 - ((ch57 + ch100)/2) Label E70',...
            'nch71 = ch71 - ((ch57 + ch100)/2) Label E71',...
            'nch72 = ch72 - ((ch57 + ch100)/2) Label E72',...
            'nch73 = ch73 - ((ch57 + ch100)/2) Label E73',...
            'nch74 = ch74 - ((ch57 + ch100)/2) Label E74',...
            'nch75 = ch75 - ((ch57 + ch100)/2) Label E75',...
            'nch76 = ch76 - ((ch57 + ch100)/2) Label E76',...
            'nch77 = ch77 - ((ch57 + ch100)/2) Label E77',...
            'nch78 = ch78 - ((ch57 + ch100)/2) Label E78',...
            'nch79 = ch79 - ((ch57 + ch100)/2) Label E79',...
            'nch80 = ch80 - ((ch57 + ch100)/2) Label E80',...
            'nch81 = ch81 - ((ch57 + ch100)/2) Label E81',...
            'nch82 = ch82 - ((ch57 + ch100)/2) Label E82',...
            'nch83 = ch83 - ((ch57 + ch100)/2) Label E83',...
            'nch84 = ch84 - ((ch57 + ch100)/2) Label E84',...
            'nch85 = ch85 - ((ch57 + ch100)/2) Label E85',...
            'nch86 = ch86 - ((ch57 + ch100)/2) Label E86',...
            'nch87 = ch87 - ((ch57 + ch100)/2) Label E87',...
            'nch88 = ch88 - ((ch57 + ch100)/2) Label E88',...
            'nch89 = ch89 - ((ch57 + ch100)/2) Label E89',...
            'nch90 = ch90 - ((ch57 + ch100)/2) Label E90',...
            'nch91 = ch91 - ((ch57 + ch100)/2) Label E91',...
            'nch92 = ch92 - ((ch57 + ch100)/2) Label E92',...
            'nch93 = ch93 - ((ch57 + ch100)/2) Label E93',...
            'nch94 = ch94 - ((ch57 + ch100)/2) Label E94',...
            'nch95 = ch95 - ((ch57 + ch100)/2) Label E95',...
            'nch96 = ch96 - ((ch57 + ch100)/2) Label E96',...
            'nch97 = ch97 - ((ch57 + ch100)/2) Label E97',...
            'nch98 = ch98 - ((ch57 + ch100)/2) Label E98',...
            'nch99 = ch99 - ((ch57 + ch100)/2) Label E99',...
            'nch100 = ch100 - ((ch57 + ch100)/2) Label E100',...
            'nch101 = ch101 - ((ch57 + ch100)/2) Label E101',...
            'nch102 = ch102 - ((ch57 + ch100)/2) Label E102',...
            'nch103 = ch103 - ((ch57 + ch100)/2) Label E103',...
            'nch104 = ch104 - ((ch57 + ch100)/2) Label E104',...
            'nch105 = ch105 - ((ch57 + ch100)/2) Label E105',...
            'nch106 = ch106 - ((ch57 + ch100)/2) Label E106',...
            'nch107 = ch107 - ((ch57 + ch100)/2) Label E107',...
            'nch108 = ch108 - ((ch57 + ch100)/2) Label E108',...
            'nch109 = ch109 - ((ch57 + ch100)/2) Label E109',...
            'nch110 = ch110 - ((ch57 + ch100)/2) Label E110',...
            'nch111 = ch111 - ((ch57 + ch100)/2) Label E111',...
            'nch112 = ch112 - ((ch57 + ch100)/2) Label E112',...
            'nch113 = ch113 - ((ch57 + ch100)/2) Label E113',...
            'nch114 = ch114 - ((ch57 + ch100)/2) Label E114',...
            'nch115 = ch115 - ((ch57 + ch100)/2) Label E115',...
            'nch116 = ch116 - ((ch57 + ch100)/2) Label E116',...
            'nch117 = ch117 - ((ch57 + ch100)/2) Label E117',...
            'nch118 = ch118 - ((ch57 + ch100)/2) Label E118',...
            'nch119 = ch119 - ((ch57 + ch100)/2) Label E119',...
            'nch120 = ch120 - ((ch57 + ch100)/2) Label E120',...
            'nch121 = ch121 - ((ch57 + ch100)/2) Label E121',...
            'nch122 = ch122 - ((ch57 + ch100)/2) Label E122',...
            'nch123 = ch123 - ((ch57 + ch100)/2) Label E123',...
            'nch124 = ch124 - ((ch57 + ch100)/2) Label E124',...
            'nch125 = ch125 - ((ch57 + ch100)/2) Label E125',...
            'nch126 = ch126 - ((ch57 + ch100)/2) Label E126',...
            'nch127 = ch127 - ((ch57 + ch100)/2) Label E127',...
            'nch128 = ch128 - ((ch57 + ch100)/2) Label E128',...
            'nch129 = ch129 - ((ch57 + ch100)/2) Label E129'...
            },...
            'ErrorMsg', 'popup',...
            'Warning', 'on',...
            'Saveas', 'off');
end


%go through and remark trials with a label that is easier to use. It
%makes life easier down the road
[~, ltrials]=size(EEG.event);

for i = 1:ltrials
    switch EEG.event(i).type
        case '1'
            EEG.event(i).type = 'cor';
        case '0'
            EEG.event(i).type = 'err';
    end
end

%Save as .set file
pop_saveset(EEG, 'filename',[savename '.set'],...
    'filepath', dirsave);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%Alternative Data Processing Pipepline%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use ERPLab for filtering
%IIR Butterworth
%4th order (24 dB/oct)
%.01 to 30 Hz half-amplitude cutoffs
EEG_alt = pop_basicfilter(EEG_alt, 1:EEG_alt.nbchan,...
    'Cutoff', [0.01 30],...
    'Design', 'butter',...
    'Filter', 'bandpass',...
    'Order', 4,...
    'Boundary','');

%eopch -400 to 800 ms
EEG_alt = pop_epoch(EEG_alt, {'0' '1'},...
    [-.4 .8], 'newname', 'epoched', 'epochinfo', 'yes');

%remove unused epoch markers
EEG_alt = pop_selectevent(EEG_alt,... 
    'type', [0 1],...
    'deleteevents','on',...
    'deleteepochs','on',...
    'invertepochs','off');

%go through and remark trials with a label that is easier to use. It
%makes life easier down the road
[~, ltrials]=size(EEG_alt.event);

for i = 1:ltrials
    switch EEG_alt.event(i).type
        case '1'
            EEG_alt.event(i).type = 'cor';
        case '0'
            EEG_alt.event(i).type = 'err';
    end
end

%Save as .set file; all other steps will be performed using the EP Toolkit
pop_saveset(EEG_alt, 'filename',[savename '_alt.set'],...
    'filepath', dirsave_alt);

end
