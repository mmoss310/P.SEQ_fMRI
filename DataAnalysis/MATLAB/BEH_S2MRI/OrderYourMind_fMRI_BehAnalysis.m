function OrderYourMind_fMRI_BehAnalysis
%==========================================================================
% UPDATES:
% NOTES:
% 
%==========================================================================
% BASIC STEPS--------------------------------------------------------------
% 1, Join response related data and parameters into one dataset
% ->behavioralDS

% 2, Prepare additional parameters
% ->refinedDS,refinedDS_1,refinedDS_2 and etc

% 3, Apply cut off 
% ->processedDS

% 4, Do statistical analysis and Plotting
% (if needed, create a data spreadsheet for R analysis!!)


%%

%BEHAVIORAL DATA-----------------------------------------------------------
%Analysis Start!===========================================================
%==========================================================================



%%

clear all
close all

% Coding Basic Info
BASIC.MASTERDIR=uigetdir;
BASIC.EXCODE='ordymS1';
BASIC.FILETOLOAD_P='PARAM_ST';
BASIC.FILETOLOAD_R='OrderYourMind_fMRI_Beh';
BASIC.ALLGRANDIR=strcat(BASIC.MASTERDIR,'/w_ALLGRAND');
cd(BASIC.MASTERDIR);addpath(pwd);

% Get all subject directories...
cd(BASIC.ALLGRANDIR);
BASIC.PARAMFILES=dir('*PARAM_ST.mat');
subs=cellfun(@(f)str2double(f(regexp(f,'\d'))),{BASIC.PARAMFILES.name},'Uni',false);
BASIC.ALLSUBS=unique(vertcat(subs{:}));

%%


%Loop through all subjects and load parameter files!! 
allDS=cell(1,length(BASIC.ALLSUBS));
for sub=BASIC.ALLSUBS'
    %Load data parameter and record files
    s=BASIC.ALLSUBS==sub;
    
    %parameter files
    pr=v2struct(load(BASIC.PARAMFILES(s).name));
    btdim=size(pr.CHUNK);
    pr.SUBID=repmat(sub,btdim);
    pr.BLOCK=repmat([1:btdim(1)]',1,btdim(2));
    pr.TRIAL=repmat([1:btdim(2)],btdim(1),1);
    pr.PRACTICE=repmat(pr.PRACTICE,1,btdim(2));
    pr.SEQTYPE=repmat(pr.SEQTYPE,1,btdim(2));
    allDS(s)={strechSt2DS(pr)};
end

% Save PARAM DS FILE
cd(BASIC.ALLGRANDIR)
paramDS=vertcat(allDS{:});
saveMatFile(paramDS,'ALL_PARAM','ALL_PARAM');

% Create RESP DS FILE
txt2dataSet(strcat(BASIC.FILETOLOAD_R,'.txt'),'ALL_RESP');
respDS=v2struct(load('ALL_RESP.mat'));

setdiff(unique(respDS.SUBID),unique(paramDS.SUBID))
setdiff(unique(paramDS.SUBID),unique(respDS.SUBID))

% Merge them into one!!
% Base is respDS because subjects may not have compelte data!
mergeKeys={'SUBID','BLOCK','TRIAL'};
behDS=join(respDS,paramDS,'key',mergeKeys);
saveMatFile(behDS,'behDataS','behDataS');

%%

% % %Process2:Preparing all variables in dataset array
cd(BASIC.ALLGRANDIR)
data=v2struct(load('behDataS.mat'));
data=OrderYourMind_fMRI_ProcessDSVars('BEH',data);
% data=OrderYourMind_fMRI_ProcessDSVars('KGROUP',data);
saveMatFile(data,'refinedDataS','refinedDataS');

% % %Process2:Preparing all variables in dataset array
cd(BASIC.ALLGRANDIR)
data=v2struct(load('refinedDataS.mat'));
data=OrderYourMind_fMRI_ProcessDSVars('CUTOFF',data);
saveMatFile(data,'processedDataS','processedDataS');

% % Process3:Aggregate data and get averages of interest!
cd(BASIC.ALLGRANDIR)
data=v2struct(load('processedDataS.mat'));
export(data,'file','OrderYourMind_fMRI_BEHP.txt');

%%

%%
%PLOTTING!=================================================================
%==========================================================================


%%
%CORRELATION!==============================================================
%==========================================================================

%%


%%
% Done!
disp('Done!');

end