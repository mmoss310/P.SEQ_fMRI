function OrderYourMindfMRI_beh
%==========================================================================
% PURPOSE:
% TO DO:
% 1, adjut feedback time for FMRI
% 2, check fMRI blocks

% SETTING:=================================================================
% Session 1 (behavior)
% SEQUENCE = 12 sequences made of 2 chunks with 3 elements
% RETRIEVAL = 1000 ms
% CUTOFF = 750 ms
% PACING = fixed

% Session 2 (fMRI)
% SEQUENCE = 12 sequences made of 2 chunks with 3 elements
% RETRIEVAL = 1000 ms
% CUTOFF = 1000 ms
% PACING = fixed
%==========================================================================
clear all;
close all;
sca;
if 0, Screen('Preference', 'SkipSyncTests', 1); end
Screen('Preference', 'SkipSyncTests', 1);
commandwindow;while KbCheck; end
%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST PREX
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL RUNID BTRIAL BLOCKTS TRIGERTS
%**************************************************************************

%GET SUBJECT ID & OPEN FILES ==============================================
%FILE PATH SETTING
FPNAME_ST.TXTCODE='ordch_beh';
FPNAME_ST.EXFILE='OrderYourMindfMRI_beh.m';
FPNAME_ST.PTLPATH=PsychtoolboxRoot;
exFilePath=which(FPNAME_ST.EXFILE);
FPNAME_ST.EXFPATH=exFilePath(1:end-length(FPNAME_ST.EXFILE));
cd(FPNAME_ST.EXFPATH);
% FPNAME_ST.PORTID_S={'DCC8'};
% FPNAME_ST.PORTID_R={'DCD8'};
% FPNAME_ST.PORDIR={strcat(pwd,'\PortCode_Supports')};

%SEQUENCE SETTING
SEQ_ST.BASE(1,:)=[1,2,3];
SEQ_ST.BASE(2,:)=[1,3,2];
SEQ_ST.BASE(3,:)=[2,1,3];
SEQ_ST.BASE(4,:)=[2,3,1];
SEQ_ST.BASE(5,:)=[3,1,2];
SEQ_ST.BASE(6,:)=[3,2,1];
SEQ_ST.NUMCHUNK=size(SEQ_ST.BASE,1);
%Ignoring chunk associations (12 sequences)
chunkpairs(1,:)=[1,2];chunkpairs(2,:)=[1,3];
chunkpairs(3,:)=[2,1];chunkpairs(4,:)=[2,5];
chunkpairs(5,:)=[3,4];chunkpairs(6,:)=[3,1];
chunkpairs(7,:)=[4,3];chunkpairs(8,:)=[4,6];
chunkpairs(9,:)=[5,6];chunkpairs(10,:)=[5,2];
chunkpairs(11,:)=[6,5];chunkpairs(12,:)=[6,4];
SEQ_ST.CHUNKPAIRS=chunkpairs;
SEQ_ST.NUMSEQTYPE=size(chunkpairs,1);
SEQ_ST.SERTYPE=[1:SEQ_ST.NUMSEQTYPE];%sequence type
SEQ_ST.MAXSEQLENGTH=[6];%smaximum elements in sequence
seqpairs(1,:)=[1,7,9];seqpairs(2,:)=[2,8,10];
seqpairs(3,:)=[3,5,11];seqpairs(4,:)=[4,6,12];
SEQ_ST.SEQPAIRS=seqpairs;
SEQ_ST.ONTRACKRATE=0.50;

%EXPERIMENT
% BEH: 12 runs = 4 runs (3 seq types) * 3
% fMRI: 8 runs = 4 runs (3 seq types) * 2
EXC_ST.NUMPHASE(1)=3;% one phase contains 12 sequences
EXC_ST.NUMPHASE(2)=2;
EXC_ST.NUMEBLOCK(1)=SEQ_ST.NUMSEQTYPE*EXC_ST.NUMPHASE(1);
EXC_ST.NUMEBLOCK(2)=SEQ_ST.NUMSEQTYPE*EXC_ST.NUMPHASE(2);
EXC_ST.NUMPBLOCK(1)=EXC_ST.NUMEBLOCK(1);
EXC_ST.NUMPBLOCK(2)=0;
EXC_ST.NUMTBLOCK(1)=EXC_ST.NUMPBLOCK(1)+EXC_ST.NUMEBLOCK(1);
EXC_ST.NUMTBLOCK(2)=EXC_ST.NUMPBLOCK(2)+EXC_ST.NUMEBLOCK(2);
EXC_ST.NUMPTRIAL(1)=SEQ_ST.MAXSEQLENGTH*2;
EXC_ST.NUMPTRIAL(2)=0;
EXC_ST.NUMETRIAL(1)=SEQ_ST.MAXSEQLENGTH*5;%SEQ_ST.MAXSEQLENGTH*5
EXC_ST.NUMETRIAL(2)=SEQ_ST.MAXSEQLENGTH*5;%SEQ_ST.MAXSEQLENGTH*5
EXC_ST.ERRORCUT(1)=[1];%error % must be less than this % value
EXC_ST.ERRORCUT(2)=[100];%always move on (no practice)
EXC_ST.NUMRUNS=EXC_ST.NUMEBLOCK./size(SEQ_ST.SEQPAIRS,2);%3 unique sequences per run!

%TIMING
TIME_ST.RCIT(1)=0.800;
TIME_ST.RCIT(2)=4.00;
TIME_ST.RETRIEVAL(1)=0.700;
TIME_ST.RETRIEVAL(2)=1.000;
TIME_ST.CUTBASE(1)=0.700;
TIME_ST.CUTBASE(2)=1.00;
TIME_ST.TRIALT_P=0.700;
TIME_ST.ADJT(1)=1;
TIME_ST.ADJT(2)=0.25;%reduce time for first cycle!
% TIME_ST.INSTT=2.50;%fMRI only, self-pace for BEH
TIME_ST.FBT=1.00;%fMRI only, self-pace for BEH

%STIMULUS SETTING
STIM_ST.RESOLUTION=[1024,768];
STIM_ST.FONTSIZE=22;
STIM_ST.SIDEOFF=200;
STIM_ST.FRAMESIZE=6;
STIM_ST.SQSIZE=300;
STIM_ST.PLSIZE=80;
STIM_ST.FCADJ=-50;
STIM_ST.STIMRECT=ceil([0 0 STIM_ST.SQSIZE*1.43 STIM_ST.SQSIZE]);
STIM_ST.PLACERECT=ceil([0 0 STIM_ST.PLSIZE*1.43 STIM_ST.PLSIZE]);

% COLORS
COLOR_ST.BLACK = [0,0,0];
COLOR_ST.GRAY = [180,180,180];
COLOR_ST.WHITE = [250,250,250];
COLOR_ST.BACKG = COLOR_ST.GRAY;
COLOR_ST.RED = [255,0,0];
COLOR_ST.GREEN = [0,250,0];
COLOR_ST.BLUE = [0,0,250];
COLOR_ST.PURPLE = [155,0,250];
COLOR_ST.YELLOW=[250,250,0];
COLOR_ST.ORANGE = [250,165,0];
COLOR_ST.PINK = [250,0,250];
COLOR_ST.LIGHTBLUE = [0,250,250];
COLOR_ST.DARKBLUE = [50 90 170];
COLOR_ST.RULE(1,:)=COLOR_ST.RED;
COLOR_ST.RULE(2,:)=COLOR_ST.GREEN;
COLOR_ST.RULE(3,:)=COLOR_ST.PURPLE;
COLOR_ST.RULE(4,:)=COLOR_ST.BLUE;
COLOR_ST.RULE(5,:)=COLOR_ST.YELLOW;
COLOR_ST.RULE(6,:)=COLOR_ST.ORANGE;

% KEYS
KbName('UnifyKeyNames');
KEY_ST.KIDX(1)=-1;
KEY_ST.KIDX(2)=-1;
KEY_ST.TRIGGER=52;
KEY_ST.YKEY = KbName('y');
KEY_ST.NKEY = KbName('n');
KEY_ST.SPACEKEY = KbName('SPACE');
KEY_ST.ESCKEY = KbName('ESCAPE');
KEY_ST.SAMEKEY(1) = KbName('/?');KEY_ST.SAMEKEY(2) = 35;%left-index
KEY_ST.DIFFKEY(1) = KbName('z');KEY_ST.DIFFKEY(2) = 32;%right-indec

% DEFINE SESSIONS and OPEN FILES===========================================
EXC_ST.SUBID = input('Enter Subject ID: ','s');
EXC_ST.EXMODE = input('BEH(0) / FMRI(1) : ','s');
EXC_ST.RUNSTART = input('Starting Run : ');
RECORD_ST.SUBID=str2double(EXC_ST.SUBID);EXC_ST.EX=str2double(EXC_ST.EXMODE)+1;
RECORD_ST.DFILE = fopen([FPNAME_ST.EXFPATH EXC_ST.SUBID FPNAME_ST.TXTCODE '.txt'],'at');

rng(RECORD_ST.SUBID);

%SCREEN SETTING============================================================
% Open a graphics window on the main screen
% using the PsychToolbox's Screen function.
STIM_ST.SCCODE=max(Screen('Screens'));
% Screen('Resolution',STIM_ST.SCCODE,STIM_ST.RESOLUTION(1),STIM_ST.RESOLUTION(2));
Screen('Resolution',STIM_ST.SCCODE);%only after Psych version (3.0.14)
[W, ~]=Screen('OpenWindow',STIM_ST.SCCODE, 0,[],32,2);%[0 0 200 200]
Screen('FillRect',W,COLOR_ST.BACKG);

%Prevent Typing in the code!
rng(str2double(EXC_ST.SUBID),'twister');
ListenChar(1);
HideCursor;

%MXY = matrix of xy coordinates(6 pairs of 12 points)
%finds screen dimensions, returns center x,y
[swidth, sheight]=Screen('WindowSize',W);
STIM_ST.XC=fix(swidth/2);
STIM_ST.YC=fix(sheight/2);

Screen('TextFont',W,'Helvetica');
Screen('TextSize',W,STIM_ST.FONTSIZE);

%Instruction===============================================================
% locShowMessege('BEGINING');
%READY?

%Setting up randomized array for a pracblock===============================
locLoadStim;
locGenHierarchy;
locSetUpEx;
KEY_ST.SDKEYS=[KEY_ST.DIFFKEY(EXC_ST.EX),KEY_ST.SAMEKEY(EXC_ST.EX)];

%YAY!

% Calculate the starting blocl/run
numBinRun=EXC_ST.NUMEBLOCK(EXC_ST.EX)/EXC_ST.NUMRUNS(EXC_ST.EX);
startB=1+(numBinRun*(EXC_ST.RUNSTART-1));

%Main Loop Structure=======================================================
for BLOCK=[startB:EXC_ST.NUMTBLOCK(EXC_ST.EX)]
    %Message at the begining every block
    PREX.GO=0;PREX.REP=1;
    PREX.NUMTRIAL=[EXC_ST.NUMPTRIAL(EXC_ST.EX),EXC_ST.NUMETRIAL(EXC_ST.EX)];
    PREX.TYPE=PARAM_ST.PRACTICE(BLOCK)+1;
    RUNID=ceil(BLOCK/size(SEQ_ST.SEQPAIRS,2));
    
    %Set up intro
    %BINTRO is regular intro for each block
    %RUNINTRO resets the trigger time!!
    switch EXC_ST.EX
        case {1},ref=BLOCK;hitV=0;introT='BINTRO';
        case {2},ref=EXC_ST.NUMRUNS(EXC_ST.EX);hitV=1;introT='RUNINTRO';
    end
    if mod(BLOCK,ref)==hitV,locShowMessege(introT);end
    
    %Timing for TaskTS and MusicTS
    BLOCKTS=GetSecs;
    TKeeper=GetSecs;
    CCTime=GetSecs-TKeeper;
    %locCodeEvent('TEMP',10+PREX.TYPE);%begining of block @@@@@@@@@@@@@@@@@
    RECORD_ST.EV_BLOCK(BLOCK,:)=BLOCKTS-TRIGERTS;%record timing for re-epoch!
    
    %Repeatedly loop trials for practice, and one loop for experiment!!
    while ~PREX.GO
        %Display rule first
        locDrawFix('NEUTRAL');
        locDrawRule('BRULE');
        Screen('Flip',W);
        %if EXC_ST.EX==1,while KbCheck;end;KbWait;else,WaitSecs(TIME_ST.INSTT);end;
        while KbCheck(KEY_ST.KIDX(EXC_ST.EX));end;
        locGetTrigger('WAIT');%self-paced

        for TRIAL=1:PREX.NUMTRIAL(PREX.TYPE)
            %Within trial computations
            RECORD_ST.EV_TRIAL(BLOCK,TRIAL)=GetSecs-BLOCKTS;
            [CCTime,TKeeper] = locDoTrial(CCTime,TKeeper);
            
            %Recording varibales of interest
            locKeepRecord(RECORD_ST,'VARS');
            
            %Exit out if esc is pressed
            [~,~,keyCode]=KbCheck;
            if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;
        end;
        
        %Clean up RECORD_ST??
        RECORD_ST.EV_BLEND(BLOCK,:)=GetSecs-TRIGERTS;
        
        %Block result (NEED THIS???)
        if BLOCK < EXC_ST.NUMTBLOCK(EXC_ST.EX) && EXC_ST.EX==1,locShowMessege('BRESULT');end;
        
        %End Screen
        if BLOCK == EXC_ST.NUMTBLOCK(EXC_ST.EX),locShowMessege('EXEND');break;end;
    end
    
end

%Putting everything back
Screen('CloseAll');
ShowCursor;
ListenChar;

end

%##########################################################################

function locSetUpEx
%locSetUpEx

%**************************************************************************
%BLOCK ORGANIZATION/%TIME REGULATION/%STIMULUS/%KEYS/%FONT/%COLOR/%PARAMETER
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST
%**************************************************************************

%PARAMST SETTING===========================================================
PARAM_ST.PARAMETER={
    'PRACEX';
    'SEQSENCE';%set up sequence properties at the block-level
    'CONTENTS';%code parameters coding contents and position of 2-levels
    'TASK';%position of stimulus (1=top_l,2=top_r,3-bot_r,4=bot_l)
    };

%RANDOMIZATION LOOP========================================================
%Parameter Sign
for p=1:length(PARAM_ST.PARAMETER)
    
    %Parameter Setting!
    fn=PARAM_ST.PARAMETER{p,1};
    numPB=EXC_ST.NUMPBLOCK(EXC_ST.EX);numPT=EXC_ST.NUMPTRIAL(EXC_ST.EX);
    numEB=EXC_ST.NUMEBLOCK(EXC_ST.EX);numET=EXC_ST.NUMETRIAL(EXC_ST.EX);
    numTB=EXC_ST.NUMTBLOCK(EXC_ST.EX);numPH=EXC_ST.NUMPHASE(EXC_ST.EX);
    
    switch fn
        case {'SEQSENCE'}
            %Practice:0=practice 1=experiment
            blType(1)={[0,1]};blType(2)={[1]};
            PARAM_ST.('PRACTICE')=BalanceTrials(numTB,0,blType{EXC_ST.EX});
            
            %Sequence Type
            seqtype=cell(numPH,1);seqP=num2cell(SEQ_ST.SEQPAIRS,2);
            for s=1:numPH
                c=cellfun(@(s)Shuffle(s)',seqP,'Uni',false);%Shuffle within cmbs
                c=Shuffle(c);seqtype(s)={vertcat(c{:})};%Shuffle across cmbs
            end
            seqtype=repmat(vertcat(seqtype{:}),1,length(blType{EXC_ST.EX}));
            PARAM_ST.('SEQTYPE')=reshape(seqtype',1,[])';
            
        case {'CONTENTS'}
            %Define parameters below for later analysis....!
            zeroD=zeros(numTB,max(numET));
            PARAM_ST.CHUNK=zeroD;%CHUNK = identity of chunk
            PARAM_ST.ACPOS=zeroD;%ACPOS = position across-chunks
            PARAM_ST.ELEMENT=zeroD;%ELEMENT = task (vertical,horizontal,diagonal)
            PARAM_ST.WCPOS=zeroD;%WCPOS = position within-chunk
            numC=SEQ_ST.NUMCHUNK;
            chp=SEQ_ST.CHUNKPAIRS;
            flatL=@(d)reshape(d,1,[]);
            
            %Loop through each block
            for b=1:numTB
                %Change where to read for numtrial,chunk codde etc!
                seqC=chp(PARAM_ST.SEQTYPE(b),:);
                
                %CHUNK
                chunkC=flatL(repmat(seqC,SEQ_ST.MAXSEQLENGTH/numC,1));
                PARAM_ST.CHUNK(b,:)=repmat(chunkC,[1,numET/length(chunkC)]);
                
                %ACPOS
                acposC=flatL(repmat(1:numC,SEQ_ST.MAXSEQLENGTH/numC,1));
                PARAM_ST.ACPOS(b,:)=repmat(acposC,[1,numET/length(acposC)]);
                
                %ELEMENT
                elmC=SEQ_ST.('BASESEQUENCE')(PARAM_ST.SEQTYPE(b),:);
                elmC=horzcat(elmC{:});
                PARAM_ST.ELEMENT(b,:)=repmat(elmC,[1,numET/length(elmC)]);
                
                %WCPOS
                wcposC=1:SEQ_ST.MAXSEQLENGTH/numC;
                PARAM_ST.WCPOS(b,:)=repmat(wcposC,[1,numET/max(wcposC)]);
                
                %BSPOS
                bsposC=1:SEQ_ST.MAXSEQLENGTH;
                PARAM_ST.BSPOS(b,:)=repmat(bsposC,[1,numET/max(bsposC)]);
            end
            
        case {'TASK'}
            % Match/mismatch probe
            prob=SEQ_ST.ONTRACKRATE;
            PARAM_ST.('ONTRACK')=reshape(CoinFlip(numTB*numET,prob),numTB,numET);
            
            % Test probe
            elm=SEQ_ST.BASE(1,:);
            probe=zeros(numTB,numET);
            ontIDX=logical(PARAM_ST.('ONTRACK'));
            probe(ontIDX)=PARAM_ST.ELEMENT(ontIDX);
            probe(~ontIDX)=cell2mat(cellfun(@(p)Sample(setdiff(elm,p)),...
                num2cell(PARAM_ST.ELEMENT(~ontIDX)),'Uni',false));
            PARAM_ST.('PROBE')=probe;
            
            %Correct response
            refK=[KEY_ST.DIFFKEY(EXC_ST.EX);KEY_ST.SAMEKEY(EXC_ST.EX);];
            PARAM_ST.('CORRESP')=refK(PARAM_ST.ONTRACK+1);
    end
end

%RECORD SETTING============================================================
% Add BLOCK and TRIAL to PARAM_ST to make it easier to Join!!
cellD=cell(numTB,numET);
zeroD=zeros(numTB,numET);
cVars={};
zVars={'SUBID','BLOCK','TRIAL','REP','RUNID','RT','RESP','ACC','RSIJIT',...
    'EV_START','EV_BLOCK','EV_TRIAL','EV_RETRIEVAL','EV_PROBE','EV_RESP','EV_BLEND'};
for p=zVars(2:end),if ~isfield(RECORD_ST,p),RECORD_ST.(p{:})=zeroD;end;end
% for p=cVars,if ~isfield(RECORD_ST,p),RECORD_ST.(p{:})=cellD;end;end
RECORD_ST.PARAMORDER=[zVars,cVars];
RECORD_ST.ACC=nan(size(zeroD));%For rule-relearning!
RECORD_ST.SUBID=str2double(EXC_ST.SUBID);

% Write header
locKeepRecord(RECORD_ST ,'header');

%RECORD PARAMST JUST IN CASE!!=============================================
paramtMat=[EXC_ST.SUBID '_' FPNAME_ST.TXTCODE '_PARAM_ST.mat'];
save(paramtMat,'PARAM_ST');
disp('YAY! Completed!!');

end

%##########################################################################

function [CCTime,TKeeper] = locDoTrial(CCTime,TKeeper)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST PREX
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL BTRIAL BLOCKTS TRIGERTS
%**************************************************************************

%PHASE1:Trial Preparation--------------------------------------------------
[~,~,keyCode]=KbCheck;
if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;

%Cover probe
Screen('FillRect', W,COLOR_ST.BACKG);
if RECORD_ST.ACC(BLOCK,max(TRIAL-1,1))==0,locDrawRule('RELEARN');end
Screen('Flip',W);

%Jittering RSI-------------------------------------------------------------
% locCodeEvent('TEMP',9);%begining of trial @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% locCodeEvent('CHAR',[]);
RECORD_ST.EV_TRIAL(BLOCK,TRIAL)=GetSecs-TRIGERTS;%record timing for re-epoch!
timeF=TIME_ST.RCIT(EXC_ST.EX);
if TRIAL<=SEQ_ST.MAXSEQLENGTH,tADJ=TIME_ST.ADJT(EXC_ST.EX);else tADJ=1;end
rsiJit=timeF*tADJ;%no jittering
% rsiJit=Sample((timeF*0.75):0.01:(timeF*1.25));jittered from 75~125%
WaitSecs(rsiJit);RECORD_ST.RSIJIT(BLOCK,TRIAL)=rsiJit*1000;

%Self-pasing
%locDrawFix('PRETRIAL');
%Screen('Flip',W);
%locGetInitiation;
%locCodeEvent('TEMP',1);%self-started trial onset@@@@@@@@@@@@@@@@@@@@@@@@@@

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PHASE1:FIXATION DISPLAY & FIXED RSI INTERVAL------------------------------
%DO NOT Pause until subject press space bar!(trial starts automatically)
locDrawFix('NEUTRAL');
if RECORD_ST.ACC(BLOCK,max(TRIAL-1,1))==0,locDrawRule('RELEARN');end
Screen('Flip',W);
% locCodeEvent('TEMP',1);%fixed pre-stimulus period onset@@@@@@@@@@@@@@@@@@
RECORD_ST.EV_RETRIEVAL(BLOCK,TRIAL)=GetSecs-TRIGERTS;%record timing for re-epoch!
WaitSecs(TIME_ST.RETRIEVAL(EXC_ST.EX));

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PHASE2:STIMULUS DISPLAY & RECORDING---------------------------------------
%Drawing a test probe
locDrawStim;
if RECORD_ST.ACC(BLOCK,max(TRIAL-1,1))==0,locDrawRule('RELEARN');end
Screen('Flip',W);
% locCodeEvent('TEMP',2);%stimulus onset@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
RECORD_ST.EV_PROBE(BLOCK,TRIAL)=GetSecs-TRIGERTS;%record timing for re-epoch!

%Measureing response!
locGetResponse;
% locCodeEvent('RESP',[]);%resp onset@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
RECORD_ST.EV_RESP(BLOCK,TRIAL)=GetSecs-TRIGERTS;%record timing for re-epoch!

%Showing feedback
if RECORD_ST.ACC(BLOCK,TRIAL)==0
    locDrawFix('NEUTRAL');
    locDrawStim;
    locDrawRule('RELEARN');
    Screen('Flip',W);
    if EXC_ST.EX==1,while KbCheck(KEY_ST.KIDX(EXC_ST.EX));end;locGetTrigger('WAIT');else,WaitSecs(TIME_ST.FBT);end
end

%Event Code :End of trial and Characteristic codes@@@@@@@@@@@@@@@@@@@@@@@@@
%--------------------------------------------------------------------------


end

%##########################################################################

% function locCodeEvent(codeT,code)
% %==========================================================================
% % locCodeEvent creates three different kind of event codes
% % temporal event codes:codes for individual events regardress of contents
% % characteristic codes:codes for the contents of one trial(all conditions)
% %==========================================================================
% %STIMULUS PROPERTY / TIME REGULATION / COLOR
% global BLOCK TRIAL BTRIAL RECORD_ST
% global EXC_ST SEQ_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST COLOR_ST KEY_ST
% %**************************************************************************
%
% switch codeT
%     case {'TEMP'},%Temporal codes
%         portCode=FPNAME_ST.PORTID_S{:};eventCode={code};
%
%     case {'CHAR'},%Characteristic codes
%         portCode=FPNAME_ST.PORTID_S{:};
%
%         %-Block (200~): (block+200)
%         eventCode(1,:)={BLOCK+200};
%
%         %-Trial (100~): (trial+100)
%         eventCode(2,:)={TRIAL+100};
%
%     case {'RESP'},%Response codes
%         portCode=FPNAME_ST.PORTID_R{:};
%
%         %-Response: 90(different),111 (same)
%         eventCode(1,:)={RECORD_ST.RESP(BLOCK,TRIAL)};
% end
%
% %Actually write event code here!
% if srtcmp(EXC_ST.EXMODE,'1'),cd(FPNAME_ST.PORDIR{:});end
% for i=1:size(eventCode,1)
%     switch EXC_ST.EXMODE
%         case {'1'},write_parallel(portCode,eventCode{i,:});
%         otherwise,disp(strcat(portCode,'_',num2str(eventCode{i})));
%     end
% end
% if srtcmp(EXC_ST.EXMODE,'1'),cd('../');end;
%
% end.


%##########################################################################

function locGetTrigger(workT)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST SEQ_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST RECORD_ST
global AUDIO_ST COLOR_ST KEY_ST
global BLOCK TRIAL BTRIAL BLOCKTS TRIGERTS W
%**************************************************************************
% Trigger using apostrophe
goF = false;

while ~goF
    [~,~,keyC]=KbCheck(KEY_ST.KIDX(EXC_ST.EX));
    
    %Switch which keys to listen to
    switch workT
        case {'TRIGGER'},keyW=KEY_ST.TRIGGER;
        case {'WAIT'},keyW=KEY_ST.SDKEYS;
    end

    if any(ismember(find(keyC),keyW)),goF=1;TRIGERTS=GetSecs;end;
    if EXC_ST.EX==1 && any(keyC),goF=1;TRIGERTS=GetSecs;end;
    if keyC(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;
end

% Record only trigger time
if strcmp(workT,'TRIGGER'),RECORD_ST.EV_START(:,:)=TRIGERTS;end;


end

%##########################################################################

function locGetResponse
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST PREX RUNID
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL BLOCKTS
%**************************************************************************

%Start measuring time------------------------------------------------------
startT = GetSecs;
coResp = PARAM_ST.CORRESP(BLOCK,TRIAL);
keyCF = true;
keyRF = false;
if ~PARAM_ST.PRACTICE(BLOCK),cutT=TIME_ST.TRIALT_P;else,cutT=TIME_ST.CUTBASE(EXC_ST.EX);end

% Loop while there is time (Loop for one trial)
while (GetSecs - startT < cutT)
    %Check keyboard when KeyCF is true
    if keyCF,[keyP, secs, keyC] = KbCheck(KEY_ST.KIDX(EXC_ST.EX));end
    if keyC(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar(0);break;end;
    if any(ismember(find(keyC),KEY_ST.SDKEYS)),keyRF=true;end;%only listen to correct keys

    %BEH only
    if EXC_ST.EX==1 && keyP==1 && length(find(keyC))==1 && keyRF,break;end;
    
    %Wait until cutT (test phase interval) passes for FMRI!
    if EXC_ST.EX==2 && keyP==1 && length(find(keyC))==1 && keyRF,keyCF=false;end;
end

%Record some variables!!
RECORD_ST.BLOCK(BLOCK,TRIAL)=BLOCK;
RECORD_ST.TRIAL(BLOCK,TRIAL)=TRIAL;
RECORD_ST.REP(BLOCK,TRIAL)=PREX.REP;
RECORD_ST.RUNID(BLOCK,TRIAL)=RUNID;
RECORD_ST.RT(BLOCK,TRIAL)=(secs-startT)*1000;
RECORD_ST.TASKTS(BLOCK,TRIAL)=(secs-BLOCKTS)*1000;
if isempty(find(keyC, 1)),keyC=99;end;
RECORD_ST.RESP(BLOCK,TRIAL)=find(keyC);
if coResp==find(keyC),acc=1; else, acc=0;end
RECORD_ST.ACC(BLOCK,TRIAL)=acc;

end

%##########################################################################

function locDrawFix(workT)
%==========================================================================
% locDrawFix draws a very nice fixation cross in the center of screen
%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST SEQ_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST RECORD_ST
global AUDIO_ST COLOR_ST KEY_ST IMAGE_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL BTRIAL BLOCKTS W
%**************************************************************************

% Fixation type:Cross
fixation='+';
Screen('TextFont',W,'Courier');
switch workT
    case {'PRETRIAL'},color=COLOR_ST.RED;fsize=STIM_ST.FONTSIZE+20;
    case{'NEUTRAL'},color=COLOR_ST.BLACK;fsize=STIM_ST.FONTSIZE+20;
end
Screen('TextSize',W,fsize);
DrawFormattedText(W,fixation,'center',STIM_ST.YC+STIM_ST.FCADJ,color);
Screen('TextSize',W,STIM_ST.FONTSIZE);%put back the textsize
Screen('TextFont',W,'Helvetica'); %dg
end

%##########################################################################

function locDrawRule(workT)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST SEQ_ST FPNAME_ST STIM_ST  TIME_ST COLOR_ST IMAGE_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST RECORD_ST BLOCK TRIAL PRACEXC
%**************************************************************************

% Drawing
for s=1:SEQ_ST.MAXSEQLENGTH
    %Specify element and XY cordinates!
    elm=PARAM_ST.ELEMENT(BLOCK,s);xadj=0;if s>SEQ_ST.MAXSEQLENGTH/2,xadj=30;end
    xp=STIM_ST.XC+((STIM_ST.PLSIZE+30)*(s-1))-290+xadj;yp=STIM_ST.YC+180;
    
    %Element----------------------->
    im=CenterRectOnPoint(STIM_ST.PLACERECT,xp,yp);
    Screen('DrawTexture',W,IMAGE_ST.IMARRAY{elm,3},[],im);
    
    %Frame------------------------->
    color=COLOR_ST.BLACK;seqL=SEQ_ST.MAXSEQLENGTH;
    if strcmp(workT,'RELEARN') && isequal(s,mod(TRIAL-1,seqL)+1),color=COLOR_ST.RED;end;
    if strcmp(workT,'BRULE'),color=COLOR_ST.BLUE;end
    fm=CenterRectOnPoint([0,0,STIM_ST.PLSIZE+25 STIM_ST.PLSIZE+25],xp,yp);
    Screen('FrameRect', W,color,fm,STIM_ST.FRAMESIZE);
    %Screen('DrawText',W,num2str(s),xp-5,yp-80,COLOR_ST.WHITE);
    
    %Message----------------------->
    %if EXC_ST.EX==1 && strcmp(workT,'BRULE')
    %    expCond={'Practice','Experiment'};
    %    msg=sprintf('%s ',expCond{PARAM_ST.PRACTICE(BLOCK)+1});
    %    DrawFormattedText(W,msg,'center',STIM_ST.YC+325,COLOR_ST.WHITE);
    %end
end

end

%##########################################################################

function locDrawStim
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST BLOCK TRIAL IMAGE_ST
%**************************************************************************

% Draw one of the task
probe=PARAM_ST.PROBE(BLOCK,TRIAL);
im=CenterRectOnPoint(STIM_ST.STIMRECT,STIM_ST.XC,STIM_ST.YC+STIM_ST.FCADJ);
Screen('DrawTexture',W,IMAGE_ST.IMARRAY{probe,3},[],im);

% Debug
% Screen('DrawText',W,num2str(PARAM_ST.ONTRACK(BLOCK,TRIAL)),STIM_ST.XC,STIM_ST.YC-250,COLOR_ST.WHITE);

end

%##########################################################################

function locShowMessege(msgtime)
%==========================================================================
% locShowMessege displays all of string based messages for a variety of
% purposes depending on msgtime in combination with locCondMessage
%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST BLOCK TRIAL PRACEXC PREX RUNID
%**************************************************************************

%Consistant settings
Screen('FillRect',W,COLOR_ST.BACKG);
Screen('TextSize',W,STIM_ST.FONTSIZE);

switch msgtime
    case {'BINTRO'}% Begining Block BREAK screen
        [ConText]=locCondMessage(msgtime);
        TEXT(1)={ConText{1}};
        TEXT(2)={ConText{2}};
        TEXT(3)={''};
        TEXT(4)={'Be as accurate and fast as possible'};
        TEXT(5)={'Press any key to proceed'};
        C={'RED','BLACK','BLACK','BLACK','BLACK'}; %dg
        
    case {'RUNINTRO'}% Begining of run
        [ConText]=locCondMessage(msgtime);
        TEXT(1)={ConText{1}};
        TEXT(2)={ConText{2}};
        TEXT(3)={''};
        TEXT(4)={''};
        TEXT(5)={'Please wait for the instructions'};
        C={'BLACK','BLACK','BLACK','BLACK','BLACK'}; %dg
                
    case {'BRESULT'}% Ending Block BREAK screen (no money)
        [ConText]=locCondMessage(msgtime);
        TEXT(1)={['You have completed block ',num2str(BLOCK),' of ',num2str(EXC_ST.NUMTBLOCK(EXC_ST.EX))]}; %dg
        TEXT(2)={ConText{1}};
        TEXT(3)={ConText{2}};
        TEXT(4)={''};
        TEXT(5)={'Press any key to proceed'};
        C=repmat({'BLACK'},1,size(TEXT,2)); %dg
        
    case {'EXEND'}% Ending Experiment screen
        TEXT(1)={'Thank you'};
        TEXT(2)={'You have completed the experiment'};
        TEXT(3)={''};
        TEXT(4)={''};
        TEXT(5)={'Press any key to proceed'};
        C=repmat({'BLACK'},1,size(TEXT,2)); %dg
end

%Drawing
for i=1:size(TEXT,2)
    if strcmp(msgtime,'BRULE'),adj=50;else adj=100;end;
    text=TEXT{i};xp=STIM_ST.XC-350;yp=STIM_ST.YC-200+(adj*(i-1));
    if i==size(TEXT,2),yp=STIM_ST.YC-200+400;end
    Screen('DrawText',W,text,xp,yp,COLOR_ST.(C{i}));
end
Screen('Flip',W);

%Prevents jumping due to key pressed ahead
%and waits for keyboard input
switch msgtime
    case {'RUNINTRO','BINTRO'},locGetTrigger('TRIGGER');
    %otherwise,if EXC_ST.EX ==1,while KbCheck;end;KbWait;else,WaitSecs(TIME_ST.INSTT);end
    otherwise,while KbCheck(KEY_ST.KIDX(EXC_ST.EX));end;locGetTrigger('WAIT');%always self-paced instruction
end

%Exit out if esc is pressed
[~,~,keyCode]=KbCheck;
if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;

end

%##########################################################################

function [ConText]=locCondMessage(msgtime)
%==========================================================================
% locShowMessege displays all of string based messages for a variety of
% purposes depending on msgtime in combination with locCondMessage
%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST BLOCK TRIAL PRACEXC PREX RUNID
%**************************************************************************

switch msgtime
    case {'BINTRO'}% Begining Block BREAK screen
        expCond={'PRACTICE!!','EXPERIMENT!!'};
        ConText(1)={sprintf('%s',expCond{PREX.TYPE})};
        ConText(2)={'Please remember a sequence of orientations'};
        
    case {'RUNINTRO'}% Begining Block BREAK screen
        ConText(1)={sprintf('Session %d is starting!',RUNID)};
        ConText(2)={'Please cycle a sequence of orientations'};
        
    case {'BRESULT'}% Begining Block BREAK screen
        %Change Setting for Practice and Experiment
        trial=PREX.NUMTRIAL(PREX.TYPE);PREX.GO=true;
        errorScore=length(find(RECORD_ST.ACC(BLOCK,1:trial)==0));
        errorCheck=(errorScore/trial*100)<=EXC_ST.ERRORCUT(PREX.TYPE);
        ConText(1)={sprintf('You missed %d out of %d trials',errorScore,trial)};
        ConText(2)={sprintf('Please go on to the next block!')};
        
        %Evaluation for practice blocks
        if PREX.TYPE==1 && ~errorCheck
            PREX.GO=false;PREX.REP=PREX.REP+1;
            ConText(2)={'Too many mistakes!! Please repeat practice!!'};
            
            %Re-randomize Stimpos & correct answers
            numPT=EXC_ST.NUMPTRIAL(EXC_ST.EX);stIDX=RandSample(1:numPT);
            PARAM_ST.ONTRACK(BLOCK,:)=circshift(PARAM_ST.ONTRACK(BLOCK,:),[0,-stIDX+1]);
            PARAM_ST.CORRESP(BLOCK,:)=circshift(PARAM_ST.CORRESP(BLOCK,:),[0,-stIDX+1]);
            elm=SEQ_ST.BASE(1,:);probe=zeros(1,numPT);
            ontIDX=logical(PARAM_ST.ONTRACK(BLOCK,:));
            probe(ontIDX)=PARAM_ST.ELEMENT(BLOCK,ontIDX);
            probe(~ontIDX)=cell2mat(cellfun(@(p)Sample(setdiff(elm,p)),...
                num2cell(PARAM_ST.ELEMENT(BLOCK,~ontIDX)),'Uni',false));
            PARAM_ST.('PROBE')(BLOCK,:)=probe;
            
            %Erase data from this attempt
            RECORD_ST.RT(BLOCK,:)=0;RECORD_ST.ACC(BLOCK,:)=NaN;RECORD_ST.RESP(BLOCK,:)=0;
        end
end

end

%##########################################################################

function locKeepRecord(dataArray,workT)
%==========================================================================
% locKeepRecord records the title of text files for the final output.
% dataArray needs to contain VARLABELS
% VARLABELS = specify all variables needs to be stored from dataArray
% NOTE:Currently, it only records the first column of info from one
% column dimension of the dataArray!!!!(may change later?)
%==========================================================================

%**************************************************************************
%BLOCK ORGANIZATION/%TIME REGULATION/%STIMULUS/%KEYS/%FONT/%COLOR/%PARAMETER
global RECORD_ST EXC_ST TASKST STIMST COLORST KEYST PARAM_ST BLOCK TRIAL
%**************************************************************************

% Setting ahead
if ~isstruct(dataArray),error('locKeepRecord assums struct array!');end;
ActArray=cell(length(dataArray.PARAMORDER),1);

for g=1:length(dataArray.PARAMORDER)
    %Is it header writing or datawriting
    switch workT
        case {'header'},var=dataArray.PARAMORDER(g);
        otherwise,var=dataArray.(dataArray.PARAMORDER{g});
    end
    
    %Checking dimensions of data
    %ECOND=static through exp,BCOND=block wise rand,TCOND=trial wise rand!
    if size(var,1)==1,varType='ECOND';ref=[1,1];end;
    if size(var,1)==EXC_ST.NUMTBLOCK(EXC_ST.EX),varType='BCOND';ref=[BLOCK,1];end
    if size(var,2)==EXC_ST.NUMETRIAL(EXC_ST.EX),varType='TCOND';ref=[BLOCK,TRIAL];end;
    
    %Put Data First
    if iscell(var),value=var{ref(1)};else value=var(ref(1),ref(2));end;ActArray(g,1)={value};
    
    switch ischar(ActArray{g,1})
        case {1},dataT='s';
        case {0},if mod(ActArray{g,1},1)==0,dataT='d';else dataT='7.4f';end;
    end
    ActArray(g,2)={dataT};
    
    if g==length(ActArray),ActArray(g,3)={'n'};else ActArray(g,3)={'t'};end
    fprintf(RECORD_ST.DFILE,strcat('%',ActArray{g,2},'\',ActArray{g,3}),ActArray{g,1});
end

end

%##########################################################################

function  locGenHierarchy
%==========================================================================
% SEQUENCES FROM ORIGINAL PILEHIGH (2 task-elements)-----------------------
% TASK ELEMENT = vertical (1), horizontal (2), diagonal (3)
% TASK CUHNK = chunk 1 (1,2,3), chunk 2 (2,3,1), chunk3 (3,1,2)
% SEQUENCE = 6 different combination of TASK CHUNKS
% EX) chunk1-chunk2-chunk3, chunk1-chunk2-chunk3,
%==========================================================================
%STIMULUS PROPERTY / TIME REGULATION / COLOR
global EXC_ST SEQ_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST COLOR_ST KEY_ST
%**************************************************************************

% Setting Rule of TaskElenent!!++++++++++++++++++++++++++++++++++++++++++++
% Creating template of hieararchy!!++++++++++++++++++++++++++++++++++++++++
baseseq=cell(size(SEQ_ST.CHUNKPAIRS,1),1);
for c=1:size(SEQ_ST.CHUNKPAIRS,1)
    c1=SEQ_ST.BASE(SEQ_ST.CHUNKPAIRS(c,1),:);
    c2=SEQ_ST.BASE(SEQ_ST.CHUNKPAIRS(c,2),:);
    baseseq(c)={horzcat({c1},{c2})};
end

SEQ_ST.BASESEQUENCE=vertcat(baseseq{:});

end

%##########################################################################

function locLoadStim
%==========================================================================
% locLoadStim loads up all of visual and auditory data. Right now
% compatible for jpg and wav data. Provide the directory path for the file
% location and exact file names, and it will create a summary cell array.
%==========================================================================
global FPNAME_ST AUDIO_ST IMAGE_ST RECORD_ST STIM_ST W
%**************************************************************************

%Image File Names
%1=filename
%2=imagedata
IMAGE_ST.IMFPATH=strcat(FPNAME_ST.EXFPATH,'OrderYourMind_Images');
IMAGE_ST.IMARRAY=cell(3,3);
IMAGE_ST.IMARRAY(1,1)={'Bar_1.jpg'};
IMAGE_ST.IMARRAY(2,1)={'Bar_2.jpg'};
IMAGE_ST.IMARRAY(3,1)={'Bar_3.jpg'};

% Visual stimuli preparation
for i=1:size(IMAGE_ST.IMARRAY,1)
    %reading jpg image data
    imaCode = char([IMAGE_ST.IMFPATH filesep IMAGE_ST.IMARRAY{i,1}]);
    IMAGE_ST.IMARRAY{i,2}=imread(strtrim(imaCode));
    IMAGE_ST.IMARRAY{i,3}=Screen('MakeTexture',W,IMAGE_ST.IMARRAY{i,2});
end

end

