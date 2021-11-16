function OrderYourRule
%==========================================================================
% PURPOSE:
% TO DO:
% SETTING:
% 1, Completely fixed RSI period (blink few as possible)
% 2, 700 ms RSI, 850 ms cutoff
%==========================================================================
clear all;
close all;
sca;
if 1, Screen('Preference', 'SkipSyncTests', 1); end
commandwindow;while KbCheck; end

%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST TASK_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL PRACEXC BTRIAL BLOCKSTS PREX
%**************************************************************************

%1=V,2=H,3=C,4=CC
SEQ_ST.NUMBASEPAIR=3;
SEQ_ST.BASE(1,:)=[1,2,3];
SEQ_ST.BASE(2,:)=[1,3,2];
SEQ_ST.BASE(3,:)=[2,1,3];
SEQ_ST.BASE(4,:)=[2,3,1];
SEQ_ST.BASE(5,:)=[3,1,2];
SEQ_ST.BASE(6,:)=[3,2,1];

% Manual code of chunls:Fixed probs at 1st level
chunkpairs(1,:)=[1,2];chunkpairs(2,:)=[1,3];
chunkpairs(3,:)=[2,1];chunkpairs(4,:)=[2,5];
chunkpairs(5,:)=[3,4];chunkpairs(6,:)=[3,1];
chunkpairs(7,:)=[4,3];chunkpairs(8,:)=[4,6];
chunkpairs(9,:)=[5,6];chunkpairs(10,:)=[5,2];
chunkpairs(11,:)=[6,5];chunkpairs(12,:)=[6,4];
SEQ_ST.CHUNKPAIRS=chunkpairs;
SEQ_ST.NUMSEQTYPE=size(chunkpairs,1);
SEQ_ST.MAXLENGTH=[6];%smaximum elements in sequence
SEQ_ST.NUMCHUNKS=[2];%Number of Chunks

% TASK SETTING
% pos  |% vert|% hori|% clock|% counter|
% 1,2  |% 4,3 |% 2,1 |% 2,3  |% 4,1|
% 4,3  |% 1,2 |% 3,4 |% 1,4  |% 3,2|
TASK_ST.TASK={'vertical','horizontal','diagonal'};
TASK_ST.DESTREF(1)={[4,3;1,2]};
TASK_ST.DESTREF(2)={[2,1;3,4]};
TASK_ST.DESTREF(3)={[3,4;2,1]};
TASK_ST.DESTREF(4)={[1,2;4,3]};
TASK_ST.DREFIDX=length(TASK_ST.DESTREF);%last one is a reference
TASK_ST.ERRORCUT=[1,10];%error % must be less than this value
TASK_ST.TIMECUT=[1,1];%0=no cuff off, 1= cutoff
TASK_ST.RTCUT=[70];%percentile cutoff


%EXPERIMENT SETTING
EXC_ST.NUMEBLOCK=SEQ_ST.NUMSEQTYPE*3;
EXC_ST.NUMPBLOCK=EXC_ST.NUMEBLOCK;%1 for each number of chunks
EXC_ST.NUMTBLOCK=EXC_ST.NUMPBLOCK+EXC_ST.NUMEBLOCK;
EXC_ST.NUMPTRIAL=SEQ_ST.MAXLENGTH*1;%use (2) based on 2-chunk seqs
EXC_ST.NUMETRIAL=SEQ_ST.MAXLENGTH*5;

% Set up the timer
TIME_ST.YEILDT=0.00001;
TIME_ST.TRIALT_P=inf;
TIME_ST.RSIT=0.0;
TIME_ST.PREPT=0.700;
TIME_ST.CUTT=[0.850];%850


%FILE PATH SETTING
FPNAME_ST.TXTCODE='oyr';
FPNAME_ST.EXFILE='OrderYourRule.m';
FPNAME_ST.PTLPATH=PsychtoolboxRoot;
exFilePath=which(FPNAME_ST.EXFILE);
FPNAME_ST.EXFPATH=exFilePath(1:end-length(FPNAME_ST.EXFILE));
cd(FPNAME_ST.EXFPATH);
FPNAME_ST.PORTID_S={'DCC8'};
FPNAME_ST.PORTID_R={'DCD8'};
FPNAME_ST.PORDIR={strcat(pwd,'\PortCode_Supports')};

%STIMULUS SETTING
%General
STIM_ST.RESOLUTION=[1024,768];
STIM_ST.FONTSIZE=22;

%Basic obejects info
STIM_ST.FRAMESIZE=5;
STIM_ST.SQSIZE=300;
STIM_ST.PLSIZE=80;

%General Stimulus
STIM_ST.STIMSIZE=60;
STIM_ST.FRSIZE=250;
STIM_ST.FRAMESIZE=5;
STIM_ST.FRAMERECT=[0 0 STIM_ST.FRSIZE STIM_ST.FRSIZE];
STIM_ST.STIMRECT=ceil([0 0 STIM_ST.STIMSIZE STIM_ST.STIMSIZE]);
OFFP=60;OFFC=+5;
STIM_ST.OFFSETS(1,:)=[-OFFP,OFFP,OFFP,-OFFP];%offset for X
STIM_ST.OFFSETS(2,:)=[-OFFP,-OFFP,OFFP,OFFP];%offset for Y
% STIM_ST.OFFCOVER(1,:)=[-OFFC,-OFFC,OFFC,OFFC];
% STIM_ST.OFFCOVER(2,:)=[OFFC,-OFFC,-OFFC,OFFC];
if ismac,STIM_ST.FXADJ=5;else STIM_ST.FXADJ=-43;end %#ok<SEPEX>

% BASIC COLORS
COLOR_ST.BLACK = [0,0,0];
COLOR_ST.GRAY = [165,165,165];
COLOR_ST.WHITE = [250,250,250];
COLOR_ST.BACKG = COLOR_ST.GRAY;
% STIM COLORS
COLOR_ST.RED = [255,0,0];
COLOR_ST.GREEN = [0,250,0];
COLOR_ST.BLUE = [0,0,250];
COLOR_ST.PURPLE = [155,0,250];
COLOR_ST.YELLOW=[250,250,0];
COLOR_ST.ORANGE = [250,165,0];
COLOR_ST.PINK = [250,0,250];
COLOR_ST.LIGHTBLUE = [0,250,250];
COLOR_ST.DARKBLUE = [50 90 170];

% Set keys
KbName('UnifyKeyNames');
KEY_ST.T_LKEY = KbName('4');%4
KEY_ST.T_RKEY = KbName('5');%5
KEY_ST.B_LKEY=KbName('1');%1
KEY_ST.B_RKEY=KbName('2');%2
% % Testing
% KEY_ST.T_LKEY = KbName('i');%4
% KEY_ST.T_RKEY = KbName('o');%5
% KEY_ST.B_LKEY=KbName('k');%1
% KEY_ST.B_RKEY=KbName('l');%2
KEY_ST.SPACEKEY = KbName('SPACE');
KEY_ST.ESCKEY = KbName('ESCAPE');


%Get subject ID and open data file for output==============================
EXC_ST.SUBID = input('Enter Subject ID: ','s');
EXC_ST.EXMODE = input('EEG(1) / Behavioral(0): ','s');
fName=[fileroot(FPNAME_ST.EXFILE) [EXC_ST.SUBID FPNAME_ST.TXTCODE '.txt']];
RECORD_ST.DFILE=fopen(fName,'at');
if exist('rng','file'),rng(str2double(EXC_ST.SUBID),'twister');else rand('twister',str2double(EXC_ST.SUBID));end


%Screen Stuff==============================================================
% Open a graphics window on the main screen
% using the PsychToolbox's Screen function.
STIM_ST.SCCODE=max(Screen('Screens'));
Screen('Resolution',STIM_ST.SCCODE,STIM_ST.RESOLUTION(1),STIM_ST.RESOLUTION(2));
[W, wRect]=Screen('OpenWindow',STIM_ST.SCCODE, 0,[],32,2);%[0 0 200 200]
Screen('FillRect',W,COLOR_ST.BACKG);

%Prevent Typing in the code!
% ListenChar(2);
HideCursor;

%MXY = matrix of xy coordinates(6 pairs of 12 points)
%finds screen dimensions, returns center x,y
[swidth, sheight]=Screen('WindowSize',W);
STIM_ST.XC=fix(swidth/2);
STIM_ST.YC=fix(sheight/2);

Screen('TextFont',W,'Helvetica');
Screen('TextSize',W,STIM_ST.FONTSIZE);

%Setting up randomized array for a pracblock===============================
locGenHierarchy;
locSetUpEx;
%YAY!

%Main Loop Structure=======================================================
for BLOCK=[1:EXC_ST.NUMTBLOCK]
    
    %Message at the begining every block
    PREX.GO=0;PREX.REP=1;PREX.NUMTRIAL=[EXC_ST.NUMPTRIAL,EXC_ST.NUMETRIAL];
    PREX.TYPE=PARAM_ST.PRACTICE(BLOCK)+1;
    locShowMessege('BINTRO');
        
    %Timing for TaskTS and MusicTS
    BLOCKSTS=GetSecs;
    TKeeper=GetSecs;
    CCTime=GetSecs-TKeeper;
    locCodeEvent('TEMP',10+PREX.TYPE);%begining of block @@@@@@@@@@@@@@@@@@
    RECORD_ST.EV_BLOCK(BLOCK,:)=BLOCKSTS;%record timing for re-epoch!
    
    %Repeatedly loop trials for practice, and one loop for experiment!!
    while ~PREX.GO
        %Display rule first
        locDrawFrame('ON');
        locDrawRule('BRULE');
        Screen('Flip',W);KbWait;
        
        for TRIAL=1:PREX.NUMTRIAL(PREX.TYPE)
            %Within trial computations
            RECORD_ST.EV_TRIAL(BLOCK,TRIAL)=GetSecs-BLOCKSTS;
            [CCTime,TKeeper] = locDoTrial(CCTime,TKeeper);
            
            %Recording varibales of interest
            locKeepRecord(RECORD_ST,'VARS');
            
            %Exit out if esc is pressed
            [~,~,keyCode]=KbCheck;
            if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;
        end;
        
        %Clean up RECORD_ST??
        RECORD_ST.EV_END(BLOCK,:)=GetSecs-BLOCKSTS;
        
        %Block result (NEED THIS???)
        if BLOCK < EXC_ST.NUMTBLOCK,locShowMessege('BRESULT');end;
        
        %End Screen
        if BLOCK == EXC_ST.NUMTBLOCK,locShowMessege('EXEND');break;end;
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
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TASK_ST COLOR_ST KEY_ST W
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
    numPB=EXC_ST.NUMPBLOCK;numPT=EXC_ST.NUMPTRIAL;
    numEB=EXC_ST.NUMEBLOCK;numET=EXC_ST.NUMETRIAL;
    numTB=EXC_ST.NUMTBLOCK;
    
    switch fn
        case {'SEQSENCE'}
            %Practice:0=practice 1=experiment
            PARAM_ST.('PRACTICE')=BalanceTrials(numTB,0, 0:1);
            
            %Sequence Type (36 combination of bases ignoring order): 
            seqtype=BalanceTrials(numEB,1,1:SEQ_ST.NUMSEQTYPE);
            PARAM_ST.('SEQTYPE')=reshape(repmat(seqtype,1,2)',1,[])';
            
            %Cutoff Type:
            PARAM_ST.('CUTTIME')=repmat(TASK_ST.TIMECUT,[1,numTB/2])';
  
        case {'CONTENTS'}

            %Define parameters below for later analysis....!
            zeroD=zeros(numTB,max(numET));
            PARAM_ST.CHUNK=zeroD;%CHUNK = identity of chunk
            PARAM_ST.ACPOS=zeroD;%ACPOS = position across-chunks
            PARAM_ST.ELEMENT=zeroD;%ELEMENT = task (vertical,horizontal,diagonal)
            PARAM_ST.WCPOS=zeroD;%WCPOS = position within-chunk
            numC=SEQ_ST.NUMCHUNKS;
            chp=SEQ_ST.CHUNKPAIRS;
            flatL=@(d)reshape(d,1,[]);
    
            %Loop through each block
            for b=1:numTB
                %Change where to read for numtrial,chunk codde etc!
                seqC=chp(PARAM_ST.SEQTYPE(b),:);

                %CHUNK
                chunkC=flatL(repmat(seqC,SEQ_ST.MAXLENGTH/numC,1));
                PARAM_ST.CHUNK(b,:)=repmat(chunkC,[1,numET/length(chunkC)]);

                %ACPOS
                acposC=flatL(repmat(1:numC,SEQ_ST.MAXLENGTH/numC,1));
                PARAM_ST.ACPOS(b,:)=repmat(acposC,[1,numET/length(acposC)]);
                
                %ELEMENT
                elmC=SEQ_ST.('BASESEQUENCE')(PARAM_ST.SEQTYPE(b),:);
                elmC=horzcat(elmC{:});
                PARAM_ST.ELEMENT(b,:)=repmat(elmC,[1,numET/length(elmC)]);
                
                %WCPOS
                wcposC=1:SEQ_ST.MAXLENGTH/numC;
                PARAM_ST.WCPOS(b,:)=repmat(wcposC,[1,numET/max(wcposC)]);
                
                %BSPOS
                bsposC=1:SEQ_ST.MAXLENGTH;
                PARAM_ST.BSPOS(b,:)=repmat(bsposC,[1,numET/max(bsposC)]);
            end

        case {'TASK'}
            %STIMPOS:1=top_l,2=top_r,3-bot_r,4=bot_l clock-wise from top-l
            PARAM_ST.('STIMPOS')=RandSample(1:length(STIM_ST.OFFSETS),[numTB,max(numET)]);
            
            %DESTINATION:1=top_l,2=top_r,3-bot_r,4=bot_l clock-wise from top-l
            dref=TASK_ST.DESTREF;drefB=TASK_ST.DESTREF{end};
            elm=PARAM_ST.ELEMENT;elm(elm==0)=5;%temporally index!
            refK=[KEY_ST.T_LKEY,KEY_ST.T_RKEY;KEY_ST.B_LKEY,KEY_ST.B_RKEY];
            
            spref=cellfun(@(s)drefB==s,num2cell(PARAM_ST.('STIMPOS')),'Uni',false);
            PARAM_ST.('DESTINATION')=cell2mat(cellfun(@(s,e)dref{e}(s),spref,num2cell(elm),'Uni',false));
            desref=cellfun(@(s)drefB==s,num2cell(PARAM_ST.('DESTINATION')),'Uni',false);
            
             %Code keyboard responses
            PARAM_ST.('CORRESP')=cell2mat(cellfun(@(s)refK(s),desref,'Uni',false));
            
            %Put back unnecessary information (this exists for practice)
            PARAM_ST.('STIMPOS')(elm==5)=0;
            PARAM_ST.('DESTINATION')(elm==5)=0;
            PARAM_ST.('CORRESP')(elm==5)=0;      
    end
end

%RECORD SETTING============================================================
% Specify the order of data!
cellD=cell(EXC_ST.NUMTBLOCK,max(EXC_ST.NUMETRIAL));
zeroD=zeros(EXC_ST.NUMTBLOCK,max(EXC_ST.NUMETRIAL));
cVars={};
zVars={'SUBID','BLOCK','TRIAL','REP','RT','RESP','ACC','TASKTS','RSIJIT',...
    'EV_BLOCK','EV_TRIAL','EV_PREP','EV_STIM','EV_RESP','EV_END','MONEY'};%skip SUBID
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
global PARAM_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL BTRIAL BLOCKSTS
%**************************************************************************

%PHASE0:Trial Preparation--------------------------------------------------
while KbCheck,end;[~,~, keyCode]=KbCheck;
if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;

%Cover probe
Screen('FillRect', W,COLOR_ST.BACKG);
if RECORD_ST.ACC(BLOCK,max(TRIAL-1,1))==0,locDrawRule('RELEARN');end
locDrawFrame('ON');
Screen('Flip',W);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PHASE1:Random RSI+ Jittering PREPARATION----------------------------------
locCodeEvent('TEMP',9);%begining of trial @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
locCodeEvent('CHAR',[]);

% Random RSI period ()
% timeF=TIME_ST.RSIT;rsiJit=Sample((timeF*0.75):0.01:(timeF*1.25));%Random
timeF=TIME_ST.RSIT;rsiJit=timeF;%Fixed
if TRIAL==1,rsiJit=rsiJit*2;end;% will be removed anyway!
WaitSecs(rsiJit);RECORD_ST.RSIJIT(BLOCK,TRIAL)=timeF*1000;
   
% Fixed Preparation period
startT=GetSecs;
if RECORD_ST.ACC(BLOCK,max(TRIAL-1,1))==0,locDrawRule('RELEARN');end
locDrawFrame('ON');
Screen('Flip',W);
locCodeEvent('TEMP',1);%preparation onset@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
RECORD_ST.EV_PREP(BLOCK,TRIAL)=GetSecs-BLOCKSTS;%record timing for re-epoch!
while GetSecs-startT<TIME_ST.PREPT,end;%wait time

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PHASE2:STIMULUS DISPLAY & RECORDING---------------------------------------

%Drawing stimulus
%Measureing response!
if RECORD_ST.ACC(BLOCK,max(TRIAL-1,1))==0,locDrawRule('RELEARN');end
locDrawFrame('ON');
locDrawStim;
Screen('Flip',W);
locCodeEvent('TEMP',2);%stimulus onset@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
RECORD_ST.EV_STIM(BLOCK,TRIAL)=GetSecs-BLOCKSTS;%record timing for re-epoch!
locGetResponse;
locCodeEvent('RESP',[]);%resp onset@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
RECORD_ST.EV_RESP(BLOCK,TRIAL)=GetSecs-BLOCKSTS;%record timing for re-epoch!

%Showing feedback
if RECORD_ST.ACC(BLOCK,TRIAL)==0
    Screen('FillRect',W,COLOR_ST.RED);Screen('Flip', W);WaitSecs(0.1);
    Screen('FillRect',W,COLOR_ST.BACKG);
    locDrawFrame('ON');
    locDrawStim;
    locDrawRule('RELEARN');
    Screen('Flip',W);WaitSecs(0.75);
    
    %Force to type in correct answer to move on!
    coResp = PARAM_ST.CORRESP(BLOCK,TRIAL);
    while 1,[~,~,keyC]=KbCheck;
        if coResp==find(keyC),break;end;
        if keyC(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;
    end
end

%Event Code :End of trial and Characteristic codes@@@@@@@@@@@@@@@@@@@@@@@@@
%--------------------------------------------------------------------------


end

%##########################################################################

function locCodeEvent(codeT,code)
%==========================================================================
% locCodeEvent creates three different kind of event codes
% temporal event codes:codes for individual events regardress of contents
% characteristic codes:codes for the contents of one trial(all conditions)
%==========================================================================
%STIMULUS PROPERTY / TIME REGULATION / COLOR
global BLOCK TRIAL BTRIAL RECORD_ST
global EXC_ST SEQ_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST COLOR_ST KEY_ST
%**************************************************************************

switch codeT
    case {'TEMP'}%Temporal codes
        portCode=FPNAME_ST.PORTID_S{:};eventCode={code};
        
    case {'CHAR'}%Characteristic codes
        portCode=FPNAME_ST.PORTID_S{:};
        
        %-Block (100~): (block+100)
        eventCode(1,:)={BLOCK+100};
        
        %-Trial (200~): (trial+100)
        eventCode(2,:)={TRIAL+200};
        
    case {'RESP'}%Response codes
        portCode=FPNAME_ST.PORTID_R{:};
        
        %-Response: 
        eventCode(1,:)={RECORD_ST.RESP(BLOCK,TRIAL)};
        
        %-Accuracy:11 (correct),10 (wrong)
        eventCode(2,:)={RECORD_ST.ACC(BLOCK,TRIAL)+10};
        
        %-Lag 1 Accuracy:21 (correct),20 (wrong)
        eventCode(3,:)={RECORD_ST.ACC(BLOCK,max(TRIAL-1,1))+20};
        
end

%Actually write event code here!
if EXC_ST.EXMODE==1,cd(FPNAME_ST.PORDIR{:});end
for i=1:size(eventCode,1)
    switch EXC_ST.EXMODE
        case {'1'},write_parallel(portCode,eventCode{i,:});
        otherwise,disp(strcat(portCode,'_',num2str(eventCode{i})));
    end
end
if EXC_ST.EXMODE==1,cd('../');end;

end


%##########################################################################

function locGetResponse
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL PREX BLOCKSTS
%**************************************************************************

%Start measuring time------------------------------------------------------
startT = GetSecs;
coResp = PARAM_ST.CORRESP(BLOCK,TRIAL);
cutT = TIME_ST.CUTT(PARAM_ST.CUTTIME(BLOCK));

% Loop while there is time (Loop for one trial)
while (GetSecs - startT < cutT)
    [keyP, secs, keyC] = KbCheck;
    
    if keyP==1 && length(find(keyC))==1,break;end
    if keyC(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar(0);break;end;
end

%Record some variables!!
RECORD_ST.BLOCK(BLOCK,TRIAL)=BLOCK;
RECORD_ST.TRIAL(BLOCK,TRIAL)=TRIAL;
RECORD_ST.REP(BLOCK,TRIAL)=PREX.REP;
RECORD_ST.RT(BLOCK,TRIAL)=(secs-startT)*1000;
RECORD_ST.TASKTS(BLOCK,TRIAL)=(secs-BLOCKSTS)*1000;
if isempty(find(keyC, 1)),keyC=99;end;
RECORD_ST.RESP(BLOCK,TRIAL)=find(keyC);
if coResp==find(keyC),acc=1; else acc=0;end
RECORD_ST.ACC(BLOCK,TRIAL)=acc;

end


%##########################################################################

function locDrawRule(workT)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST SEQ_ST FPNAME_ST STIM_ST  TIME_ST COLOR_ST TASK_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST RECORD_ST BLOCK TRIAL PREX
%**************************************************************************

%Code Task Element in Sequence Advance! Keep position info, too.-----------
switch workT,case {'BRULE'},posSig=false;case {'RELEARN'},posSig=true;end
seqL=SEQ_ST.MAXLENGTH;
elmT=PARAM_ST.ELEMENT(BLOCK,1:seqL);
seqRule=TASK_ST.TASK(elmT);

%Drawing
for i=1:length(seqRule)
    %Change color for highligting!
    color=COLOR_ST.WHITE;
    if posSig && isequal(i,mod(TRIAL-1,seqL)+1),color=COLOR_ST.RED;end

    %Adjust Y cordinate to creat a chunk structure!
    if i<=3,yp=STIM_ST.YC-350+(30*(i-1));else yp=STIM_ST.YC-325+(30*(i-1));end;
    DrawFormattedText(W,seqRule{i},'center',yp,color);
end

end

%##########################################################################

function locDrawFrame(workT)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL IMAGE_ST
%**************************************************************************
% Frame
offHolder=CenterRectOnPoint(STIM_ST.FRAMERECT,STIM_ST.XC,STIM_ST.YC);
Screen('FrameRect',W,COLOR_ST.WHITE,offHolder,STIM_ST.FRAMESIZE);

% Fixation cross
switch workT,case{'OFF'},c='BLACK';case{'ON'},c='WHITE';end
fixation='+';Screen('TextSize',W,STIM_ST.FONTSIZE+25);
DrawFormattedText(W,fixation,'center',STIM_ST.YC+STIM_ST.FXADJ,COLOR_ST.(c));
Screen('TextSize',W,STIM_ST.FONTSIZE);%put back the textsize

end

%##########################################################################

function locDrawStim
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL IMAGE_ST
%**************************************************************************

% Setting up XY cordinates
locC=PARAM_ST.STIMPOS(BLOCK,TRIAL);
xLoc=STIM_ST.OFFSETS(1,locC);
yLoc=STIM_ST.OFFSETS(2,locC);

% Drawing white square and cover by black one
% to create a directional point!
stim=CenterRectOnPoint(STIM_ST.STIMRECT,STIM_ST.XC+xLoc,STIM_ST.YC+yLoc);
Screen('FillOval',W,COLOR_ST.WHITE,stim);

% Debug
% DrawFormattedText(W,sprintf('stimpos:%d',PARAM_ST.STIMPOS(BLOCK,TRIAL)),'center',STIM_ST.YC+200,COLOR_ST.WHITE);
% DrawFormattedText(W,sprintf('element:%d',PARAM_ST.ELEMENT(BLOCK,TRIAL)),'center',STIM_ST.YC+250,COLOR_ST.WHITE);
% DrawFormattedText(W,sprintf('resp:%d',PARAM_ST.CORRESP(BLOCK,TRIAL)),'center',STIM_ST.YC+300,COLOR_ST.WHITE);


end

%##########################################################################

function locShowMessege(msgtime)
%==========================================================================
% locShowMessege displays all of string based messages for a variety of
% purposes depending on msgtime in combination with locCondMessage
%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST SEQ_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL PREX
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
        TEXT(4)={'Be accurate and fast as much as you can'};
        TEXT(5)={'Press any key to proceed'};
        C={'RED','WHITE','WHITE','WHITE','WHITE'};

    case {'BRESULT'}% Ending Block BREAK screen (no money)
        [ConText]=locCondMessage(msgtime);
        TEXT(1)={strcat('You have completed block ',num2str(BLOCK),' of ',num2str(EXC_ST.NUMTBLOCK))};
        TEXT(2)={ConText{1}};
        TEXT(3)={ConText{2}};
        TEXT(4)={ConText{3}};
        TEXT(5)={'Press any key to proceed'};
        C=repmat({'WHITE'},1,size(TEXT,2));

        
    case {'EXEND'}% Ending Experiment screen
        TEXT(1)={'Thank you'};
        TEXT(2)={'You have completed the experiment'};
        TEXT(3)={''};
        TEXT(4)={''};
        TEXT(5)={'Press any key to proceed'};
        C=repmat({'WHITE'},1,size(TEXT,2));

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
while KbCheck; end;KbWait;WaitSecs(1.5);

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
global EXC_ST TASK_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL PREX
%**************************************************************************

switch msgtime
    case {'BINTRO'}% Begining Block BREAK screen
        expCond={'PRACTICE','EXPERIMENT'};
        ConText(1)={sprintf('ATTENTION!! This is an %s!!',expCond{PREX.TYPE})};
        ConText(2)={'Please remember a sequence of rules'};

    case {'BRESULT'}% Begining Block BREAK screen
        %Change Setting for Practice and Experiment
        trial=PREX.NUMTRIAL(PREX.TYPE);PREX.GO=true;
        errorScore=length(find(RECORD_ST.ACC(BLOCK,1:trial)==0));
        errorCheck=(errorScore/trial*100)<=TASK_ST.ERRORCUT(PREX.TYPE);
        ConText(1)={sprintf('You missed %d out of %d trials',errorScore,trial)};
        ConText(2)={sprintf('Please go on to the next block!')};

        %Evaluation for practice blocks
        if PREX.TYPE==1 && ~errorCheck
            PREX.GO=false;PREX.REP=PREX.REP+1;
            ConText(2)={'Too many mistakes!! Please repeat practice!!'};
            
            %Re-randomize Stimpos & correct answers
            stIDX=RandSample(1:EXC_ST.NUMPTRIAL:EXC_ST.NUMETRIAL);
            PARAM_ST.STIMPOS(BLOCK,:)=circshift(PARAM_ST.STIMPOS(BLOCK,:),[0,-stIDX+1]);
            PARAM_ST.DESTINATION(BLOCK,:)=circshift(PARAM_ST.DESTINATION(BLOCK,:),[0,-stIDX+1]);
            PARAM_ST.CORRESP(BLOCK,:)=circshift(PARAM_ST.CORRESP(BLOCK,:),[0,-stIDX+1]);
            
            %Erase data from this attempt
            RECORD_ST.RT(BLOCK,:)=0;RECORD_ST.ACC(BLOCK,:)=NaN;RECORD_ST.RESP(BLOCK,:)=0;
        end
        
        %Accuracy and RT performance check
        prevRT=RECORD_ST.RT(1:max(BLOCK-1,1),:);
        prevRT=prevRT(logical(PARAM_ST.PRACTICE(1:size(prevRT,1))),:);
        if isempty(prevRT),prevRT=RECORD_ST.RT(1,:);end;%very first block!
        useIDX=prevRT>0;prevRT=prevRT(useIDX);
        cutRT=prctile(prevRT,TASK_ST.RTCUT);
        fastRTs=find(RECORD_ST.RT(BLOCK,1:TRIAL)<cutRT);   
        
        %if accuarcy is okay, give money for all fast RTs
        if all(errorCheck)&&PREX.TYPE==2,moneyE=length(fastRTs);else moneyE=0;end
        RECORD_ST.MONEY(BLOCK,:)=moneyE;        
        moneyT=sum(RECORD_ST.MONEY(1:BLOCK,1));
        ConText(3)={sprintf('You earned %d cents, and the total is %d cents.',moneyE,moneyT)};
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
    if size(var,1)==max(EXC_ST.NUMTBLOCK),varType='BCOND';ref=[BLOCK,1];end
    if size(var,2)==max(EXC_ST.NUMETRIAL),varType='TCOND';ref=[BLOCK,TRIAL];end;
    
    %Put Data First
    if iscell(var),value=var{ref(1)};else value=var(ref(1),ref(2));end;ActArray(g,1)={value};
    
    switch ischar(ActArray{g,1});
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

function [path]=fileroot(fileName)
path=which(fileName); %path is to file
i=find(filesep==path); %file seperator
path=path(1:i(end));
end
