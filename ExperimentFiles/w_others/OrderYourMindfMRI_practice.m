function OrderYourChunks_practice
%==========================================================================
% PURPOSE:
% 1,By using forward encoding model, try to extract a neural signature of
% position in the task
% TO DO:
% 1,

% SETTING:=================================================================
% BASELINE = 1000 ms
% CUTOFF = inf ms (used to be 750 ms)
% SEQUENCE = 6 sub-sequences
% SELF PACING = OFF
%==========================================================================
clear all;
close all;
sca;
if 1, Screen('Preference', 'SkipSyncTests', 1); end
commandwindow;while KbCheck; end

%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST TASK_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST PARAM_ALL
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL PRACEXC REPCOUNT BTRIAL BLOCKSTS
%**************************************************************************

%=================== Setting Standard variables============================

%TASK SETTING
TASK_ST.LAG2SEQENCES='on';
TASK_ST.BASECHUNK(1,:)=[1,2,3];
TASK_ST.BASECHUNK(2,:)=[1,3,2];
TASK_ST.BASECHUNK(3,:)=[2,1,3];
TASK_ST.BASECHUNK(4,:)=[2,3,1];
TASK_ST.BASECHUNK(5,:)=[3,1,2];
TASK_ST.BASECHUNK(6,:)=[3,2,1];
TASK_ST.NUMCHUNK=size(TASK_ST.BASECHUNK,1);
switch TASK_ST.LAG2SEQENCES
    case {'on'},%Ignoring chunk associations (12 sequences)
        chunkpairs(1,:)=[1,2];chunkpairs(2,:)=[1,3];
        chunkpairs(3,:)=[2,1];chunkpairs(4,:)=[2,5];
        chunkpairs(5,:)=[3,4];chunkpairs(6,:)=[3,1];
        chunkpairs(7,:)=[4,3];chunkpairs(8,:)=[4,6];
        chunkpairs(9,:)=[5,6];chunkpairs(10,:)=[5,2];
        chunkpairs(11,:)=[6,5];chunkpairs(12,:)=[6,4];
    case {'off'},%Full counterbalancing (30 sequences)
        chunkpairs=allcomb([1:TASK_ST.NUMCHUNK],[1:TASK_ST.NUMCHUNK]);
        % 1,when order matters->allcomb! 2,when order does not matter->combnk!
        chunkpairs=chunkpairs(find(chunkpairs(:,1)-chunkpairs(:,2)~=0),:);
end
TASK_ST.CHUNKPAIRS=chunkpairs;
TASK_ST.NUMSEQTYPE=size(chunkpairs);
TASK_ST.SERTYPE=[1:TASK_ST.NUMSEQTYPE];%sequence type
TASK_ST.MAXSEQLENGTH=[3];%smaximum elements in sequence
TASK_ST.ERRORCUT=1;%error percentage cut off


%EXPERIMENT SETTING
EXC_ST.HIDEFB='off';
EXC_ST.NUMPBLOCK=size(TASK_ST.BASECHUNK,1);
EXC_ST.NUMPTRIAL=[];
EXC_ST.NUMREP=1;
PARAM_ALL=struct;

% Set up the timer
TIME_ST.YEILDT=0.00001;
TIME_ST.TRIALT_P=inf;
TIME_ST.CUTBASE =inf;
TIME_ST.RCIT=0.750;
TIME_ST.BASELINE=1.000;
TIME_ST.CHOICET=.750;

%FILE PATH SETTING
FPNAME_ST.TXTCODE='ordch_prac';
FPNAME_ST.EXFILE='OrderYourChunks_practice.m';
FPNAME_ST.PTLPATH=PsychtoolboxRoot;
exFilePath=which(FPNAME_ST.EXFILE);
FPNAME_ST.EXFPATH=exFilePath(1:end-length(FPNAME_ST.EXFILE));
cd(FPNAME_ST.EXFPATH);
FPNAME_ST.PORTID_S={'DCC8'};
FPNAME_ST.PORTID_R={'DCD8'};
FPNAME_ST.PORDIR={strcat(pwd,'\PortCode_Supports')};


%STIMULUS SETTING++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%General
%General
STIM_ST.RESOLUTION=[1024,768];
STIM_ST.FONTSIZE=22;
STIM_ST.RADSIZE=120;

%Basic obejects info
STIM_ST.YCOFF=-65;
STIM_ST.ADJUSTX=[0,-50,-120,-190,-260];
STIM_ST.FRAMESIZE=5;
STIM_ST.SQSIZE=80;
STIM_ST.PLSIZE=85;
STIM_ST.STIMRECT=ceil([0 0 STIM_ST.SQSIZE*1.43 STIM_ST.SQSIZE]);
STIM_ST.PLACERECT=ceil([0 0 STIM_ST.PLSIZE*1.43 STIM_ST.PLSIZE]);


%LocationHolders
% NOTE:
% 1, window size should be odd number!!
% 2, reference does not include "a center window"
STIM_ST.WINDMAP=cell(2,2);
STIM_ST.WINDMIDDLE=5;
allwind=1:size(STIM_ST.WINDMAP,1)*size(STIM_ST.WINDMAP,2);
STIM_ST.WINDOPEN=setdiff(allwind,STIM_ST.WINDMIDDLE);
STIM_ST.BOXRANGE=[0.12];%up to -BOXRANGE to BOXRANGE
STIM_ST.JITRANGE=[0.00];

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
COLOR_ST.RULE(1,:)=COLOR_ST.RED;
COLOR_ST.RULE(2,:)=COLOR_ST.GREEN;
COLOR_ST.RULE(3,:)=COLOR_ST.PURPLE;
COLOR_ST.RULE(4,:)=COLOR_ST.BLUE;
COLOR_ST.RULE(5,:)=COLOR_ST.YELLOW;
COLOR_ST.RULE(6,:)=COLOR_ST.ORANGE;

% Set keys
KbName('UnifyKeyNames');
KEY_ST.T_LEFTKEY = KbName('4');%4
KEY_ST.T_RIGHTKEY = KbName('5');%5
KEY_ST.B_LEFTKEY=KbName('1');%1
KEY_ST.B_RIGHTKEY=KbName('2');%2
KEY_ST.YKEY = KbName('y');
KEY_ST.NKEY = KbName('n');
KEY_ST.SPACEKEY = KbName('SPACE');
KEY_ST.ESCKEY = KbName('ESCAPE');
KEY_ST.SAMEKEY = KbName('UpArrow');
KEY_ST.DIFFKEY = KbName('DownArrow');


try
    %Get subject ID and open data file for output==========================
    EXC_ST.SUBID = input('Enter Subject ID: ','s');
    EXC_ST.EXMODE = input('EEG(1) / Behavioral(0): ','s');
    path=fileroot(FPNAME_ST.EXFILE);
    RECORD_ST.DFILE = fopen([path [EXC_ST.SUBID FPNAME_ST.TXTCODE '.txt']],'at');
    %cd(strcat(path,'PortCode_Supports'));
    if exist('rng','file'),rng(str2double(EXC_ST.SUBID),'twister');else rand('twister',str2double(EXC_ST.SUBID));end
    
    
    %Screen Stuff==========================================================
    % Open a graphics window on the main screen
    % using the PsychToolbox's Screen function.
    STIM_ST.SCCODE=max(Screen('Screens'));
    Screen('Resolution',STIM_ST.SCCODE,STIM_ST.RESOLUTION(1),STIM_ST.RESOLUTION(2),[],32);
    [W, wRect]=Screen('OpenWindow',STIM_ST.SCCODE, 0,[],32,2);%[0 0 200 200]
    Screen('FillRect',W,COLOR_ST.BACKG);
    
    %Prevent Typing in the code!
    ListenChar(2);
    HideCursor;
    
    %MXY = matrix of xy coordinates(6 pairs of 12 points)
    %finds screen dimensions, returns center x,y
    [swidth, sheight]=Screen('WindowSize',W);
    STIM_ST.XC=fix(swidth/2);
    STIM_ST.YC=fix(sheight/2);
    
    Screen('TextFont',W,'Helvetica');
    Screen('TextSize',W,STIM_ST.FONTSIZE);
    
    %Instruction===========================================================
    %locShowInstruction;
    %READY?
    
    %Setting up randomized array for a pracblock===========================
    locLoadStim;
    locRadialCSN;
    locGenHierarchy;
    TASK_ST.PRACORDER=shuffle(TASK_ST.BASECHUNK,1);
    %YAY!
    
    
    %Main Loop Structure===================================================
    for BLOCK=[1:EXC_ST.NUMPBLOCK]
        
        %Message at the begining every block
        PRACEXC=false;GOFLAG=0;REPCOUNT=0;
        locShowMessege('BINTRO');
        
        %Display rule first
        locDrawFix('NEUTRAL');
        locDrawRule('BRULE');
        Screen('Flip',W);KbWait;
        
        %Timing for TaskTS and MusicTS
        BLOCKSTS=GetSecs;
        TKeeper=GetSecs;
        CCTime=GetSecs-TKeeper;
        
        %Repeatedly loop trials for practice, and one loop for experiment!!
        while ~GOFLAG 
            %Update PARAM_ST
            locSetUpBlock;

            for TRIAL=1:EXC_ST.NUMPTRIAL
                %Within trial computations
                BTRIAL=((BLOCK-1)*EXC_ST.NUMPTRIAL)+TRIAL;
                [CCTime,TKeeper] = locDoTrial(CCTime,TKeeper);
                
                %Recording varibales of interest
                locKeepRecordM(RECORD_ST,'VARS');
                
                %Exit out if esc is pressed
                [~,~,keyCode]=KbCheck;
                if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;
            end;
            
            %Block result (NEED THIS???)
            if BLOCK < EXC_ST.NUMPBLOCK,locShowMessege('BRESULT');end;
            
            %End block Evaluation
            %Evaluate if subject is ready to move on after practice!
            accScore=length(find(RECORD_ST.ACC(:)==0));
            if ~PRACEXC && (accScore/EXC_ST.NUMPTRIAL*100)<TASK_ST.ERRORCUT,GOFLAG=1;end;%practice
            REPCOUNT=REPCOUNT+1;
        end
        
        %End Screen
        if BLOCK == EXC_ST.NUMPBLOCK,locShowMessege('EXEND');end;
    end
    
catch ME
    ME.message
    rethrow(ME);
    fprintf('error\n');
    ShowCursor;
    ListenChar;
    Screen('CloseAll');
end

%Putting everything back
Screen('CloseAll');
ShowCursor;
ListenChar;


end


%##########################################################################


function y = locSetUpBlock
%locSetUpBlock

%**************************************************************************
%BLOCK ORGANIZATION/%TIME REGULATION/%STIMULUS/%KEYS/%FONT/%COLOR/%PARAMETER
global EXC_ST TASK_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST PARAM_ALL BLOCK TRIAL REPCOUNT
%**************************************************************************

% NOTE:
% 1,'WUPDATE' needs to be randomly 1 or 0
% color, orientation, xy cordinates need to be kept same as long as it's 0
% think some clever way to do this without for loop

% 2,'INTERRUPT' needs to be randomly 1 or 0 again
% if it is 1, interpose interruption at the begining of trial....
% exclude first trial????


%PARAMST SETTING===========================================================
PARAM_ST.PARAMETER={
    'SEQELM_POSITION','ManV';
    'RESPCODE','ManV';
    'RESPLOC','ManV';};

%RANDOMIZATION LOOP========================================================
for p=1:length(PARAM_ST.PARAMETER)
    %Process1:Parameter Setting!+++++++++++++++
    fn=PARAM_ST.PARAMETER{p,1};
    system=PARAM_ST.PARAMETER{p,2};
    varName=fn(1:strfind(fn,'Code')-1);
    
    %Process2:System Setting!++++++++++++++++++
    switch system
        case {'LBal'},sysLBalance(fn);%Loose balance
        case {'RBal_s'},sysRBalance(fn);%Rigid balance
        case {'Alt'},sysAlternate(fn);%Alternate
        case {'Prob'},sysProbability(fn);
        case {'RanF'},sysRandSample(fn);
        case {'CouF'},sysCountFactor(fn);
        case {'ManV'},sysManualV(fn);
        case {'Shuf'},sysShuffle(fn); %Not Vectorized!
        case {'VecC'},sysVecCopy(fn);
        case {'Jitt'},sysRandJitter(fn);
        case {'GenEq'},sysGenEquation(fn);
        case {'Rept'},sysGenRepeats(fn);
        otherwise,
    end
    %PARAMST
end

%RECORD SETTING============================================================
% Add BLOCK and TRIAL to PARAM_ST to make it easier to Join!!
EXC_ST.NUMPTRIAL=length(PARAM_ST.SEQELM)
PARAM_ST.BLOCK=repmat([BLOCK]',1,EXC_ST.NUMPTRIAL);
PARAM_ST.TRIAL=1:EXC_ST.NUMPTRIAL;

% Specify the order of data!
RECORD_ST.PARAMORDER={'SUBID','REPCOUNT','BLOCK','TRIAL',...
    'IT','INITIALTS','RT','RESP','ACC','TASKTS','RSIJIT'};

% Add place holders to all recording variables!
for i=1:length(RECORD_ST.PARAMORDER)
    pName=RECORD_ST.PARAMORDER{i};%renew RECORD_ST
    RECORD_ST.(pName)=zeros(1,EXC_ST.NUMPTRIAL);
end

% Write header
RECORD_ST.SUBID=str2double(EXC_ST.SUBID);
if BLOCK==1,locKeepRecordM(RECORD_ST ,'header');end

%RECORD PARAMST JUST IN CASE!!=============================================
PARAM_ALL.(strcat('b',num2str(BLOCK),'r',num2str(REPCOUNT)))=PARAM_ST;
paramtMat=[EXC_ST.SUBID '_' FPNAME_ST.TXTCODE '_PARAM_ST.mat'];
save(paramtMat,'PARAM_ALL');

disp('YAY! Completed!!');

end

%##########################################################################

function sysManualV(fn)
%Manual system
%channelT=dependency type
%Note:
% 1,This function is not flexible
% 2,It is basically do whatever computation! BUT, they should all depend on
% other factors that were computed earlier....Otherwise,use other funcs!
%==========================================================================
global EXC_ST TASK_ST PARAM_ST BLOCK STIM_ST TIME_ST COLOR_ST KEY_ST
%**************************************************************************

%Parameter Sign
switch fn
    case {'SEQELM_POSITION'}
        %randomize chunks and duplicate!
        seq=num2cell(TASK_ST.PRACORDER(1:BLOCK,:),2);
        seq=shuffle(reshape(repmat(seq,1,EXC_ST.NUMREP),1,[]));
        chunkcode=cellfun(@(s) find(ismember(TASK_ST.BASECHUNK,s,'rows')),seq,'UniformOutput',false);
        
        %Sequence & Chunk Level
        csize=size(TASK_ST.BASECHUNK,2);
        PARAM_ST.SEQELM=horzcat(seq{:});
        PARAM_ST.POSITION_CH=repmat(1:csize,1,length(PARAM_ST.SEQELM)/csize);
        PARAM_ST.CHUNKCODE=cell2mat(reshape(repmat(chunkcode,csize,1),1,[]));
        
    case {'RESPCODE'},
        respcode=repmat({TASK_ST.BASECHUNK(1,:)},1,length(PARAM_ST.SEQELM));
        PARAM_ST.(fn)=cellfun(@(se) Shuffle(se),respcode,'UniformOutput',false);
        
    case {'RESPLOC'},
        %key assingment
        respInd=cellfun(@(rc,elm) find(rc==elm),PARAM_ST.RESPCODE,num2cell(PARAM_ST.SEQELM),'UniformOutput',false);
        PARAM_ST.(fn)=cell2mat(respInd);
end

end

%##########################################################################

function [CCTime,TKeeper] = locDoTrial(CCTime,TKeeper)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
global PARAM_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL BTRIAL BLOCKSTS
%**************************************************************************


%PHASE1:Trial Preparation----------------------------------------------
while KbCheck,end;[~,~, keyCode]=KbCheck;
if keyCode(KEY_ST.ESCKEY),Screen('CloseAll');ListenChar;end;

%Cover probe
Screen('FillRect', W,COLOR_ST.BACKG);
Screen('Flip',W);

%Jittering RSI---------------------------------------------------------
timeF=TIME_ST.RCIT;
rsiJit=Sample((timeF*0.75):0.01:(timeF*1.25));
WaitSecs(rsiJit);RECORD_ST.RSIJIT(TRIAL)=rsiJit*1000;

%Self-pasing
locDrawFix('PRETRIAL');
Screen('Flip',W);
locGetInitiation;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PHASE1:FIXATION DISPLAY & FIXED RSI INTERVAL--------------------------
%DO NOT Pause until subject press space bar!(trial starts automatically)
locDrawFix('NEUTRAL');
Screen('Flip',W);
WaitSecs(TIME_ST.BASELINE);


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PHASE2:STIMULUS DISPLAY & RECORDING-----------------------------------

%Drawing stimulus
%Measureing response!
locDrawRFrame;
locDrawFix('NEUTRAL');
Screen('Flip',W);
[CCTime,TKeeper] = locGetResponse('TEST',CCTime,TKeeper);


%Showing feedback (show chosen elements regardless?)
locDrawChoice;
locDrawRFrame;
Screen('Flip',W);
WaitSecs(TIME_ST.CHOICET);

if RECORD_ST.ACC(TRIAL)==0,
    while KbCheck,end;
    locDrawRFrame;
    locDrawChoice;
    locDrawRule('RELEARN');
    Screen('Flip',W);WaitSecs(1);KbWait;
end


%Event Code :End of trial and Characteristic codes@@@@@@@@@@@@@@@@@@@@@
%----------------------------------------------------------------------
HideCursor;
% if BLOCK==2,sca;keyboard;end

end


%##########################################################################

function y =locGetInitiation
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST TASK_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST RECORD_ST
global AUDIO_ST COLOR_ST KEY_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL BTRIAL BLOCKSTS W
%**************************************************************************

%Start measuring time
startT = GetSecs;

% Loop while there is time (Loop for one trial)
while (GetSecs - startT < inf)
    [keyIsPressed, secs, keyCode] = KbCheck;
    if keyIsPressed==1 && length(find(keyCode))==1,
        switch find(keyCode)
            case {KEY_ST.SPACEKEY},break;
            case {KEY_ST.ESCKEY},Screen('CloseAll');ListenChar;break;
        end
    end;
end

if keyIsPressed
    %Record some variables!!
    RECORD_ST.IT(TRIAL)=(secs-startT)*1000;
    RECORD_ST.INITIALTS(TRIAL)=(secs-BLOCKSTS)*1000;
end

end

%##########################################################################

function [CCTime,TKeeper] = locGetResponse(workT,CCTime,TKeeper)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST TASK_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST RECORD_ST PRACEXC
global AUDIO_ST COLOR_ST KEY_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL REPCOUNT BLOCKSTS W
%**************************************************************************

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PHASE1:Trial Preparation--------------------------------------------------
[~,~,buttons] = GetMouse;
while any(buttons),[~,~,buttons] = GetMouse;end
while KbCheck,end;[~, ~, keyCode]=KbCheck;
ShowCursor;
confirmF=0;
captureF=0;
SetMouse(STIM_ST.XC,STIM_ST.YC-STIM_ST.YCOFF);

%Start measuring time------------------------------------------------------
startT = GetSecs;
coResp = PARAM_ST.RESPLOC(TRIAL);

% Loop and track the mouse, drawing the contour
while confirmF==0
    %Update keyboard and mouse button status
    [x,y,buttons] = GetMouse(W);
    
    %While any button is pressed
    if any(buttons) && captureF==0,
        %check if the mouse is in template
        secs=GetSecs;
        tempcheck=IsInTemplates([1:TASK_ST.MAXSEQLENGTH],x,y);
        chosenItem=find(tempcheck);
        if ~isempty(chosenItem),captureF=1;end
    end
    
    %After buttons a stimulus is chosen
    if captureF==1,
        if ~isempty(chosenItem),
            %Record some variables!!
            RECORD_ST.RT(TRIAL)=(secs-startT)*1000;
            RECORD_ST.TASKTS(TRIAL)=(secs-BLOCKSTS)*1000;
            RECORD_ST.RESP(TRIAL)=chosenItem;confirmF=1;
            if coResp==chosenItem
                RECORD_ST.ACC(TRIAL)=1;
            end
        end
    end
    if keyCode(KEY_ST.ESCKEY),
        Screen('CloseAll');ListenChar;
        if ispc,ShowHideWinTaskbarMex(1);end;
    end;
end

%Record other variables====================================================
if strcmp(workT,'TEST')
    RECORD_ST.REPCOUNT(TRIAL)=REPCOUNT;
    RECORD_ST.BLOCK(TRIAL)=BLOCK;
    RECORD_ST.TRIAL(TRIAL)=TRIAL;
end

end

%##########################################################################

function reply = IsInTemplates(range,x,y)
%**************************************************************************
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST TASK_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST RECORD_ST
global AUDIO_ST COLOR_ST KEY_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL BTRIAL BLOCKSTS W
%**************************************************************************

% initialize
reply=nan(1,size(range,1));

for t=range
    %#FLAT FORMAT: Change reference!!
    %xOff=STIM_ST.ADJUSTX(length(TASK_ST.BASEELM));
    %tx=STIM_ST.XC+xOff+((STIM_ST.PLSIZE+30)*(t-1));ty=STIM_ST.YC+180;
    
    %#RADIAL FORMAT: Change reference!!
    tx=STIM_ST.RCXYREF(1,t);ty=STIM_ST.RCXYREF(2,t);
    NImage=CenterRectOnPoint([0,0,STIM_ST.PLSIZE+25 STIM_ST.PLSIZE+25],tx,ty);
    reply(t)=IsInRect(x,y,NImage);
end

end


%##########################################################################

function locDrawFix(workT)
%==========================================================================
% locDrawFix draws a very nice fixation cross in the center of screen
%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
%STIMULUS PROPERTY / TIME REGULATION / DATA RECORDING /AUDIO CUE / COLORS / KEYS
global EXC_ST TASK_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST RECORD_ST
global AUDIO_ST COLOR_ST KEY_ST IMAGE_ST
%EXPERIMENT ONLINE GLOBALS
global BLOCK TRIAL BTRIAL BLOCKSTS W
%**************************************************************************

Screen('TextFont',W,'Courier');
switch workT,
    case {'PRETRIAL'},f='+';chunkcode=PARAM_ST.CHUNKCODE(TRIAL);color=COLOR_ST.RULE(chunkcode,:);fsize=STIM_ST.FONTSIZE+20;
    case{'NEUTRAL'},f='+';color=COLOR_ST.BLACK;fsize=STIM_ST.FONTSIZE+20;
    case{'WRONGCHOICE'},f='?';color=COLOR_ST.BLACK;fsize=STIM_ST.FONTSIZE+100;
end

Screen('TextSize',W,fsize);
DrawFormattedText(W,f,'center',STIM_ST.YC-STIM_ST.YCOFF-25,color);
Screen('TextSize',W,STIM_ST.FONTSIZE);%put back the textsize


end
%##########################################################################

function locDrawRule(workT)
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPNAME_ST STIM_ST  TIME_ST COLOR_ST IMAGE_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST RECORD_ST BLOCK TRIAL PRACEXC
%**************************************************************************

% =========================================================================
%Code Task Element in Sequence Advance! Keep position info, too.-----------
xOff=STIM_ST.ADJUSTX(length(TASK_ST.BASECHUNK(1,:)));maxseq=TASK_ST.MAXSEQLENGTH;
switch workT,case {'BRULE'},posSig=false;tadj=1;case {'RELEARN'},posSig=true;end

% Finally draw all task rules for a given sequence!!!YAY!!-----------------
for s=1:maxseq
    %Specify task and location!
    xp=STIM_ST.XC+xOff+((STIM_ST.PLSIZE+40)*(s-1));yp=STIM_ST.YC-200;
    
    %Color------------------------->
    switch posSig
        case {0},%initial presentation (show the new added chunk)
            text='Please remember this new chunk!!';
            seq=TASK_ST.PRACORDER(BLOCK,:);%get elements first
            newchunkcode=find(ismember(TASK_ST.BASECHUNK,seq,'rows'));
            color=COLOR_ST.RULE(newchunkcode,:);
        case {1},%re-learning
            text='Please remember this chunk again!!';
            chunkcode=PARAM_ST.CHUNKCODE(TRIAL);%get chunkcode first
            seq=TASK_ST.BASECHUNK(chunkcode,:);
            color=COLOR_ST.WHITE;if s==mod(TRIAL,maxseq),color=COLOR_ST.RULE(chunkcode,:);end;
            if s==maxseq && ~mod(TRIAL,maxseq),color=COLOR_ST.RULE(chunkcode,:);end;
    end
    
    %Task-------------------------->
    im=CenterRectOnPoint(STIM_ST.PLACERECT,xp,yp);
    Screen('DrawTexture',W,IMAGE_ST.IMARRAY{seq(s),3},[],im);
    %Frame------------------------->
    fm=CenterRectOnPoint([0,0,STIM_ST.PLSIZE+25 STIM_ST.PLSIZE+25],xp,yp);
    Screen('FrameRect', W,color,fm,STIM_ST.FRAMESIZE);
    %Screen('DrawText',W,num2str(s),xp-5,yp-80,COLOR_ST.WHITE);
end

% Put some instruction
DrawFormattedText(W,text,'center',STIM_ST.YC-300,COLOR_ST.BLACK);

end

%##########################################################################

function locDrawRFrame
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL IMAGE_ST
%**************************************************************************

% Update response coding for each elements every trials
rscode=PARAM_ST.RESPCODE{TRIAL};

for s=1:length(rscode)
    %xOff=STIM_ST.ADJUSTX(length(TASK_ST.BASEELM));
    %xp=STIM_ST.XC+xOff+((STIM_ST.PLSIZE+30)*(s-1));yp=STIM_ST.YC+180;
    xp=STIM_ST.RCXYREF(1,s);yp=STIM_ST.RCXYREF(2,s);
    
    %Element-------------------------->
    im=CenterRectOnPoint(STIM_ST.PLACERECT,xp,yp);
    Screen('DrawTexture',W,IMAGE_ST.IMARRAY{rscode(s),3},[],im);
    %Frame------------------------->
    fm=CenterRectOnPoint([0,0,STIM_ST.PLSIZE+25 STIM_ST.PLSIZE+25],xp,yp);
    %Screen('FrameOval', W,COLOR_ST.BLACK,fm,STIM_ST.FRAMESIZE);
    %Screen('DrawText',W,num2str(s),xp-5,yp-80,COLOR_ST.WHITE);
end

% Debug
% Screen('DrawText',W,num2str(PARAM_ST.ONTRACK(BLOCK,TRIAL)),STIM_ST.XC,STIM_ST.YC-250,COLOR_ST.WHITE);

end

%##########################################################################

function locDrawChoice
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL IMAGE_ST
%**************************************************************************

%Element-------------------------->
rcode=PARAM_ST.RESPCODE{1,TRIAL};
choice=rcode(RECORD_ST.RESP(TRIAL));
switch isempty(choice)
    case {true},%non-coed response!
        locDrawFix('WRONGCHOICE');
        
    case {false},%coded response!
        im=CenterRectOnPoint(STIM_ST.STIMRECT,STIM_ST.XC,STIM_ST.YC-STIM_ST.YCOFF);
        Screen('DrawTexture',W,IMAGE_ST.IMARRAY{choice,3},[],im);
end

%Frame------------------------->
% fm=CenterRectOnPoint([0,0,STIM_ST.PLSIZE+25 STIM_ST.PLSIZE+25],STIM_ST.XC,STIM_ST.YC);
% Screen('FrameRect', W,COLOR_ST.BLACK,fm,STIM_ST.FRAMESIZE);
%Screen('DrawText',W,num2str(s),xp-5,yp-80,COLOR_ST.WHITE);

end


%##########################################################################

function locShowMessege(msgtime)
%==========================================================================
% locShowMessege displays all of string based messages for a variety of
% purposes depending on msgtime in combination with locCondMessage
%==========================================================================
%BLOCK ORGANIZATION / RANDOMIZATION / DATA STORAGE / FILE PATH & NAMES
global EXC_ST TASK_ST FPNAME_ST STIM_ST RECORD_ST TIME_ST COLOR_ST KEY_ST W
%BLOCK ORGANIZATION
global PARAM_ST BLOCK TRIAL PRACEXC
%**************************************************************************

%Consistant settings
Screen('FillRect',W,COLOR_ST.BACKG);
Screen('TextSize',W,STIM_ST.FONTSIZE);

switch msgtime
    case {'BINTRO'},% Begining Block BREAK screen
        [ConText]=locCondMessage(msgtime);
        TEXT(1)={ConText{1}};
        TEXT(2)={ConText{2}};
        TEXT(3)={''};
        TEXT(4)={'Be accurate and fast as much as you can'};
        TEXT(5)={'Press any key to proceed'};
        
    case {'MONEY'},% At the end of practice blocks
        TEXT(1)={'From this block on you can earn MONEY!!'};
        TEXT(2)={'To do so you just have to get 90 correct'};
        TEXT(3)={'and beat your overall average response time'};
        TEXT(4)={'Each trial you beat your time we wll give you a penny!'};
        TEXT(5)={'Press any key to proceed'};
        
    case {'BHAND'},% Begining Block BREAK screen
        TEXT(1)={'ATTENTION: Please wait!'};
        TEXT(2)={'PLEASE SWITCH HAND to the other side!!'};
        TEXT(3)={''};
        TEXT(4)={''};
        TEXT(5)={'Press any key to proceed'};
        
    case {'BRESULT'},% Ending Block BREAK screen (no money)
        [ConText]=locCondMessage(msgtime);
        TEXT(1)={strcat('You have completed block ',num2str(BLOCK),' of ',num2str(EXC_ST.NUMPBLOCK))};
        TEXT(2)={ConText{1}};
        TEXT(3)={ConText{2}};
        TEXT(4)={''};
        TEXT(5)={'Press any key to proceed'};
        
    case {'EXEND'},% Ending Experiment screen
        TEXT(1)={'Thank you'};
        TEXT(2)={'You have completed the practice phase'};
        TEXT(3)={''};
        TEXT(4)={''};
        TEXT(5)={'Press any key to proceed'};
end

%Drawing
for i=1:size(TEXT,2)
    if strcmp(msgtime,'BRULE'),adj=50;else adj=100;end;
    text=TEXT{i};xp=STIM_ST.XC-350;yp=STIM_ST.YC-200+(adj*(i-1));
    if i==size(TEXT,2),yp=STIM_ST.YC-200+400;end
    Screen('DrawText',W,text,xp,yp,COLOR_ST.BLACK);
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
global PARAM_ST BLOCK TRIAL PRACEXC
%**************************************************************************

switch msgtime
    case {'BINTRO'},% Begining Block BREAK screen
        if ~PRACEXC
            ConText(1)={'ATTENTION!! This is a PRACTICE!!'};
            ConText(2)={'Please remember a sequence of orientations'};
        else
            ConText(1)={'ATTENTION!! This is an EXPERIMENT!!'};
            ConText(2)={'Please remember a sequence of chunks'};
        end
        
    case {'BRESULT'},% Begining Block BREAK screen
        %Change Setting for Practice and Experiment
        trial=EXC_ST.NUMPTRIAL;
        accScore=length(find(RECORD_ST.ACC(:)==0));
        ConText(1)={strcat('You missed ',num2str(accScore),' of ',num2str(trial),'trials')};ConText(2)={''};
        if ~PRACEXC && (accScore/trial*100)<TASK_ST.ERRORCUT,ConText(2)={'Please go on to next practice'};end
        if ~PRACEXC && (accScore/trial*100)>=TASK_ST.ERRORCUT,ConText(2)={'Too many mistakes!! Please repeat practice!!'};end
end


end

%##########################################################################

function locKeepRecordM(dataArray,workT)
%==========================================================================
% locKeepRecordM records the title of text files for the final output.
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
if ~isstruct(dataArray),error('locKeepRecordM assums struct array!');end;
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
    if size(var,1)==EXC_ST.NUMPBLOCK,varType='BCOND';ref=[1,1];end
    if size(var,2)==EXC_ST.NUMPTRIAL,varType='TCOND';ref=[1,TRIAL];end;
    
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
global EXC_ST TASK_ST PARAM_ST FPNAME_ST STIM_ST TIME_ST COLOR_ST KEY_ST
%**************************************************************************

% Setting Rule of TaskElenent!!++++++++++++++++++++++++++++++++++++++++++++
% Creating template of hieararchy!!++++++++++++++++++++++++++++++++++++++++

for c=1:size(TASK_ST.CHUNKPAIRS,1)
    scode=strcat('S',num2str(c));
    c1=TASK_ST.BASECHUNK(TASK_ST.CHUNKPAIRS(c,1),:);
    c2=TASK_ST.BASECHUNK(TASK_ST.CHUNKPAIRS(c,2),:);
    TASK_ST.RCODET.(strcat('LV2_',scode))=horzcat(c1,c2);
end

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
IMAGE_ST.IMARRAY(1,1)={'Bar_A.jpg'};
IMAGE_ST.IMARRAY(2,1)={'Bar_B.jpg'};
IMAGE_ST.IMARRAY(3,1)={'Bar_C.jpg'};

% Visual stimuli preparation
for i=1:size(IMAGE_ST.IMARRAY,1)
    %reading jpg image data
    imaCode = char([IMAGE_ST.IMFPATH filesep IMAGE_ST.IMARRAY{i,1}]);
    IMAGE_ST.IMARRAY{i,2}=imread(strtrim(imaCode));
    IMAGE_ST.IMARRAY{i,3}=Screen('MakeTexture',W,IMAGE_ST.IMARRAY{i,2});
end

end

%##########################################################################

function  locRadialCSN
%==========================================================================
%Radial Constallation
%Note: specific sets of constalleations dependin on the number of tasksets
%VERY TEDIOUS FOR LOOP STRUCTURE! THIS SHOULD BE VECTORIZED COMPUTATIONALLY!
%==========================================================================
global FPNAME_ST EXC_ST TASK_ST IMAGE_ST COLOR_ST RECORD_ST STIM_ST W
%**************************************************************************

% Creating position XY Codes
numelm=length(TASK_ST.BASECHUNK(1,:));
for i=1:numelm
    angle=(i-1)*(360/numelm);
    STIM_ST.RCXYREF(1,i)=ceil(STIM_ST.XC+sin(((-180+angle)*pi/180))*STIM_ST.RADSIZE);
    STIM_ST.RCXYREF(2,i)=ceil(STIM_ST.YC-STIM_ST.YCOFF+cos(((-180+angle)*pi/180))*STIM_ST.RADSIZE);
end

end


%##########################################################################

function [path]=fileroot(fileName)
path=which(fileName); %path is to file
i=find(filesep==path); %file seperator
path=path(1:i(end));
end
