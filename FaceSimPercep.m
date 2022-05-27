function FaceSimPercep
% Perception-based similarity rating task

% Participants are presented with pairs of faces and are asked to rate the
% visual similarity between them using a 7-point scale (1 means very
% different and 7 means very similar). The pairwise similarity ratings will
% be used as input for the perception-based image reconstruction. 

% Note: the trial order file for a subject must be generated prior to
% running this script (see tpOrder.m)

%% 
close all
clc
sca

try
%% Open debug mode or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debugMode = 0; % change it to 1 for debugging mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize the random number generator 
% rand('state', sum(100*clock)); 
rng('shuffle')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect experiment information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~debugMode 
%------------------------ If Debug Mode is OFF ---------------------------%
    % Input box
    prompt = {'Group:','Subject ID:','Session','Run ID:','Practice (Y/N)?','Study'};
    boxTitle = 'ExpInfo - Similarity Rating';
    boxSize = [1 100];
    defaultAnswer = {'ips','99','1','1','Y','ips'};
    opts.Resize = 'on';
    answ = inputdlg(prompt,boxTitle,boxSize,defaultAnswer,opts);
    
    % Group label: different integers can be assigned for different groups
    % of participants
    group = answ{1,1};
    checkError = 1;
    while checkError
        if isempty(group)
            disp('Not input a group label. Please try again!')
            answ_group = inputdlg('Input a group label:','Group',boxSize,defaultAnswer(1),opts);
            group = answ_group{1,1};           
        else
            checkError = 0;
        end
    end
%     if strcmp(group,'c') % converting group label to a number for coding covenience
%         groupNum = 1;
%         groupID = 'c';
%     elseif strcmp(group,'p')
%         groupNum = 2;
%         groupID = 'p';        
%     end
    groupNum = 0;
    groupID = 'ips'; % in-person subject
        
    % Subject ID 
    subID = answ{2,1};
    checkError = 1;
    while checkError
        if isempty(subID)
            disp('Not input a subID. Please try again!')
            answ_subID = inputdlg('Input subject ID again:','Subject ID',boxSize,defaultAnswer(2),opts);
            subID = answ_subID{1,1};
        elseif str2num(subID)<10 && ~strcmp(subID(1),'0')
            disp('You input a single digit number. Please add zero in the tens'' place')
            answ_subID = inputdlg('Please add zero in the tens'' place:','Subject ID',boxSize,defaultAnswer(2),opts);
            subID = answ_subID{1,1};
        else
            checkError = 0;
        end
    end
    subNum = str2num(subID); % converting subID to a number for coding covenience
    
    % Which session this task is. This and the run ID will determine which
    % blocks will be run 
    sessionID = answ{3,1};
    checkError = 1;
    while checkError
        if isempty(sessionID)
            disp('Not input session number. Please try again!')
            answ_sessionID = inputdlg('Input a session number again:','Session',boxSize,defaultAnswer(3),opts);
            sessionID = answ_sessionID{1,1};
        else
            checkError = 0;
        end
    end
    sessionNum = str2num(sessionID); % converting sessionID to a number for coding covenience
    
    % Which run this task is. This and the session number will determine
    % which blocks will be run
    runID = answ{4,1};
    checkError = 1;
    while checkError
        if isempty(runID)
            disp('Not input runID. Please try again!')
            answ_runID = inputdlg('Input a run ID again:','Run ID',boxSize,defaultAnswer(4),opts);
            runID = answ_runID{1,1};
        else
            checkError = 0;
        end
    end
    runNum = str2num(runID); % converting runID to a number for coding covenience
    
    % Run practice block or not
    practice = answ{5,1};
    checkError = 1;
    while checkError
        if isempty(practice)
            disp('Please input Y/N to decide if running practice trials')
            answ_practice = inputdlg('Input Y/N to decide practice or not:','Practice',boxSize,defaultAnswer(5),opts);
            practice = answ_practice{1,1};
        elseif ~strcmp(practice,'Y') && ~strcmp(practice,'N')
            disp('Please input Y or N only')
            answ_practice = inputdlg('Input Y/N to decide practice or not:','Practice',boxSize,defaultAnswer(5),opts);
            practice = answ_practice{1,1};
        elseif runNum==1 && strcmp(practice,'N') % default is that only 1st run in each session requires practice
            disp('Note: Default is that the first run in each session requires practice.')
            strResponse_prac = input('Proceed without practice trials for the first run? Y/N: ','s');
            if strcmp(strResponse_prac,'Y')
                disp('Not running practice trials for the first run in this session.')
                checkError = 0;
            elseif strcmp(strResponse_prac,'N')
                disp('Will run practice trials for this run.')
                practice = 'Y';
                fprintf('\n')
            end
        elseif runNum~=1 && strcmp(practice,'Y')
            disp('Note: Default is that ONLY the first run in each session requires practice.')
            strResponse_prac = input('Do practice trials in other runs? Y/N: ','s');
            if strcmp(strResponse_prac,'Y')
                disp('Running practice trials for this run in this session.')
                checkError = 0;
            elseif strcmp(strResponse_prac,'N')
                disp('Will NOT run practice trials for this run.')
                practice = 'N';
                fprintf('\n')
            end
        else
            checkError = 0;
        end
    end
    
    % Study name
    study = answ{6,1};
    checkError = 1;
    while checkError
        if isempty(study)
            disp('Not input the study''s name. Please try again!')
            answ_study = inputdlg('Input the study''s name:','Study',boxSize,defaultAnswer(6),opts);
            study = answ_study{1,1};           
        else
            checkError = 0;
        end
    end
     
%-------------------------------------------------------------------------%

else
%------------------------ If Debug Mode is ON ----------------------------%
    strResp_debug = input('Execute Debug Mode? Y/N: ','s'); % confirm if running debug mode
    while ~strcmp(strResp_debug,'Y') && ~strcmp(strResp_debug,'N')
        disp('Please input Y/N: ')
        strResp_debug = input('Execute Debug Mode? Y/N: ','s');
    end
    
    if strcmp(strResp_debug,'Y') % if YES, load exp info for debug mode
        group = 'debug';
        groupNum = 0;
        groupID = 'debug';
        subID = '00';
        subNum = 0;
        sessionNum = 1;
        runID = '1';
        runNum = 1;
        practice = 'Y';
        study = 'ips';
        
        disp(['Group ID: ' groupID])
        disp(['Subject ID: ' subID])
        disp(['sessionNum: ' num2str(sessionNum)])
        disp(['runNum: ' num2str(runNum)])
        disp(['Practice or not: ' practice])
        
        strResp_debugInfo = input('Are all experiment information for the debug mode okay? Y/N: ','s');
        while ~strcmp(strResp_debugInfo,'Y') && ~strcmp(strResp_debugInfo,'N')
            disp('Please input Y/N: ')
            strResp_debugInfo = input('Are all experiment information for the debug mode okay? Y/N: ','s');
        end
        if strcmp(strResp_debugInfo,'N')
            disp('Change the experiment information for the debug mode in the script.')
            return
        end
           
    else % otherwise, quit the experiment
        disp('Please change the debug variable to 0 in the script.')
        return
    end

end

% Summarize subject and exp info
expInfo.groupNum = groupNum;
expInfo.groupID = groupID; 
expInfo.subID = subID;
expInfo.subNum = subNum;
expInfo.sessionNum = sessionNum;
expInfo.runNum = runNum;
expInfo.practice = practice;
expInfo.study = study;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up stimuli lists 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set directories --------------------------------------------------------%
scriptName = mfilename;
p = mfilename('fullpath');
exp_folder = p(1:end-length(scriptName)); 
cd(exp_folder);

% Get stimuli set for the experiment -------------------------------------%
% A cell array storing face images (nImages x 1)
% Each cell is a face image (in my study, it's a 172 x 115 x 3 matrix
% (uint8)) (i.e., a color face image with height and width of 172 and 115 
% pixels, respectively)
stim_file= fullfile(exp_folder, 'ims_new.mat'); 
load(stim_file); % variable name: ims_new

% A cell array storing the names of face images (nImages x 1)
% Each cell is the name of a face image (in my study, it's a number
% representing the ID of the face image)
name_file= fullfile(exp_folder, 'name_new.mat'); 
load(name_file) % variable name: name_new

nTrials = length(ims_new); % number of faces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Order of Trials during the actual experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The order file must be generated in advance using the script tpOrder.m
% Load the order file
order_file = fullfile(exp_folder,['order_tp_' study],['tp_' groupID subID '.mat']); % mat name: order_blkord
if ~exist(order_file,'file')
    disp('order file does not exist.')
    % [~]=tpOrder(subID);
    return
else
    load(order_file)
end

% Set the start and the end block, get the pairs that will be presented
if sessionNum==1 && runNum==1 %part==1
    startBlock = 1;
    endBlock = 3;
elseif sessionNum==1 && runNum==2
    startBlock = 4;
    endBlock = 7;
elseif sessionNum==2 && runNum==1 
    startBlock = 8;
    endBlock = 10;
elseif sessionNum==2 && runNum==2 % repeated block (trials are the same as 4th block)
    startBlock = 4;
    endBlock = 4;
elseif sessionNum==2 && runNum==3
    startBlock = 11;
    endBlock = 14;
end

% Trial order for the blocks being conducted in the current session and
% run
nPairs = size(pair,1); % number of all pairs of faces
breakAfterTrials = nPairs/14; % number of trials in one block; the total number of pairs is divided by the total number of blocks
randomizedTrials = 1+(startBlock-1)*breakAfterTrials : endBlock*breakAfterTrials;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultsDir = fullfile(['results_' study],[groupID subID]);
% if the subject result folder doesn't exist, then create one
if ~exist(resultsDir,'dir') 
    mkdir(resultsDir);
end
% if the subject archive folder doesn't exit, then create one
archiveDir = fullfile(['results_' study],'archive',[groupID subID]);
if ~exist(archiveDir,'dir') 
    mkdir(archiveDir);
end

% Define information for the output files
fileHeader = ['Group\t' 'subID\t' 'session\t' 'runID\t' 'trial\t' 'block\t' ...
    'imageFileName1\t' 'imageFileName2\t' 'response\t' 'RT\n']; 
% note: the image file name corresponds to the stimulus file name (i.e.,
% from the name_new matrix; there won't be 23, 31 and 44 b/c they are used
% as learned faces for memory reconstruction in my MA project)
fileFormatSpec = ['%d\t' '%d\t' '%d\t' '%d\t' '%d\t' '%d\t' ...
    '%d\t' '%d\t' '%d\t' '%f\n'];
matColNames.resVarName = {'Group', 'subID', 'sessionNum', 'runID', 'trial', 'block', ...
    'imageFileName1', 'imageFileName2', 'response', 'RT'};
matColNames.durVarName = {'Fix', 'Images', 'TrialLength'};

%---------------------------- For Practice -------------------------------%
if strcmp(practice,'Y') % if running the practice 
    disp('Setting up output files for the practice block...')
    filename_prac = ['resultfile_tp_' groupID subID '_sess' num2str(sessionNum) '_run' runID '_Prac'];
    [outputfile_prac,output_prac] = setOutFile(resultsDir,archiveDir,filename_prac,...
        fileHeader,expInfo,matColNames);
    filenameTxt_prac = fullfile(resultsDir,[filename_prac '.txt']); % full path of the output text file
    filenameMat_prac = fullfile(resultsDir,[filename_prac '.mat']); % full path of the output mat file 
    fprintf('\n')
end % end of the if statement whether running practice or not

%------------------------- For Actual Blocks -----------------------------%
disp('Setting up output files for the actual blocks...')
filename = ['resultfile_tp_' groupID subID '_sess' num2str(sessionNum) '_run' runID];
[outputfile,output] = setOutFile(resultsDir,archiveDir,filename,...
    fileHeader,expInfo,matColNames);
filenameTxt = fullfile(resultsDir,[filename '.txt']); % full path of the output text file
filenameMat = fullfile(resultsDir,[filename '.mat']); % full path of the output mat file
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input devices setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Keyboard setup %%%%
KbName('UnifyKeyNames');
% Response keys (optional; for no subject response use empty list {})
% Type KbName in command window to find out key's name 
responseKeys = {'1!','2@','3#','4$','5%','6^','7&'}; 

% Keys used for non-task-related purpose (e.g., proceed, quit, redo
% buttons)
KbCheckList = [KbName('space'),KbName('ESCAPE'),KbName('/?')];
for i = 1:length(responseKeys)
    KbCheckList = [KbName(responseKeys{i}),KbCheckList];
end
RestrictKeysForKbCheck(KbCheckList);
[keyboardIndices, productNames, allInfos] = GetKeyboardIndices;

% %%%% Mouse setup %%%%
% [mouseIndices, productNames, allInfo] = GetMouseIndices;
% 
% %%%% Touchscreen setup %%%%
% [touchIndices, productNames, allInfo] = GetTouchDeviceIndices;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Task parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Duration
duration_fix = 0.5; % 500 ms
duration_img = Inf;

% Fixation cross color
fixColor = [255, 255, 255]; % white

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up Psychtoolbox (don't modify this section)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1); % skip Sync Test (for OSX)
AssertOpenGL;

% Screen setup
nScreen = length(Screen('Screens')); % number of monitors
if nScreen==1
    % When there is only 1 monitor
    whichScreen = Screen('Screens');
else
    % When there is more than 1 monitor
    
    % ScreenNumber==0 is the screen with the menu bar (you can change it in
    % System Preferences/Display). Screen with other numbers are considered
    % as external screens. Usually experiments are done on an external
    % screen, so the maximum of screen numbers is usually a good guess to
    % select the external screen you want to run your experiment
    % If the program is run on the computer's screen rather than the
    % external screen, make sure the computer's screen is the one with menu
    % bar
    whichScreen = max(Screen('Screens')); % Define which screen to display
    
%     % Alternatively, you can call the ScreenTest function to test the
%     % performance of each monitor and determine which monitor you want to
%     % use
%     ScreenTest
%     disp('Please select which screen you want to use to run the experiment.')
%     screenNum = input('0/1/2/...?: ');
%     whichScreen = screenNum;
end

% Define background colour
blackBkg = BlackIndex(whichScreen);
whiteBkg = WhiteIndex(whichScreen);
grayBkg = GrayIndex(whichScreen);
backgroundColor = blackBkg; % define background colour

if ~debugMode % No debug mode   
    [window1, rect] = PsychImaging('OpenWindow', whichScreen, backgroundColor);% Open the screen
    % [window1, rect] = Screen('Openwindow',whichScreen,backgroundColor,[],[],2);
    W = rect(RectRight); % screen width in pixels (=screenXpixels)
    H = rect(RectBottom); % screen height in pixels (=screenYpixels)
    % [screenXpixels, screenYpixels] = Screen('WindowSize', window1);
    [xCenter, yCenter] = RectCenter(rect); % Get the centre coordinate of the window (xCenter=W/2; yCenter=H/2)
    
    % Set priority level
    Priority(MaxPriority(window1));
    % Priority(2);
    
    % Fill the screen with background colour
    Screen(window1,'FillRect',backgroundColor);
    Screen('Flip', window1);
    
    % Enable listening, additionally any output of keypresses to Matlabs or 
    % Octaves windows is suppressed
    ListenChar(2);
    
    % Hide the cursor
    HideCursor; 

%------------------------ If Debug Mode is ON ----------------------------%
else
    PsychDebugWindowConfiguration
    [window1, rect] = PsychImaging('OpenWindow', whichScreen, backgroundColor);% Open the screen
    W = rect(RectRight); % screen width in pixels (=screenXpixels)
    H = rect(RectBottom); % screen height in pixels (=screenYpixels)
    % [screenXpixels, screenYpixels] = Screen('WindowSize', window1);
    [xCenter, yCenter] = RectCenter(rect); % Get the centre coordinate of the window (xCenter=W/2; yCenter=H/2)
    
    % Set priority level
    Priority(MaxPriority(window1));
    
%     sca
%     sreenSize_debugMode = [0 0 W/2 H/2];
%     [window1, rect] = PsychImaging('OpenWindow', whichScreen, backgroundColor,sreenSize_debugMode);% Open smaller screen
%     W = rect(RectRight); % screen width in pixels (=screenXpixels)
%     H = rect(RectBottom); % screen height in pixels (=screenYpixels)
%     [xCenter, yCenter] = RectCenter(rect); % Get the centre coordinate of the window
    
    % Fill the screen with background colour
    Screen(window1,'FillRect',backgroundColor);
    Screen('Flip', window1);
end
%-------------------------------------------------------------------------%

% Interframe interval (1/refreshing rate)
ifi = Screen('GetFlipInterval', window1);
% slack = Screen('GetFlipInterval', window1)/2; % getflipinterval: know the interval btw. flip

% Text font and size for instruction page
Screen('Preference','TextRenderer', 0); % 0: use system legacy text renderer
Screen('TextFont',window1, 'Arial');
Screen('TextSize',window1, 18);

textColor = [1 1 1]; %[0 0 0]; %255; % Text color: choose a number from 0 (black) to 255 (white)
textColor_important = [1 0 0];%[1 1 0]; %[255 255 85]; % used when displaying important message (yellow)

inst_xpos = 50;
top_ypos = 50;
inst_ypos = 30;

% Save the screen information
screenInfo.W = W;
screenInfo.H = H;
screenInfo.xCenter = xCenter;
screenInfo.yCenter = yCenter;
screenInfo.ifi = ifi;
output.screenInfo = screenInfo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run experiment - Instruction page
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DrawFormattedText(window1,'In this task you will be presented with pairs of faces.', inst_xpos, top_ypos+inst_ypos, textColor);
DrawFormattedText(window1,'You will be asked to rate how similar the two faces are on a 7-point scale,', inst_xpos, top_ypos+inst_ypos*2, textColor);
DrawFormattedText(window1,'for which 1 means that the two faces are very dissimilar in appearance', inst_xpos, top_ypos+inst_ypos*3, textColor);
DrawFormattedText(window1,'and 7 means that the two faces are very similar.', inst_xpos, top_ypos+inst_ypos*4, textColor);

drawScale(window1,W/2,H/2-25); %drawScale(window1,W/4,H/2-25);

DrawFormattedText(window1,'Please respond by speaking out your rating,', inst_xpos, top_ypos+inst_ypos*10, textColor);
DrawFormattedText(window1,'and try to use the whole range of the scale.', inst_xpos, top_ypos+inst_ypos*11, textColor);

if strcmp(practice,'Y')
    DrawFormattedText(window1,'Press the space bar when you are ready to begin the practice!', inst_xpos, top_ypos+inst_ypos*13, textColor);
else
    DrawFormattedText(window1,'Press the space bar when you are ready to begin!', inst_xpos, top_ypos+inst_ypos*13, textColor);
end

Screen('Flip',window1);

% Wait for subject to press spacebar
while 1
    [keyIsDown,secs,keyCode] = KbCheck(keyboardIndices);
%     [keyIsDown,secs,keyCode]=PsychHID('KbCheck',6);
    if keyCode(KbName('space'))==1
        break
    end
    
    % ESC key quits the experiment
    if keyCode(KbName('ESCAPE')) == 1
        %             clear all
        %             close all
        sca
        disp('Experiment Terminated')
        ListenChar; % enable listening again
        ShowCursor;
        return;
    end
end
WaitSecs(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run experiment - Practice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only execute the following codes when practice trials are required
if strcmp(practice,'Y') 
    % pick up 12 trials from the last block for practice
    if sessionNum==1
        pracStart = 14;
    elseif sessionNum==2
        pracStart = 7; 
    end
    pracord = 1+(pracStart-1)*breakAfterTrials : 1+(pracStart-1)*breakAfterTrials+11; 
    doPrac = 1; % indicator to start practice while loop
    nPrac = 1; % number of practice block
    cTrial = 1; % count trials
    blockNum = (nPrac-1)+pracStart; % block ID among all blocks (14 blocks in total)
else
    doPrac = 0;
end

%%%% Practice starts %%%%
while doPrac
    if nPrac==1
        Screen('BlendFunction', window1, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    end
    
    %---------------------------------------------------------------------%
    %%%% Trail start %%%%
    %---------------------------------------------------------------------%
    for t = pracord % go through every practice trials 
        % Get face images to be presented for the current trial ----------------%
        ind = pair(t,:); %pair(randomizedTrials(t),:);
        file1 = name_new{ind(1),1}; %num2str(name_new{ind(1),1}); %['stim' num2str(name{t,1})];
        file2 = name_new{ind(2),1}; %num2str(name_new{ind(2),1}); %['stim' num2str(name{t,2})];
        img1 = ims_new{ind(1),1};
        img2 = ims_new{ind(2),1};
        imageSize1 = size(img1);
        imageSize2 = size(img2);
        imgPos1 = [((W-imageSize1(2))/2)-100 (H-imageSize1(1))/2 ((W+imageSize1(2))/2)-100 (H+imageSize1(1))/2];
        imgPos2 = [((W-imageSize2(2))/2)+100 (H-imageSize2(1))/2 ((W+imageSize2(2))/2)+100 (H+imageSize2(1))/2];
        %-----------------------------------------------------------------------%
        
        % Show fixation cross --------------------------------------------------%
        acceptResp_fix = 0;
        [fix_start, fix_end, ~, ~] = showCross(window1, W, H,...
            fixColor, ifi, duration_fix, acceptResp_fix);
        % Blank screen
        Screen(window1, 'FillRect', backgroundColor);
        Screen('Flip', window1, fix_end+(1-0.5)*ifi);
        %-----------------------------------------------------------------------%
        
        % Show two face images and a 7-point scale -----------------------------%
        acceptResp_img = 1;
        goAfterResp_img = 1;
        [imgTwoStart, imgTwoEnd, rt_img, resp_img] = showStimAndScale(window1, img1, img2,...
            imgPos1, imgPos2, ifi, duration_img, acceptResp_img, ...
            responseKeys, keyboardIndices, goAfterResp_img, screenInfo);
        % recode the response output
        if isempty(resp_img)
            respCode_img = 0;
        else
            keyPress = nan(1,length(responseKeys));
            for i = 1:length(responseKeys)
                keyPress(i) = strcmp(resp_img,responseKeys{i});
            end
            if sum(keyPress)>1 % press more than 1 response key
                respCode_img = 0;
            else
                if find(keyPress)==1
                    respCode_img = 1;
                elseif find(keyPress)==2
                    respCode_img = 2;
                elseif find(keyPress)==3
                    respCode_img = 3;
                elseif find(keyPress)==4
                    respCode_img = 4;
                elseif find(keyPress)==5
                    respCode_img = 5;
                elseif find(keyPress)==6
                    respCode_img = 6;
                elseif find(keyPress)==7
                    respCode_img = 7;
                end
            end
        end
        % Blank screen
        Screen(window1, 'FillRect', backgroundColor);
        Screen('Flip', window1, imgTwoEnd+(1-0.5)*ifi);
        %-----------------------------------------------------------------------%
        
        % Save results to file--------------------------------------------------%
        % save results to text file
        fprintf(outputfile_prac, fileFormatSpec,...
            groupNum, subNum, sessionNum, runNum, cTrial, blockNum, ...
            file1, file2, respCode_img, rt_img);
        % save results to a mat file
        output_prac.results(cTrial,:) = [groupNum, subNum, sessionNum, runNum, cTrial, blockNum, ...
            file1, file2, respCode_img, rt_img];
        % save duration of fixation and images to a mat file
        trlLngth = imgTwoEnd - fix_start; % trial length
        output_prac.Duration(cTrial,:) = [fix_end-fix_start, ...
            imgTwoEnd-imgTwoStart, ...
            trlLngth];
        % save into a mat file
        save(filenameMat_prac,'output_prac');
        %-----------------------------------------------------------------------%
        
        % update trial counts
        cTrial = cTrial+1;
        
    end % end of all practice trials
    
    %---------------------------------------------------------------------%
    %%%% Practice End Page %%%%
    %---------------------------------------------------------------------%
    DrawFormattedText(window1,'The end of practice trials. Any question?', inst_xpos, top_ypos+inst_ypos, textColor);
    Screen('Flip',window1);
    
    % Wait for subject to press spacebar
    while 1
        [keyIsDown,secs,keyCode] = KbCheck(keyboardIndices);
        if keyCode(KbName('space'))==1
            doPrac = 0;
            break
        end
        
        if keyCode(KbName('/?'))==1
%             outputfile_prac = fopen(filename_prac,'a'); %append output to the end of previous practice data
            nPrac = nPrac+1; % number of practice block
            break
        end
        
        % ESC key quits the experiment
        if keyCode(KbName('ESCAPE')) == 1
%             clear all
%             close all
            ListenChar(1);
            sca
            return;
        end
    end
    WaitSecs(1);
    
end % end of while statement of running practice

if strcmp(practice, 'Y')
    fclose(outputfile_prac); % close the output file for actual trials
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run experiment - Actual task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Run experimental trials %%%%
% if practice trials are not required, skip the above practice block and
% execute the following codes after the instruction page

%%%% Actual task starts %%%%
if strcmp(practice,'N')
    Screen('BlendFunction', window1, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
end

cTrial = 1; % count trials
cBlk = 1; % count blocks
blockNum = (cBlk-1)+startBlock; % block ID among all blocks (14 blocks in total)
for t = randomizedTrials % go through every trials in the selected blocks
  % Get face images to be presented for the current trial ----------------%
  ind = pair(t,:); %pair(randomizedTrials(t),:);
  file1 = name_new{ind(1),1}; %num2str(name_new{ind(1),1}); %['stim' num2str(name{t,1})];
  file2 = name_new{ind(2),1}; %num2str(name_new{ind(2),1}); %['stim' num2str(name{t,2})];
  img1 = ims_new{ind(1),1};
  img2 = ims_new{ind(2),1};
  imageSize1 = size(img1);
  imageSize2 = size(img2);
  imgPos1 = [((W-imageSize1(2))/2)-100 (H-imageSize1(1))/2 ((W+imageSize1(2))/2)-100 (H+imageSize1(1))/2];
  imgPos2 = [((W-imageSize2(2))/2)+100 (H-imageSize2(1))/2 ((W+imageSize2(2))/2)+100 (H+imageSize2(1))/2];
  %-----------------------------------------------------------------------%
  
  % Show fixation cross --------------------------------------------------%
  acceptResp_fix = 0;
  [fix_start, fix_end, ~, ~] = showCross(window1, W, H,...
      fixColor, ifi, duration_fix, acceptResp_fix);
  % Blank screen
  Screen(window1, 'FillRect', backgroundColor);
  Screen('Flip', window1, fix_end+(1-0.5)*ifi);
  %-----------------------------------------------------------------------%
  
  % Show two face images and a 7-point scale -----------------------------%
  acceptResp_img = 1;
  goAfterResp_img = 1;
  [imgTwoStart, imgTwoEnd, rt_img, resp_img] = showStimAndScale(window1, img1, img2,...
      imgPos1, imgPos2, ifi, duration_img, acceptResp_img, ...
      responseKeys, keyboardIndices, goAfterResp_img, screenInfo);
  % recode the response output
  if isempty(resp_img)
      respCode_img = 0;
  else
      keyPress = nan(1,length(responseKeys));
      for i = 1:length(responseKeys)
          keyPress(i) = strcmp(resp_img,responseKeys{i});
      end
      if sum(keyPress)>1 % press more than 1 response key
          respCode_img = 0;
      else
          if find(keyPress)==1
              respCode_img = 1; 
          elseif find(keyPress)==2
              respCode_img = 2; 
          elseif find(keyPress)==3
              respCode_img = 3; 
          elseif find(keyPress)==4
              respCode_img = 4; 
          elseif find(keyPress)==5
              respCode_img = 5; 
          elseif find(keyPress)==6
              respCode_img = 6; 
          elseif find(keyPress)==7
              respCode_img = 7; 
          end
      end
  end
  % Blank screen
  Screen(window1, 'FillRect', backgroundColor);
  Screen('Flip', window1, imgTwoEnd+(1-0.5)*ifi);
  %-----------------------------------------------------------------------%
  
  % Save results to file--------------------------------------------------% 
  % save results to text file  
  fprintf(outputfile, fileFormatSpec,...
      groupNum, subNum, sessionNum, runNum, cTrial, blockNum, ...
      file1, file2, respCode_img, rt_img);
  % save results to a mat file
  output.results(cTrial,:) = [groupNum, subNum, sessionNum, runNum, cTrial, blockNum, ...
      file1, file2, respCode_img, rt_img];
  % save duration of fixation and images to a mat file
  trlLngth = imgTwoEnd - fix_start; % trial length
  output.Duration(cTrial,:) = [fix_end-fix_start, ...
      imgTwoEnd-imgTwoStart, ...
      trlLngth];
  % save into a mat file
  save(filenameMat,'output');
  %-----------------------------------------------------------------------%
  
  % Provide a short bread if the current block is not the end of the task
  if mod(cTrial,breakAfterTrials)==0 && cTrial~=length(randomizedTrials) % if not the end of the task
      % Update block counts
      cBlk = cBlk+1;
      blockNum = (cBlk-1)+startBlock;
      
      % Provide a short break after a certain number of trials
      DrawFormattedText(window1,'Break time. Please let me know when you are ready to continue', inst_xpos, top_ypos+inst_ypos, textColor);
      Screen('Flip',window1);
      
      % Wait for subject to press spacebar
      while 1
          % ESC key quits the experiment
          if keyCode(KbName('ESCAPE')) == 1
              %     clear all
              %     close all
              ListenChar(1);
              sca
              disp('Experiment Terminated')
              return;
          end
          
          [keyIsDown,secs,keyCode] = KbCheck(keyboardIndices);
          if keyCode(KbName('space')) == 1
              break
          end
      end
      WaitSecs(1);     
  end
  
  % update trial counts
  cTrial = cTrial+1; 
       
end % end of all trials in the to-be presented blocks

fclose(outputfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Thank you page
DrawFormattedText(window1,'The task is over. Thank you.', inst_xpos, top_ypos+inst_ypos, textColor);
Screen('Flip',window1);
WaitSecs(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% End the experiment (don't change anything in this section)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Priority(0);
RestrictKeysForKbCheck([]);
ListenChar; % enable listening again
ShowCursor;
Screen(window1,'Close');
close all
sca;
return

catch ME
    disp('Something wrong happened')
    Priority(0);
    ShowCursor;
    RestrictKeysForKbCheck([]);
    ListenChar(1);
    sca;
    %pfp_ptb_cleanup
    rethrow(ME)
    return
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%-------------------------------------------------------------------------%
% Set up output files
%-------------------------------------------------------------------------%
function [outputfile,output] = setOutFile(resultsDir,archiveDir,filename,...
    fileHeader,expInfo,matColNames)
% Input:
%   resultDir - the folder where the output files will be saved
%   archiveDir - the foler where the archived files will be saved
%   filename - the name of the output text and mat files (without filename
%              extension, e.g., .txt)
%   fileHeader - the header in the output text file
%   expInfo - a struct containing the current subject and experiment 
%             information
%   matColNames - column names of matrices in the output struct
% Outputs:
%   outputfile - fileID of the output text file
%   output - a struct containing all results and task information

filenameTxt = fullfile(resultsDir,[filename '.txt']); % path of the output text file
% filenameMat = fullfile(resultsDir,[filename '.mat']); % path of the output mat file
groupID = expInfo.groupID; 
subID = expInfo.subID;

% Check if the subject output text file already exists
if exist(filenameTxt,'file')~=0 % if the subject output already exists
    disp(['Datafile already exists for ' groupID num2str(subID)]);
    disp('Press 1 to overwrite, 2 to move to archive, or 3 to exit (followed by enter)!')
    strResp = input('1/2/3: ');
    while strResp~=1 && strResp~=2 && strResp~=3
        strResp = input('Please input 1/2/3:');
    end
    if strResp==3
        return        
    elseif strResp==2
        try
            movefile(filenameTxt,archiveDir);
        catch
            disp(['Failed to move the file ' filename]);
        end
    end    
end

% Open the output text file
outputfile = fopen(filenameTxt,'wt');
fprintf(outputfile, '%s %s \n','Test Date:', datestr(now));
fprintf(outputfile,'%s %s \n','Matlab version ', version);
PTversion = PsychtoolboxVersion;
fprintf(outputfile,'%s %s \n','Psychtoolbox version ', PTversion);
fprintf(outputfile,'%s \n','<><><><><><><><>');
fprintf(outputfile,fileHeader);

% Prepare the output mat file
output.expInfo = expInfo;
output.testDate = clock;
output.matColNames = matColNames;
% output.resVarName = matColNames.resVarName;
% output.durVarName = matColNames.durVarName;
% output.eventOnsetName = matColNames.eventOnsetName;
% save(filenameMat,'output')

end

%%
%-------------------------------------------------------------------------%
% Record responses
%-------------------------------------------------------------------------%
function [rt,resp] = recordResp(startTime, responseKeys, keyboardIndices)
% input
%   startTime - time point of the stimulus onset
%   stiDuration - length of stimulus presentation
%   respDuration - length of interval participants can respond after
%                  stimulus disappears
% output
%   rt - reaction time (default = 0)
%   resp - response key (default = [])

rt = 0;
resp = [];

[keyIsDown,secs,keyCode] = KbCheck(keyboardIndices);
respTime = secs; %GetSecs;
pressedKeys = find(keyCode);

% ESC key quits the experiment
if keyCode(KbName('ESCAPE')) == 1
%     clear all
%     close all
    sca
    disp('Experiment Terminated')
    return;
end

% Check for response keys
if ~isempty(pressedKeys)
    for i = 1:length(responseKeys)
        if KbName(responseKeys{i}) == pressedKeys(1)
            resp = responseKeys{i};
            rt = respTime - startTime; %respTime - scaleStart;
        end
    end
end

end

%%
%-------------------------------------------------------------------------%
% Draw a fixation cross (overlapping horizontal and vertical bar)
% (presented before target presentation)
% And record response if need be
%-------------------------------------------------------------------------%
function [crossStart, crossEnd, rt, resp] = showCross(window, W, H, barColor,...
    ifi, duration, acceptResp, responseKeys, keyboardIndices, goAfterResp)
% input
%   window - window of the experiment
%   W - width of the screen
%   H - height of the screen
%   barColor - the color of fixation cross (e.g., black:[0 0 0]; white:[255
%              255 255]; red:[255 0 0]; yellow:[255 255 0]
%   ifi - interframe interval (1/refreshing rate)
%   duration - length of the cross presentation
%   acceptResp - whether a response is required to be made (=1) or not (=0)
%   responseKeys - designated response keys
%   keyboardIndices - index of the input keyboard
%   goAfterResp - whether move on to the next trial once a response is made
%                 (1-yes; 0-no)
% output
%   crossStart - starting time point of the fixation cross presentation
%   crossEnd -  ending time point of the fixation cross presentation
%   rt - reaction time
%   resp - which key is pressed

% check input arguments
if acceptResp==0 && nargin<10 
    responseKeys = [];
    keyboardIndices = [];
elseif acceptResp==1 && nargin<10
    error('input arguments are missing or acceptResp is wrong')
elseif duration==Inf && acceptResp==0
    % when the stimulus is remained on the screen until a response is made,
    % the goAfterResp value must be 1 (i.e., move to the next trial once a
    % response is made)
    error('Stimulus presentation length cannot be Inf when no response is required.')
end

barLength = 16; % in pixels
barWidth = 2; % in pixels
% barColor = 255; % number from 0 (black) to 1 (white)
Screen('FillRect', window, barColor,[ (W-barLength)/2 (H-barWidth)/2 (W+barLength)/2 (H+barWidth)/2]);
Screen('FillRect', window, barColor ,[ (W-barWidth)/2 (H-barLength)/2 (W+barWidth)/2 (H+barLength)/2]);
vblCross = Screen('Flip', window); 
crossStart = vblCross; % starting time point of the fixation cross presentation

nFrames = round(duration / ifi); % number of frames based on duration of fixation cross presentation
for frame = 1:nFrames - 1 % minus one b/c the first frame is flipped outside the loop
    Screen('FillRect', window, barColor,[ (W-barLength)/2 (H-barWidth)/2 (W+barLength)/2 (H+barWidth)/2]);
    Screen('FillRect', window, barColor ,[ (W-barWidth)/2 (H-barLength)/2 (W+barWidth)/2 (H+barLength)/2]);
    vblCross = Screen('Flip', window, vblCross + (1 - 0.5) * ifi, 0); % keep updating frames
    
    if acceptResp % if a response is required to be made here 
        % record reaction time and response keys (default rt and resp are
        % zeros)
        [rt,resp] = recordResp(crossStart,responseKeys,keyboardIndices);  
        if rt>0 && ~isempty(resp) % when a response is made   
            KbWait(keyboardIndices, 1);
            acceptResp = 0; % change to 0 so that the recorded rt and resp won't be overwritten

            if goAfterResp 
                % move on to the next trial once a resp is made
                break % break the for loop (stimulus disappears and move on)
            else
                % Don't break the loop because the stimulus should remain
                % on the screen for the same duration even a response is made
                % so that we can have equal length of EEG data for each
                % stimulus.
            end
            
%             % Flash the screen so subjects know a response was made
%             Screen(window, 'FillRect', BlackIndex(max(Screen('Screens')))); % Screen(window, 'FillRect', backgroundColor);
%             Screen('Flip', window, vblCross+(1-0.5)*ifi);
        end
    end   
end

if acceptResp==0 && ~exist('rt','var') && ~exist('resp','var') % not accept response initially
    rt = 0;
    resp = [];
% elseif acceptResp==1 && rt==0 && isempty(resp) % not respond
% elseif acceptResp==0 && rt>0 && ~isempty(resp) % a response is made
end

crossEnd = vblCross; % ending time point of the fixation cross presentation

end

%%
%-------------------------------------------------------------------------%
% Draw a 7 point scale for similarity rating
%-------------------------------------------------------------------------%
function drawScale(window,Wstart,Hstart)
    HlineLength = 200; %length of horizontal line
    VlineLength = 10; %length of vertical line
    distance = 60;
    lineWidth = 3;
    lineColor = 255;
    textColor = 255;
    HLines = [-HlineLength 0; HlineLength 0]'; %[x1 y1; x2 y2]
    VLines = [-3*distance VlineLength; -3*distance -VlineLength;
        -2*distance VlineLength; -2*distance -VlineLength; 
        -distance VlineLength; -distance -VlineLength;
        0 VlineLength; 0 -VlineLength;
        distance VlineLength; distance -VlineLength;
        2*distance VlineLength; 2*distance -VlineLength;
        3*distance VlineLength; 3*distance -VlineLength]';
    %position of texts on the right
    Rx1 = Wstart+2.7*distance; 
    Ry1 = Hstart-60;
    Rx2 = Wstart+2.5*distance;
    Ry2 = Hstart-40;
    %position of texts on the left
    Lx1 = Wstart-3.3*distance;
    Ly1 = Hstart-60;
    Lx2 = Wstart-3.5*distance;
    Ly2 = Hstart-40;
    %position of each label
    Px1 = Wstart-3.1*distance;
    Py1 = Hstart+25;
    Px2 = Wstart-2.1*distance;
    Py2 = Hstart+25;
    Px3 = Wstart-1.1*distance;
    Py3 = Hstart+25;
    Px4 = Wstart-0.1*distance;
    Py4 = Hstart+25;
    Px5 = Wstart+0.9*distance;
    Py5 = Hstart+25;
    Px6 = Wstart+1.9*distance;
    Py6 = Hstart+25;
    Px7 = Wstart+2.9*distance;
    Py7 = Hstart+25;
    
    Screen('DrawLines', window, HLines, lineWidth, lineColor, [Wstart Hstart]);
    Screen('DrawLines', window, VLines, lineWidth, lineColor, [Wstart Hstart]);
%     Screen('TextFont',window, 'Arial');
%     Screen('TextSize',window, 18);
    Screen('DrawText',window,'very',Rx1,Ry1,textColor);
    Screen('DrawText',window,'similar',Rx2,Ry2,textColor);
    Screen('DrawText',window,'very',Lx1,Ly1,textColor);
    Screen('DrawText',window,'different',Lx2,Ly2,textColor);
    Screen('DrawText',window,'1',Px1,Py1,textColor);
    Screen('DrawText',window,'2',Px2,Py2,textColor);
    Screen('DrawText',window,'3',Px3,Py3,textColor);
    Screen('DrawText',window,'4',Px4,Py4,textColor);
    Screen('DrawText',window,'5',Px5,Py5,textColor);
    Screen('DrawText',window,'6',Px6,Py6,textColor);
    Screen('DrawText',window,'7',Px7,Py7,textColor);
    %Screen('Flip', window);
end

%%
%-------------------------------------------------------------------------%
% Present two stimulus images side by side and a 7-point scale for
% similarity rating
%-------------------------------------------------------------------------%
function [imgTwoStart, imgTwoEnd, rt, resp] = showStimAndScale(window, img1, img2,...
    imgPos1, imgPos2, ifi, duration, acceptResp, ...
    responseKeys, keyboardIndices, goAfterResp, screenInfo)
% input
%   window - window of the experiment
%   img1 - the matrix of image1 (e.g., target)
%   img2 - the matrix of image2 (e.g., foil)
%   img1Pos - where to present the image1 on the screen
%   img2Pos - where to present the image2 on the screen
%   ifi - interframe interval (1/refreshing rate)
%   duration - length of the image presentation
%   acceptResp - whether a response is required to be made (=1) or not (=0)
%   responseKeys - designated response keys
%   keyboardIndices - index of the input keyboard 
%   goAfterResp - whether move on to the next trial once a response is made
%                 (1-yes; 0-no)
%   screenInfo - Psychtoolbox screen information 
% output
%   imgStart - starting time point of image presentation
%   imgEnd -  ending time point of image presentation
%   rt - reaction time
%   resp - which key is pressed

% check input arguments
if acceptResp==0 && nargin<12
    responseKeys = [];
    keyboardIndices = [];
    goAfterResp = 0;
elseif acceptResp==1 && nargin<12
    error('input arguments are missing when acceptResp is set to 1')
elseif duration==Inf && acceptResp==0
    % when the stimulus is remained on the screen until a response is made,
    % the goAfterResp value must be 1 (i.e., move to the next trial once a
    % response is made)
    error('Stimulus presentation length cannot be Inf when no response is required.')
end

% Screen width and height
W = screenInfo.W;
H = screenInfo.H;

% draw image 
imageTexture1 = Screen('MakeTexture', window, img1);
imageTexture2 = Screen('MakeTexture', window, img2);
Screen('DrawTexture',window, imageTexture1, [], imgPos1);
Screen('DrawTexture',window, imageTexture2, [], imgPos2);
% draw a 7-point scale with images together
drawScale(window,W/2,H*3/4);

vblImg = Screen('Flip', window); 
imgTwoStart = vblImg; % starting time point of image presentation

if ~isinf(duration) % when the duration of stimulus presentation is not Inf
    nFrames = round(duration / ifi); % number of frames based on duration of image presentation
    for frame = 1:nFrames - 1 % minus one b/c the first frame is flipped outside the loop
        Screen('DrawTexture',window, imageTexture1, [], imgPos1);
        Screen('DrawTexture',window, imageTexture2, [], imgPos2);
        drawScale(window,W/2,H*3/4);

        vblImg = Screen('Flip', window, vblImg + (1 - 0.5) * ifi, 0);
        
        if acceptResp % if a response is required to be made here
            % record reaction time and response keys (default rt and resp are
            % zeros)
            [rt,resp] = recordResp(imgTwoStart,responseKeys,keyboardIndices);
            if rt>0 && ~isempty(resp) % when a response is made
                KbWait(keyboardIndices, 1);
                acceptResp = 0; % change to 0 so that the recorded rt and resp won't be overwritten
                
                if goAfterResp
                    % move on to the next trial once a resp is made
                    break % break the for loop (stimulus disappears and move on)
                else
                    % Don't break the loop because the stimulus should remain
                    % on the screen for the same duration even a response is made
                    % so that we can have equal length of EEG data for each
                    % stimulus.
                end
                
                %             % Flash the screen so subjects know a response was made
                %             Screen(window, 'FillRect', BlackIndex(max(Screen('Screens')))); % Screen(window, 'FillRect', backgroundColor);
                %             Screen('Flip', window, vblImg+(1-0.5)*ifi);
            end
        end
    end
    
    if acceptResp==0 && ~exist('rt','var') && ~exist('resp','var') % not accept response initially
        rt = 0;
        resp = [];
        % elseif acceptResp==1 && rt==0 && isempty(resp) % not respond
        % elseif acceptResp==0 && rt>0 && ~isempty(resp) % a response is made
    end
    
    imgTwoEnd = vblImg; % ending time point of image presentation
    
else % when the duration of stimulus presentation is Inf (wait until response is made)
    rt = 0;
    resp = [];
    while rt==0 || isempty(resp) % wait until a response is made
        [rt,resp] = recordResp(imgTwoStart,responseKeys,keyboardIndices);
    end
    imgTwoEnd = GetSecs;
    %             % Flash the screen so subjects know a response was made
    %             Screen(window, 'FillRect', BlackIndex(max(Screen('Screens')))); % Screen(window, 'FillRect', backgroundColor);
    %             Screen('Flip', window, vblImg+(1-0.5)*ifi);
end

% % Not move to next trial until a response is made 
% while rt==0 && isempty(resp)
%     [rt,resp] = recordResp(imgTwoStart,responseKeys,keyboardIndices); 
% end

% Clear textures
Screen(imageTexture1,'Close');
Screen(imageTexture2,'Close');

end