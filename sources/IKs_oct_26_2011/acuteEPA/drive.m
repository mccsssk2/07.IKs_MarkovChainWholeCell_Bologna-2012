% Inputs to this program are:
% CaseCase: see below for explanation.
% The voltage clamp protocol has be to accurately set. If lower voltages are changed, it should work.
% if the max voltage of 40 mV is changed, then the normalisation by max(data(:,10)) needs changing appropriately.
%  
% The experimental data has to be in a dir called exptData to avoid errenous deletion.

clear all;
%close all;

global filecounter;
global CaseCase;
global CostCost;
global printOutput;
global useSingle;

global tracesIV;
global trace5p5;
global trace12s;
global defGKs;
global fixedGKs;
global weightIV;
global weightDataIV;

global vcc;
global FoRT;
global erev;

% if print output and ploting or not. 0 for not, 1 for yes

printOutput = 0;

% CaseCase = 0 means control
% CaseCase = 1 means acute DHA
% CaseCase = 2 means acute EPA
% CaseCase = 3 means chronic EPA

CaseCase = 2;

% initGuess = 0 restarting from previous guess
% initGuess =  1 use Silva and Rudy parameters
% initGuess = 2 use Abraham parameters
% initGues = 3 use Abraham variant parameters

initGuess = 0;

% CostCost = 0 no weighting
% CostCost = 1 weighted cost function - weighted IV traces, weighted time
% points (see below), weighted tail currents
% CostCost = 2 weighted by derivative, weighted IV traces, tail currents

CostCost = 0;

% useSingle = 0 dont use single traces
% useSingle = 1 fit also the single traces

useSingle = 0;

% fixedGKs = 0 dont fix GKs
% fixedGKs = 1 fix GKs to the last value

fixedGKs = 1;

% DataType = 0 use Median traces
% DataType = 1 use Average traces

DataType = 1;


filecounter = 1;

% define the voltage clamp protocol from Restier data.
vcc(1) = -80.0; % vhold.
vcc(2) = 0.0; % lowest value of voltage clamp. CHANGE THIS MOSTLY.
vcc(3) = 10.0; % voltage step interval, so then the clamp goes -60, -40, ... This has to be a constant as of now.
vcc(4) = 60.0; %  highest value of voltage clamp. IF YOU CHANGE THIS, MAKE ALTERATIONS TO HIKSIV.M
vcc(5) = -40.0; % voltage when it comes back for the deactivation.

% some constants.
F=96485;
T=295;
R=8314;
FoRT=F/R/T;
pnak=.01833;
k_o=5.4;          % mM
na_o=136.0;       % mM
na_i=10;        % mM
k_i=141.1;        % mM
erev=(1/FoRT)*log((k_o+pnak*na_o)/(k_i+pnak*na_i));  % Permeability Ratio

if DataType == 0
if CaseCase == 0
    pathExp = sprintf('../exptData/CtrlMed.txt');
    pathExp5p5 = sprintf('../exptData/Ctrl5p5Med.txt');
    pathExp12 = sprintf('../exptData/Ctrl12Med.txt');
end
if CaseCase == 1
    pathExp = sprintf('../exptData/AcuteDHAMed.txt');
    pathExp5p5 = sprintf('../exptData/AcuteDHA5p5Med.txt');
    pathExp12 = sprintf('../exptData/AcuteDHA12Med.txt');
end
if CaseCase == 2
    pathExp = sprintf('../exptData/AcuteMed.txt');
    pathExp5p5 = sprintf('../exptData/Acute5p5Med.txt');
    pathExp12 = sprintf('../exptData/Acute12Med.txt');
end
if CaseCase == 3
    pathExp = sprintf('../exptData/ChronicMed.txt');
    pathExp5p5 = sprintf('../exptData/Chronic5p5Med.txt');
end
end

if DataType == 1
    if CaseCase == 0
        pathExp = sprintf('../exptData/CtrlAvg.txt');
        pathExp5p5 = sprintf('../exptData/Ctrl5p5Avg.txt');
        pathExp12 = sprintf('../exptData/Ctrl12Avg.txt');
    end
    if CaseCase == 1
        pathExp = sprintf('../exptData/AcuteDHAAvg.txt');
        pathExp5p5 = sprintf('../exptData/AcuteDHA5p5Avg.txt');
        pathExp12 = sprintf('../exptData/AcuteDHA12Avg.txt');
    end
    if CaseCase == 2
        pathExp = sprintf('../exptData/AcuteAvg.txt');
        pathExp5p5 = sprintf('../exptData/Acute5p5Avg.txt');
        pathExp12 = sprintf('../exptData/Acute12Avg.txt');
    end
    if CaseCase == 3
        pathExp = sprintf('../exptData/ChronicAvg.txt');
        pathExp5p5 = sprintf('../exptData/Chronic5p5Avg.txt');
    end    
end

% load experimental data
tracesIV = load(pathExp); 
% if useSingle == 1
    maximum=max(max(tracesIV(:,2:end)));
    trace5p5 = load(pathExp5p5);
    trace5p5(:,2) =  trace5p5(:,2)*maximum/max(trace5p5(:,2));
    if CaseCase ~= 3
    activationstep = 21978;
    trace12s = load(pathExp12);
    trace12s(:,2) =  trace12s(:,2)*maximum/max(trace12s(activationstep,2));
    end
% end


% WEIGHTING SETTINGS
% if CostCost
    % weight of IV traces 
    %         60 50 40 30 20 10  0 mV trace 
    weightIV=[2:-.1:.2];
    %weightIV=[.5 .5 .5 1 1.5 2 1.5];

    weightDataIV=zeros(length(tracesIV),1);
    activationstep = 21979;

for i=1:length(tracesIV)
    % weighting of data point
    if tracesIV(i,1)<600 %time in ms
        weightDataIV(i)=5; 
    elseif i>activationstep
        weightDataIV(i)=10;
    else
        weightDataIV(i)=1;
    end
end

% end


pariks = [];

% first estimate for the Restier WT data from an unconstrained fminsearch
% estimation. Note that no constraints or bounds were applied. Further, the
% cost function did not go below 1.3. The estimation was done in 3 steps,
% first an estimate based on trace at 40 mV was obtained. The estimate was
% used as starting point for estimating parameter vector by simulataneously
% using traces at 0,20, and 40 mV. Finally, the complete Restier slowpulse
% dataset for WT was used. Some intelligent optimisation will improve the
% estimate. Notice how there have been a few zero crossings.
% 3.2742598e-04  4.2745833e-01   1.7852655e-04   1.7622375e-02   6.2013976e-03   1.2892320e+00   2.5430853e-03   
%5.7081376e-01   2.1436489e-01   2.2353195e-01   1.4570502e-01  -6.4107196e-04   5.2795317e-01   6.9488933e-03  -4.2281252e-02


% estimate of the V141M data. This is how far it got.
% -1.6352481e-05  -2.2141027e-02   6.0084918e-06  -3.4226191e-01   2.7998983e-03   6.1438896e-01   7.5460757e-04   7.5914014e-03   2.6885823e-03   3.1396744e-04  -1.9140012e-01   8.4152794e-02   1.1634038e+00   2.4976813e-06  -1.1143678e+00

% best estimate for S140G so far.
%    3.5009739e-04   4.6224623e-02   1.1740509e-04   1.6503290e-02   1.0835359e-03   5.8053955e-01   4.2653523e-04  -1.4973130e-02   2.2429311e-01   3.2470147e-02   7.1023454e-01  -2.4634911e-04  -1.7168513e-01   6.8016568e-04   1.4403641e-04
if(initGuess==0)
pariks = load('restart.dat');
end
if(initGuess==1) %SR
pariks(1) = 3.98e-4; % multiplier of alpha POSITIVE
pariks(2) = 0.361; % exponential term of alpha 
pariks(3) = 5.74e-5; % multiplier of beta POSITIVE
pariks(4) = -9.23e-2; % exp. of beta
pariks(5) = 3.41e-3;  % multiplier of gamma POSITIVE
pariks(6) = 8.68e-1; % exponential of gamma
pariks(7) = 1.20e-3; % multiplier delta POSITIVE
pariks(8) = -3.30e-1; % exp. of delta
pariks(9) = 6.47e-3; % theta = has to be positive
pariks(10) = 1.25e-2; % eta multi. POSITIVE
pariks(11) = -4.81e-1; % eta exp.
pariks(12) = 6.33e-3; % psi multi. POSITIVE
pariks(13) = 1.27; % psi exp.
pariks(14) = 4.91e-3; % omega multiplier POSITIVE
pariks(15) = -6.79e-1; % omega exp.
pariks(16)=79.8; % GKs [pS] % 75.8
end
if initGuess==2 % Abraham
pariks(1) = 6.4038e-4; % multiplier of alpha POSITIVE
pariks(2) = 2.3279e-1; % exponential term of alpha 
pariks(3) = 2.0242e-4; % multiplier of beta POSITIVE
pariks(4) = -4.1445e-1; % exp. of beta
pariks(5) = 1.5171e-3;  % multiplier of gamma POSITIVE
pariks(6) = 6.9472e-1; % exponential of gamma
pariks(7) = 6.4186e-4; % multiplier delta POSITIVE
pariks(8) = 3.1183e-1; % exp. of delta
pariks(9) = 2.3767e-2; % theta = has to be positive
pariks(10) = 3.4970e-2; % eta multi. POSITIVE
pariks(11) = 1.5836e-1; % eta exp.
pariks(12) = 2.1447e-4; % psi multi. POSITIVE
pariks(13) = -3.5853e-2; % psi exp.
pariks(14) = 9.3089e-9; % omega multiplier POSITIVE
pariks(15) = -7.9873; % omega exp.
pariks(16)=79.8; % GKs [pS] % 75.8
end
if initGuess==3 % Abraham variant
pariks(1) = 1.9568*6.4038e-4; % multiplier of alpha POSITIVE
pariks(2) = 1.3573*2.3279e-1; % exponential term of alpha 
pariks(3) = 2.0094*2.0242e-4; % multiplier of beta POSITIVE
pariks(4) = .47032*(-4.1445e-1); % exp. of beta
pariks(5) = 1.0656*1.5171e-3;  % multiplier of gamma POSITIVE
pariks(6) = 1.4075*6.9472e-1; % exponential of gamma
pariks(7) = 4.8069*6.4186e-4; % multiplier delta POSITIVE
pariks(8) = -5.0756e-6*3.1183e-1; % exp. of delta
pariks(9) = .54262*2.3767e-2; % theta = has to be positive
pariks(10) = .39850*3.4970e-2; % eta multi. POSITIVE
pariks(11) = 5.8032*1.5836e-1; % eta exp.
pariks(12) = 2432.3*2.1447e-4; % psi multi. POSITIVE
pariks(13) = 3.1419e-6*(-3.5853e-2); % psi exp.
pariks(14) = 4.1363e7*9.3089e-9; % omega multiplier POSITIVE
pariks(15) = 1.3447*(-7.9873); % omega exp.
end

if printOutput == 1
    [ans,ans]=mkdir('png');
end

defGKs=pariks(16);


% so this cost function works. I think we are set for the fminsearch.
% hiksiv(pariks, vc);


% the next step can be fmincon.
if CostCost == 0 && useSingle == 0
    display('Estimation without weighting');
    for k = vcc(4):-30:0
        vcc(2)=k;
        hlaska=sprintf('include trace %d mV',k);
        display(char(hlaska))
        options=optimset('FunValCheck','on','TolFun',0.9,'Display','off','MaxIter',1000);
        [parameters, fval, exitflag, outputs] = fminsearch('hiksiv',pariks,options);
        filename=sprintf('pariks_noWeight_from_%dmV.dat',k);
        pariks=parameters;
        save(filename,'pariks','-ASCII');
    end
    CostCost = 1;
end

if CostCost == 1 && useSingle == 0
     display('The weighting is used')
    options=optimset('FunValCheck','on','TolFun',0.9,'Display','off','MaxIter',1000);
    [parameters, fval, exitflag, outputs] = fminsearch('hiksiv',pariks,options);
    pariks=parameters;
    save('pariks_weighttail.dat','pariks','-ASCII');
    useSingle = 1;
    display('Using single traces')
end


options=optimset('FunValCheck','on','TolFun',0.9,'Display','off','MaxIter',1700);
[parameters, fval, exitflag, outputs] = fminsearch('hiksiv',pariks,options);
pariks=parameters;
save('pariks_final1.dat','pariks','-ASCII');
[parameters, fval, exitflag, outputs] = fminsearch('hiksiv',pariks,options);
pariks=parameters;
save('pariks_final2.dat','pariks','-ASCII');
[parameters, fval, exitflag, outputs] = fminsearch('hiksiv',pariks,options);
pariks=parameters;
save('pariks_final3.dat','pariks','-ASCII');
