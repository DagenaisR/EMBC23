%-------------------------------------------------------------------------%
% Script: preModeling_EMBC23.m
% Author: Remi Dagenais
% Date:   07/02/2023
% Descr:  Extract features from preProcessed PPG. 
%         Run one section at a time and follow the instructions. 
% Ref:    Dagenais R., Mitsis G. D., Non-invasive estimation of arterial 
%         blood pressure fluctuations using a peripheral photoplethysmograph
%         inside the MRI scanner. EMBC23
%-------------------------------------------------------------------------%
%% Open data
clear variables; clc; clf;

filename = 'SUB0000_SES01_rest_processed.mat';
load(filename);

%% Find pressure values for ABP
[Sys,Mean,Dia] = extractPressure(abp,time,pp_idx,3*Fs/5);

%% MAP resampled at 4Hz
% ppMEAN
[MAP,tmap] = resample(Mean.pres,Mean.time,4,1,1);
%% get Low passed PPG
[LPnum,LPdenum] = butter(2,0.008/Fs);

% forward-backward filter to cancel delay
ppgLP = filtfilt(LPnum,LPdenum,ppg_filt);
% Down-sample at 4Hz both signals
[ppgLP,~] = resample(ppgLP,time,4,1,1);
ppg_filt = normalize(ppg_filt);

%% Extract PPG pulses

[P1,Pm,P0] = extractPressure(ppg_filt,time,p1_idx,3*Fs/5);
idx_P0 = zeros(length(P0.time),1);
for q = 1:length(P0.time(:))
    idx_P0(q) = find(time == P0.time(q));
end
% Longest pulse
size_pulse = max(diff(idx_P0));

pulse = zeros(length(idx_P0)-1,size_pulse);
min_pulse = zeros(length(idx_P0)-1,1);
max_pulse = zeros(length(idx_P0)-1,1);
length_pulse = diff(idx_P0);

for q = 1:length(idx_P0)-1 %Separate each pulses
   pulse(q,1:(idx_P0(q+1)-idx_P0(q))+1) = ppg_filt(idx_P0(q):idx_P0(q+1));
   pulse(q,(idx_P0(q+1)-idx_P0(q))+2:end) = pulse(q,(idx_P0(q+1)-idx_P0(q))+1);
   min_pulse(q) = min(ppg_filt(idx_P0(q):idx_P0(q+1))); %minimum envelope
   max_pulse(q) = max(ppg_filt(idx_P0(q):idx_P0(q+1))); %maximum envelope 
   pulse(q,:) = pulse(q,:)-min_pulse(q);
end

%% SCSA decomposition
[yh,res,K,phi,Beta,X] = scsa(pulse,length_pulse,10); %semi-classical signal analysis decomposition
INVS1 = (4./X).*sum(K(:,1:3),2); %First systolic invariant
INVS2 = (16./(3*X)).*sum(K(:,1:3).^3,2); %Second systolic invariant
INVD1 = (4./X).*sum(K(:,4:10),2); %First diastolic invariant
INVD2 = (16./(3*X)).*sum(K(:,4:10).^3,2); %Second systolic invariant
Kt    = (4./X).*sum(K(:,1:10),2); %Eigenvalue sum

%% Interpolate Features at 4Hz
timeFeatures = time(p1_idx(2:end-1));
timeP3 = time(p3_idx(2:end-1));
HR = diff(time(p1_idx));
min_pulse_resampled = zeros(size(tmap,1),1);
max_pulse_resampled = zeros(size(tmap,1),1);
K_resampled = zeros(size(tmap,1),size(K,2));
INVS1_resampled = zeros(size(tmap,1),1);
INVS2_resampled = zeros(size(tmap,1),1);
INVD1_resampled = zeros(size(tmap,1),1);
INVD2_resampled = zeros(size(tmap,1),1);
Kt_resampled    = zeros(size(tmap,1),1);
Beta_resampled  = zeros(size(tmap,1),2);
res_resampled   = zeros(size(tmap,1),1);
HR_resampled    = zeros(size(tmap,1),1);
P1_resampled    = zeros(size(tmap,1),1);
P3_resampled    = zeros(size(tmap,1),1);
Pm_resampled    = zeros(size(tmap,1),1);
P0_resampled    = zeros(size(tmap,1),1);
T1Tm_resampled  = zeros(size(tmap,1),1);
TmT0_resampled  = zeros(size(tmap,1),1);
T1T0_resampled  = zeros(size(tmap,1),1);
amp             = zeros(size(tmap,1),18);
temp            = zeros(size(tmap,1),4);
for q = 1:size(K,2)  
   K_resampled(:,q) = interp1(timeFeatures,K(:,q),tmap,'linear','extrap'); 
   if q < 3
      Beta_resampled(:,q) = interp1(timeFeatures,Beta(:,q),tmap,'linear','extrap'); 
   end
   
   if q < 2
      min_pulse_resampled = interp1(time(idx_P0(1:end-1)),min_pulse,tmap,'linear','extrap');
      max_pulse_resampled = interp1(time(idx_P0(1:end-1)),max_pulse,tmap,'linear','extrap');
      res_resampled   = interp1(timeFeatures,res,tmap,'linear','extrap'); 
      HR_resampled    = interp1(time(p1_idx(1:end-1)),HR,tmap,'linear','extrap');
      P1_resampled    = interp1(P1.time,P1.pres,tmap,'linear','extrap');
      P3_resampled    = interp1(timeP3,ppg_filt(p3_idx(2:end-1)),tmap,'linear','extrap');
      Pm_resampled    = interp1(Pm.time,Pm.pres,tmap,'linear','extrap');
      P0_resampled    = interp1(P0.time,P0.pres,tmap,'linear','extrap');
      T1Tm_resampled  = interp1(P1.time,P1.time-Pm.time,tmap,'linear','extrap');
      TmT0_resampled  = interp1(Pm.time,Pm.time-P0.time,tmap,'linear','extrap');
      T1T0_resampled  = interp1(P1.time,P1.time-P0.time,tmap,'linear','extrap');
      INVS1_resampled = interp1(timeFeatures,INVS1,tmap,'linear','extrap');
      INVS2_resampled = interp1(timeFeatures,INVS2,tmap,'linear','extrap');
      INVD1_resampled = interp1(timeFeatures,INVD1,tmap,'linear','extrap');
      INVD2_resampled = interp1(timeFeatures,INVD2,tmap,'linear','extrap');
      Kt_resampled    = interp1(timeFeatures,Kt,tmap,'linear','extrap');
   end
end

% Amplitude features
% P1-P0
amp(:,1) = P1_resampled-P0_resampled;

% Pm-P0
amp(:,2) = Pm_resampled-P0_resampled;

% P3-P0
amp(:,3) = P3_resampled-P0_resampled;

% P1-P3
amp(:,4) = P1_resampled-P3_resampled;

% P1-Pm
amp(:,5) = P1_resampled-Pm_resampled;

% P3-Pm
amp(:,6) = P3_resampled-Pm_resampled;

% P0/P1
amp(:,7) = P0_resampled./P1_resampled;

% Pm/P1
amp(:,8) = Pm_resampled./P1_resampled;

% P3/P1
amp(:,9) = P3_resampled./P1_resampled;

% P0/P3
amp(:,10) = P0_resampled./P3_resampled;

% pm/p3
amp(:,11) = Pm_resampled./P3_resampled;

% P0/Pm
amp(:,12) = P0_resampled./Pm_resampled;

% P3/(P1-P0)
amp(:,13) = P3_resampled./amp(:,1);

% Pm/(P1-P0)
amp(:,14) = Pm_resampled./amp(:,1);

% Pm/(P3-P0)
amp(:,15) = Pm_resampled./amp(:,3);

% (P3-P0)/(P1-P0)
amp(:,16) = amp(:,3)./amp(:,1);

% (Pm-P0)/(P1-P0)
amp(:,17) = amp(:,2)./amp(:,1);

% (Pm-P0)/(P3-P0)
amp(:,18) = amp(:,2)./amp(:,3);

% Temporal features
% T1T0
[ppg_t1t0,tt1t0] = resample(time(p1_idx(2:end))'-P0.time,time(p1_idx(2:end)),4,1,1);
temp(:,1) = interp1(tt1t0,ppg_t1t0,tmap,'linear');

% T3T0 
[ppg_t3t0,tt3t0] = resample(time(p3_idx(2:end-1))'-P0.time(1:end-1),time(p3_idx(2:end-1)),4,1,1);
temp(:,2) = interp1(tt3t0,ppg_t3t0,tmap,'linear');

% TmT0
[ppg_tmt0,ttmt0] = resample(Pm.time-P0.time,Pm.time,4,1,1);
temp(:,3) = interp1(ttmt0,ppg_tmt0,tmap,'linear');

% T3T1
[ppg_t3t1,tt3t1] = resample(ct.time(p3_idx(2:end-1))-time(p1_idx(2:end-1)),time(p3_idx(2:end-1)),4,1,1);
temp(:,4) = interp1(tt3t1,ppg_t3t1,tmap,'linear');

%% Categorize predictors
% size(SCSA,2) = 19
% size(amplitude,2) = 22
% size(temp,2) = 7
% size(drift,2) = 3
SCSA      = [K_resampled INVS1_resampled INVS2_resampled INVD1_resampled INVD2_resampled Kt_resampled Beta_resampled res_resampled HR_resampled];
amplitude = [P1_resampled P3_resampled Pm_resampled P0_resampled amp];
temporal  = [T1Tm_resampled TmT0_resampled T1T0_resampled temp];
drift     = [min_pulse_resampled max_pulse_resampled LP];

%% Export data for modelling
predictors = [SCSA amplitude temporal drift]; 
output = MAP;

filename=sprintf('SUB%i_SES0%i_%s',sub,ses,sequence);
save([filename '_modeling.mat'],'sub','ses','sequence','predictors','output');

%% Functions
%-------------------------------------------------------------------------%
% Function: Find ABP (Sys, Mean, Dia) from BP signal and systolic idx
% Written by: Rémi Dagenais
% Date: 2021-12-03
% INPUT -> BP signal
%       -> time vector
%       -> systolic idx
%       -> window size
%       -> Resampling rate (optional)
% OUTPUT -> Sys.___
%        ->    .time
%        ->    .pres
%        -> Mean.___
%        ->     .time
%        ->     .pres
%        -> Dia.___
%        ->    .time
%              .pres
% DESCRIPTION -> Will output the pressure values along with time points
% Modification:
% 05/01/2022 -> dia_idx is now a int64 variable.
% 12/10/2022 -> outputs are column vectors 
%-------------------------------------------------------------------------%

function [sys,mea,dia] = extractPressure(BP,time,sys_idx,window,varargin)

sys_idx = int64(sys_idx);
sys_idx = sys_idx(:);
% Systole 
sys.pres = BP(sys_idx); sys.pres = sys.pres(:);
sys.time = time(sys_idx); sys.time = sys.time(:);


tmp_window = window;
% Find diastole from systole
for q = 1:length(sys_idx)
    if sys_idx(q)<window
        tmp_window = window;
        window = sys_idx(q)-1;
    end
    idx = find(min(BP(sys_idx(q)-window:sys_idx(q))) == BP(sys_idx(q)-window:sys_idx(q))); 
    dia_idx(q) = idx(end);
    clear idx
    window = tmp_window;
end

%Idx is an int vector
dia_idx = int64(dia_idx);
dia_idx = dia_idx(:);

dia.time = time(sys_idx-window+dia_idx); dia.time = dia.time(:);
dia.pres = BP(sys_idx-window+dia_idx); dia.pres = dia.pres(:);

time = time(:)'; %Force row vector
% Compute Mean
mea.time = diff([0 time(sys_idx)])./2 + [0 time(sys_idx(1:end-1))]; mea.time = mea.time(:);
mea.pres = (2*dia.pres + sys.pres)./3;
mea.pres = mea.pres(:);

if nargin>4
   [sy,sty] = resample(sys.pres,sys.time,varargin{1},1,1); sys.pres = sy(:); sys.time = sty(:);
   [my,mty] = resample(mea.pres,mea.time,varargin{1},1,1); mea.pres = my(:); mea.time = mty(:);
   [dy,dty] = resample(dia.pres,dia.time,varargin{1},1,1); dia.pres = dy(:); dia.time = dty(:);
end

% figure
hold on
plot(sys.time,sys.pres,'r-');
plot(mea.time,mea.pres,'k-');
plot(dia.time,dia.pres,'b-');
hold off
end

%-------------------------------------------------------------------------%
% Function: Semi-classical signal analysis using Schrödinger Operator
% qritten by: Rémi Dagenais
% Date: 2023-01-18
% INPUT -> y (size(N_pulse,max(length_pulse))
%       -> length_pulse (size(N_pulse,1))
%       -> Nh (Integer = Desired # of eigenvalues)
% OUTPUT -> yh (size(N_pulse,max(length_pulse))
%        -> res (size(N_pulse))
%        -> K (size(N_pulse,Nh))
%        -> phi {N_pulse}(size(max(length_pulse),Nh))
%        -> Beta (size(N_pulse,2))
%        -> X_f (size(N_pulse,1))
% DESCRIPTION -> Decomposes the signal y into a discrete spectrum of the
% Schrödinger operator: -d^2(phi)/dt^2 - X*y*phi = lambda*phi
% Returns the estimate of the function y, the residuals, the eigenvalues,
% the eigenvectors, and the scalling coefficients.
%-------------------------------------------------------------------------%
function [yh,res,K,Phi,Beta,X_f] = scsa(y,length_pulse,Nh)
% Initialize variables
K = zeros(size(y,1),Nh);
res = zeros(size(y,1),1);
Beta = zeros(size(y,1),2);

yh_size = max(length_pulse);
yh = zeros(size(y,1),yh_size);
phi{size(y,1)} = zeros(length_pulse(end),Nh);
X_f = zeros(size(y,1),1);
for q = 1:size(y,1)
    if mod(length_pulse(q),2) == 0
        length_pulse(q) = length_pulse(q);
    else
        length_pulse(q) = length_pulse(q)-1;
    end
    
    N = 0;
    flag = 0;
    while flag == 0
        N = N+1; X = N*(N+1);
        M = length_pulse(q);
        dt = 0.002;
        d2 = D2(M,dt);
        [Phi,lambda] = eig(-d2-(X*diag(y(q,1:length_pulse(q))))); %Solve eigenvalue problem
        lambda = diag(lambda);
        if lambda(Nh) < 0
            K(q,:) = (-lambda(1:Nh)).^0.5;
            yh(q,1:length_pulse(q)) = 4/X*sum(K(q,:)'.*Phi(:,1:Nh)'.^2); %Compute 
            Beta(q,:) = polyfit(yh(q,1:length_pulse(q)),y(q,1:length_pulse(q)),1);
            res(q) = sum((y(q,1:length_pulse(q))-polyval(Beta(q,:),yh(q,1:length_pulse(q)))).^2);
            phi{q} = Phi(:,1:Nh);
            flag = 1;
            X_f(q) = X;
        else
            continue
        end
    end
end

    function [d2] = D2(M,dt) %Diagonal symmetric difference operator
        %Initialise variables
        d2 = zeros(M,M);
        delta = 2*pi/M;
        if mod(M,2) == 0 %M is Even
            for k = 1:M
                for j = 1:M
                    if k==j
                        d2(k,j) = (delta^2/dt^2)*(-pi^2/(3*delta^2)-1/6);
                    else
                        d2(k,j) = (delta^2/dt^2)*(-(-1)^(k-j)*0.5*1/(sin((k-j)*delta/2))^2);
                    end
                end
            end
            
        else %M is Odd
            for k = 1:M
                for j = 1:M
                    if k==j
                        d2(k,j) = (delta^2/dt^2)*(-pi^2/(3*delta^2)-1/12);
                    else
                        d2(k,j) = (delta^2/dt^2)*(-(-1)^(k-j)*cot((k-j)*delta/2)*0.5/sin((k-j)*delta/2));
                    end
                end
            end
        end
    end

end