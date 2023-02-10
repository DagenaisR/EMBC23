%-------------------------------------------------------------------------%
% Script: ABP_sim_EMBC23.m
% Author: Remi Dagenais
% Date:   07/02/2023
% Descr:  Train the stacked model and then predict MAP fluctuations. 
%         Fill the first section and then run the full code at once.
% Ref:    Dagenais R., Mitsis G. D., Non-invasive estimation of arterial 
%         blood pressure fluctuations using a peripheral photoplethysmograph
%         inside the MRI scanner. EMBC23
%-------------------------------------------------------------------------%
%% Fill this section with the proper parameters
clear variables; clc; clf;

pdir = pwd;

%Inputs
subject = 0000;
ses = 1;
calib_sequence = {'preScan','postScan'}; %File used to train model
test_sequence = {'rest'}; %File use for model-predicted MAP fluct.
%filtering parameters (Hz)
LP_F = 0.15;
HP_F = 0.008;
%Lag to include
max_lag = 2; %Number of previous heartbeats to consider in model

mem = 20; %Number of continuous block for bootstrapping

pls_comp_scsa = 10; %PLS componont for eigenvalue dPLS base model
pls_comp_amp  = 8; %PLS components for amplitude dPLS base model

alpha = 0.005; %L1-L2 ratio normalization for elasticNet (alpha)

cd C:\Working\EMBC23\data %Folder containing modeling files
ndir = pwd;
listing = dir(ndir);

%% Extract Data
%Calibration data
cal = zeros(length(calib_sequence),1);
for q = 1:length(calib_sequence)
    calib_file = sprintf('SUB%i_SES0%i_%s_modeling.mat',subject,ses,calib_sequence{q});
    for w = 1:length(listing)
        if strcmpi(listing(w).name,calib_file) == 1
            load(listing(w).name);
            cal(q) = 1;

            if q == 1
                
                X1_features = normalize(predictors(:,1:end));
                Y1_features = normalize(output);
               
            else
                X2_features = normalize(predictors(:,1:end));
                Y2_features = normalize(output);
            end
        end
    end
    if cal(q) == 0
        fprintf('The calibration file %s could not be found.\n',calib_file);
    end
end

%Validation data
val = zeros(length(test_sequence),1);
for q = 1:length(test_sequence)
    test_file = sprintf('SUB%i_SES0%i_%s_modeling.mat',subject,ses,test_sequence{q});
     for w = 1:length(listing)
        if strcmpi(listing(w).name,test_file) == 1
            load(listing(w).name);
            cal(q) = 1;
            if q == 1
                X1_features_val = normalize(predictors(:,1:end));
				BP_true = output;
            else
                fprintf("Do not input more than one test sequence for the simulation.\n");
            end
        end
    end
    if cal(q) == 0
        fprintf('The validation file %s could not be found.\n',calib_file);
    end
end

%Concatenate training data
if exist("X2_features","var")
    X = [X1_features; X2_features]; 
    Y = [Y1_features; Y2_features];
else
    X = X1_features;
    Y = Y1_features;
end

%% Filtering and Outlier removal on data
[LPnum,LPdenum] = butter(2,LP_F/2,'low');
[HPnum,HPdenum] = butter(2,HP_F/2,'high');

X_new = filloutliers(X,'nearest'); X_new(isnan(X_new))= 0;
Y_new = filloutliers(Y,'nearest');
BP_true_new = filloutliers(BP_true,'nearest');

X_filt = filtfilt(LPnum,LPdenum,X_new);
X_filt = filtfilt(HPnum,HPdenum,X_filt); X = normalize(X_filt);

Y_filt = filtfilt(LPnum,LPdenum,Y_new);
Y_filt = filtfilt(HPnum,HPdenum,Y_filt); Y = normalize(Y_filt);

BP_true_filt = filtfilt(LPnum,LPdenum,BP_true_new);
BP_true_filt = filtfilt(HPnum,HPdenum,BP_true_filt); BP_true_filt = normalize(BP_true_filt);

X1_features_val = filloutliers(X1_features_val,'nearest'); X1_features_val(isnan(X1_features_val))=0;
X1_features_val = filtfilt(LPnum,LPdenum,X1_features_val);
X1_features_val = normalize(filtfilt(HPnum,HPdenum,X1_features_val));


cd(pdir) 

%% Include lag in regression matrix
if ~exist("prediction","var")
    prediction = 1; %MAP
end

Y = Y(:,prediction);

X_lag = zeros(size(X,1),size(X,2)*max_lag);
X1_lag = zeros(size(X1_features_val,1),size(X,2)*max_lag);

count = 1;
for w=1:size(X,2)
    for q = 1:max_lag
        X_lag(:,count) = circshift(X(:,w),(q-1)*4);
        X1_lag(:,count) = circshift(X1_features_val(:,w),(q-1)*4);
        count = count+1;
    end
end

X = X_lag;
X1_features_val = X1_lag;


%% Calibrate model

%Apply bootstrapping for calibration
BS_bloc = ceil(size(X,1)/mem);

X_train = zeros(mem,size(X,2));
Y_train = zeros(mem,1);
X_test  = zeros(mem,size(X,2));
Y_test  = zeros(mem,1);
VAF_test_pred = zeros(20,1);
Y1_features_pred = zeros(size(X1_features_val,1),20);

B = zeros(size(X1_features,2)*max_lag+1,20);
for cv = 1:20
    idx_BS = randi(BS_bloc,BS_bloc,1);

    X_cat = [X;X]; Y_cat = [Y;Y];
    for q = 1:BS_bloc
        X_train = [X_train; X_cat(idx_BS(q)*mem:mem*(idx_BS(q)+1)-1,:)];
        Y_train = [Y_train; Y_cat(idx_BS(q)*mem:mem*(idx_BS(q)+1)-1,:)];
    end
    X_train = X_train(mem+1:end,:);
    Y_train = Y_train(mem+1:end);
    
    full_idx = 1:BS_bloc;

    %Test on out of sample data
    out_of_sample_idx = ismember(full_idx,idx_BS); 
    idx_test = full_idx(out_of_sample_idx);
    for q = 1:length(idx_test)
        X_test = [X_test; X_cat(idx_test*mem:mem*(idx_test+1)-1,:)];
        Y_test = [Y_test; Y_cat(idx_test*mem:mem*(idx_test+1)-1)];
    end
    X_test = X_test(mem+1:end,:);
    Y_test = Y_test(mem+1:end);

    %Train base models
    [~,~,~,~,betaSCSA,~,~,~] = plsregress(X_train(:,1:19*max_lag),Y_train,pls_comp_scsa);
    BSCSA(:,cv) = betaSCSA(:); XSCSA = [ones(size(X_train,1),1) X_train(:,1:19*max_lag)]*betaSCSA;
    
    [~,~,~,~,betaAMP,~,~,~] = plsregress(X_train(:,(19*max_lag)+1:(22*max_lag)+(19*max_lag)),Y_train,pls_comp_amp);
    BAMP(:,cv) = betaAMP(:); XAMP = [ones(size(X_train,1),1) X_train(:,(19*max_lag)+1:(22*max_lag)+(19*max_lag))]*betaAMP;
    
    %Stack model
    [BetaRidge,S] = lasso([XSCSA XAMP X_train(:,(22*max_lag)+(19*max_lag)+1:end)],Y_train,'alpha',K);
    idx_b = find(S.MSE == min(S.MSE)); BetaRidge = BetaRidge(:,idx_b);
    
    %Prediction
    y_pred_test = normalize([[ones(size(X_test,1),1) X_test(:,1:19*max_lag)]*betaSCSA [ones(size(X_test,1),1) X_test(:,(19*max_lag)+1:(22*max_lag)+(19*max_lag))]*betaAMP X_test(:,(22*max_lag)+(19*max_lag)+1:end)],1)*BetaRidge; 
    y_pred_val  = normalize([[ones(size(X1_features_val,1),1) X1_features_val(:,1:19*max_lag)]*betaSCSA [ones(size(X1_features_val,1),1) X1_features_val(:,(19*max_lag)+1:(22*max_lag)+(19*max_lag))]*betaAMP X1_features_val(:,(22*max_lag)+(19*max_lag)+1:end)],1)*BetaRidge; 
    

    VAF_test_pred(cv) = 1-var(y_pred_test-Y_test)/var(Y_test);
    Y1_features_pred(:,cv) = y_pred_val;

    %Re-initialise values
    X_train = zeros(mem,size(X,2));
    Y_train = zeros(mem,1);
    X_test  = zeros(mem,size(X,2));
    Y_test  = zeros(mem,1);
end

%% Visualize result & Interpolate beginning and end
BP_sim = mean(Y1_features_pred,2); 
time = [0:0.25:0.25*(size(BP_true_filt,1)-1)]';

figure(1)
plot(time,normalize(BP_true_filt),'b-','linewidth',2)
hold on
plot(time,normalize(BP_sim),'r-','linewidth',2);
xlabel('Time (s)','fontweight','bold','fontsize',12);
ylabel('Blood Pressure (a.u.)','fontweight','bold','fontsize',12);
legend('BP measured','BP simulation','fontweight','bold');
set(gca,'box','off');

figure(2)
mscohere(BP_true_filt,BP_sim,[],[],[],4);
xlim([0 0.2]);

for q = 1:cv
    R2(q) = 1-var(normalize(BP_true_filt(200:end))-normalize(Y1_features_pred(200:end,q)))/var(normalize(BP_true_filt(200:end)));
    r(q) = corr(BP_true_filt(200:end),Y1_features_pred(200:end,q));
end

R2_mean = mean(R2)
R2_std = std(R2)

r_mean = mean(r)
r_std = std(r)