%% Texture Discrimination processing
% updated on 2/15/2024

clear all;
close all;
clc;
%% Load in 26 texture files

sample_size = 322000; % 40 samples * window length of 8050;
    % Based on the experimental setup

% Soft Textures

% Soft Flat
file_1 = load('Multilayer_sensor_soft_flat.mat');
sample_offset = 1000;
Texture_1 = file_1.finalData(2:19,sample_offset:sample_size+sample_offset-1);


% Soft Bumps
file_2 = load('Multilayer_sensor_soft_3_bumps.mat');
sample_offset = 8000;
Texture_2 = file_2.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_3 = load('Multilayer_sensor_soft_4_bumps.mat');
sample_offset = 4000;
Texture_3 = file_3.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_4 = load('Multilayer_sensor_soft_6_bumps.mat');
sample_offset = 1;
Texture_4 = file_4.finalData(2:19,sample_offset:sample_size+sample_offset-1);


% Soft Circular ridges
file_5 = load('Multilayer_sensor_soft_3_circular_ridges.mat');
sample_offset = 1;
Texture_5 = file_5.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_6 = load('Multilayer_sensor_soft_4_circular_ridges.mat');
sample_offset = 1;
Texture_6 = file_6.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_7 = load('Multilayer_sensor_soft_6_circular_ridges.mat');
sample_offset = 1;
Texture_7 = file_7.finalData(2:19,sample_offset:sample_size+sample_offset-1);


% Soft Triangular ridges
file_8 = load('Multilayer_sensor_soft_3_triangular_ridges.mat');
sample_offset = 1;
Texture_8 = file_8.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_9 = load('Multilayer_sensor_soft_4_triangular_ridges.mat');
sample_offset = 1;
Texture_9 = file_9.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_10 = load('Multilayer_sensor_soft_6_triangular_ridges.mat');
sample_offset = 1;
Texture_10 = file_10.finalData(2:19,sample_offset:sample_size+sample_offset-1);


% Soft Sinusoidal waves
file_11 = load('Multilayer_sensor_soft_3_sinusoidal.mat');
sample_offset = 1;
Texture_11 = file_11.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_12 = load('Multilayer_sensor_soft_4_sinusoidal.mat');
sample_offset = 1;
Texture_12 = file_12.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_13 = load('Multilayer_sensor_soft_6_sinusoidal.mat');
sample_offset = 1;
Texture_13 = file_13.finalData(2:19,sample_offset:sample_size+sample_offset-1);

% Hard Textures

% Hard Flat
file_14 = load('Multilayer_sensor_hard_flat.mat');
sample_offset = 1;
Texture_14 = file_14.finalData(2:19,sample_offset:sample_size+sample_offset-1);


% Hard Bumps
file_15 = load('Multilayer_sensor_hard_3_bumps.mat');
sample_offset = 1;
Texture_15 = file_15.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_16 = load('Multilayer_sensor_hard_4_bumps.mat');
sample_offset = 1000;
Texture_16 = file_16.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_17 = load('Multilayer_sensor_hard_6_bumps.mat');
sample_offset = 7000;
Texture_17 = file_17.finalData(2:19,sample_offset:sample_size+sample_offset-1);


% Hard Circular ridges
file_18 = load('Multilayer_sensor_hard_3_circular_ridges.mat');
sample_offset = 1;
Texture_18 = file_18.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_19 = load('Multilayer_sensor_hard_4_circular_ridges.mat');
sample_offset = 1;
Texture_19 = file_19.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_20 = load('Multilayer_sensor_hard_6_circular_ridges.mat');
sample_offset = 1;
Texture_20 = file_20.finalData(2:19,sample_offset:sample_size+sample_offset-1);


% Hard Triangular ridges
file_21 = load('Multilayer_sensor_hard_3_triangular_ridges.mat');
sample_offset = 8500;
Texture_21 = file_21.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_22 = load('Multilayer_sensor_hard_4_triangular_ridges.mat');
sample_offset = 1;
Texture_22 = file_22.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_23 = load('Multilayer_sensor_hard_6_triangular_ridges.mat');
sample_offset = 1;
Texture_23 = file_23.finalData(2:19,sample_offset:sample_size+sample_offset-1);


% Hard Sinusoidal waves
file_24 = load('Multilayer_sensor_hard_3_sinusoidal.mat');
sample_offset = 1;
Texture_24 = file_24.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_25 = load('Multilayer_sensor_hard_4_sinusoidal.mat');
sample_offset = 1;
Texture_25 = file_25.finalData(2:19,sample_offset:sample_size+sample_offset-1);

file_26 = load('Multilayer_sensor_hard_6_sinusoidal.mat');
sample_offset = 1;
Texture_26 = file_26.finalData(2:19,sample_offset:sample_size+sample_offset-1);


%% Normalize voltage data for all textures

Texture_1_N = fs_dict(Texture_1);
Texture_2_N = fs_dict(Texture_2);
Texture_3_N = fs_dict(Texture_3);
Texture_4_N = fs_dict(Texture_4); 
Texture_5_N = fs_dict(Texture_5); 
Texture_6_N = fs_dict(Texture_6); 
Texture_7_N = fs_dict(Texture_7); 
Texture_8_N = fs_dict(Texture_8); 
Texture_9_N = fs_dict(Texture_9); 
Texture_10_N = fs_dict(Texture_10); 
Texture_11_N = fs_dict(Texture_11); 
Texture_12_N = fs_dict(Texture_12); 
Texture_13_N = fs_dict(Texture_13); 
Texture_14_N = fs_dict(Texture_14);
Texture_15_N = fs_dict(Texture_15);
Texture_16_N = fs_dict(Texture_16);
Texture_17_N = fs_dict(Texture_17); 
Texture_18_N = fs_dict(Texture_18); 
Texture_19_N = fs_dict(Texture_19); 
Texture_20_N = fs_dict(Texture_20); 
Texture_21_N = fs_dict(Texture_21); 
Texture_22_N = fs_dict(Texture_22); 
Texture_23_N = fs_dict(Texture_23); 
Texture_24_N = fs_dict(Texture_24); 
Texture_25_N = fs_dict(Texture_25); 
Texture_26_N = fs_dict(Texture_26); 

%% Combine data into one matrix

All_textures = [Texture_1_N Texture_2_N Texture_3_N Texture_4_N Texture_5_N Texture_6_N Texture_7_N Texture_8_N Texture_9_N Texture_10_N Texture_11_N Texture_12_N Texture_13_N Texture_14_N Texture_15_N Texture_16_N Texture_17_N Texture_18_N Texture_19_N Texture_20_N Texture_21_N Texture_22_N Texture_23_N Texture_24_N Texture_25_N Texture_26_N];

save('All_textures', 'All_textures')

%% Neuromorphic encoding
% clc;
% clear all;

load('Hybrid_finger-All_textures');


% Neuron model constants

% Slowly adapting 1 - Tonic Spiking   
a_SA_1= 0.02; % decay rate -> time constant
b_SA_1= 0.2; % spike sensitivity
c_SA_1= -65; % resting potential
d_SA_1= 8; % reset value -> adaption

% Rapidly adapting 1 - Phasic Bursting
a_RA_1= 0.02; % decay rate -> time constant
b_RA_1= 0.25; % spike sensitivity
c_RA_1= -55; % resting potential
d_RA_1= 0.05; % reset value -> adaption

% Slowly adapting 2 - Fast Spiking
a_SA_2= 0.1; % decay rate -> time constant
b_SA_2= 0.2; % spike sensitivity
c_SA_2= -65; % resting potential
d_SA_2= 2; % reset value -> adaption

% Rapidly adapting 2 - Phasic Spiking
a_RA_2= 0.02; % decay rate -> time constant
b_RA_2= 0.25; % spike sensitivity
c_RA_2= -65; % resting potential
d_RA_2= 6; % reset value -> adaption

scalefactor_SA_1 = 75;
scalefactor_RA_1 = 7;
scalefactor_SA_2 = 75;
scalefactor_RA_2 = 0.005;

% Separate data for encoding

 % Outer layer encoded with SA_1 and RA_1
 % Middle layer encoded with SA_2
 % Inner layer encoded with RA_2

Outer_layer = All_textures(10:18,:); 
Inner_layer = All_textures([2,3,5,6,8,9],:); 
Piezoelectric = All_textures(7,:); % Streamed data is repeated as 1,4,&7. Only one Pieoelectric channel used 


% Transform Data from Time Domain to Spike Information Domain
 
trial_window = 8050;

% (1) Outer layer SA1
    for i=trial_window:trial_window:length(All_textures)    % for each trial for each texture - 40 total trials 
        vec_outer = Outer_layer(:,i-trial_window+1:i);
        [v_SA_1,u_SA_1]=genspikes(vec_outer,scalefactor_SA_1,a_SA_1,b_SA_1,c_SA_1,d_SA_1); % generate neuron output i.e. spikes
        sr_SA_1=compspikerate(v_SA_1,100); % bin length of 100 ms
        for j=1:9
            SRm_SA_1(j)=mean(sr_SA_1(j,sr_SA_1(j,:)>0),2)*1000; % average nonzero spike rate
            if isnan(SRm_SA_1(j))
                SRm_SA_1(j)=0;
            end
        end
        [ISIm_SA_1,ISI]=compISI(v_SA_1,1000);  % ISIs within 1000 ms 
        Spiking_data_SA_1(i/trial_window,:)=[SRm_SA_1 ISIm_SA_1];
    end


% (2) Outer layer RA1
    for i=trial_window:trial_window:length(All_textures)    % for each trial for each texture - 40 total trials right now
        vec_outer = Outer_layer(:,i-trial_window+1:i);
        [v_RA_1,u_RA_1]=genspikes(vec_outer,scalefactor_RA_1,a_RA_1,b_RA_1,c_RA_1,d_RA_1); % generate neuron output i.e. spikes
        sr_RA_1=compspikerate(v_RA_1,100); % bin length of 100 ms
        for j=1:9
            SRm_RA_1(j)=mean(sr_RA_1(j,sr_RA_1(j,:)>0),2)*1000; % average nonzero spike rate
            if isnan(SRm_RA_1(j))
                SRm_RA_1(j)=0;
            end
        end
        [ISIm_RA_1,ISI]=compISI(v_RA_1,1000);  % ISIs within 1000 ms
        Spiking_data_RA_1(i/trial_window,:)=[SRm_RA_1 ISIm_RA_1];
    end

    
% (3) Middle layer SA2
    for i=trial_window:trial_window:length(All_textures)    % for each trial for each texture - 40 total trials right now
        vec_inner = Inner_layer(:,i-trial_window+1:i);
        [v_SA_2,u_SA_2]=genspikes(vec_inner,scalefactor_SA_2,a_SA_2,b_SA_2,c_SA_2,d_SA_2); % generate neuron output i.e. spikes
        sr_SA_2=compspikerate(v_SA_2,100); % bin length of 100 ms
        for j=1:6
            SRm_SA_2(j)=mean(sr_SA_2(j,sr_SA_2(j,:)>0),2)*1000; % average nonzero spike rate
            if isnan(SRm_SA_2(j))
                SRm_SA_2(j)=0;
            end
        end
        [ISIm_SA_2,ISI]=compISI(v_SA_2,1000);  % ISIs within 1000 ms 
        Spiking_data_SA_2(i/trial_window,:)=[SRm_SA_2 ISIm_SA_2];
    end


% (4) Inner layer RA2
    for i=trial_window:trial_window:length(All_textures)    % for each trial for each texture - 40 total trials right now
        vec_peizoelectric = Piezoelectric(:,i-trial_window+1:i);
        [v_RA_2,u_RA_2]=genspikes(vec_peizoelectric,scalefactor_RA_2,a_RA_2,b_RA_2,c_RA_2,d_RA_2); % generate neuron output i.e. spikes
        sr_RA_2=compspikerate(v_RA_2,100); % bin length of 100 ms
        for j=1
            SRm_RA_2(j)=mean(sr_RA_2(j,sr_RA_2(j,:)>0),2)*1000; % average nonzero spike rate
            if isnan(SRm_RA_2(j))
                SRm_RA_2(j)=0;
            end
        end
        [ISIm_RA_2,ISI]=compISI(v_RA_2,1000);  % ISIs within 1000 ms
        Spiking_data_RA_2(i/trial_window,:)=[SRm_RA_2 ISIm_RA_2];
    end

% Choosing which layers will be classified
    
Spiking_data = [Spiking_data_SA_1, Spiking_data_RA_1, Spiking_data_SA_2, Spiking_data_RA_2]; % all layers
% Spiking_data = [Spiking_data_SA_1, Spiking_data_RA_1]; % Only Outer layer
% Spiking_data = Spiking_data_SA_1; % Only Outer layer SA1
% Spiking_data = Spiking_data_RA_1; % Only Outer layer RA1
% Spiking_data = Spiking_data_SA_2; % Only Middle layer
% Spiking_data = Spiking_data_RA_2; % Only Inner layer
% Spiking_data = [Spiking_data_SA_1, Spiking_data_RA_1, Spiking_data_SA_2]; % Outer and middle layers
%Spiking_data = [Spiking_data_SA_1, Spiking_data_RA_1, Spiking_data_RA_2]; % Outer and inner layers
% Spiking_data = [Spiking_data_SA_2, Spiking_data_RA_2]; % Middle and inner layers

%save('Spiking data', 'Spiking_data')


%% k-fold Cross Validation and Classification

Accuracy = [];
Classifier_performance = [];
SVM_output_all = [];
ConfMat_all = [];
Labels_all = [];

for x = 1:1000

        k=4;
        members=40;
        trainlabels=[ones(1,members) ones(1,members)*2 ones(1,members)*3 ones(1,members)*4 ones(1,members)*5 ones(1,members)*6 ones(1,members)*7 ones(1,members)*8 ones(1,members)*9 ones(1,members)*10 ones(1,members)*11 ones(1,members)*12 ones(1,members)*13 ones(1,members)*14 ones(1,members)*15 ones(1,members)*16 ones(1,members)*17 ones(1,members)*18 ones(1,members)*19 ones(1,members)*20 ones(1,members)*21 ones(1,members)*22 ones(1,members)*23 ones(1,members)*24 ones(1,members)*25 ones(1,members)*26]';
        
        cvFolds = crossvalind('Kfold', trainlabels, k);
        cp = classperf(trainlabels);
        
        for i = 1:k                                  %# for each fold
            testIdx = (cvFolds == i);                %# get indices of test instances
            trainIdx = ~testIdx;                     %# get indices training instances
        
            %# train an SVM model over training instances
            SVM_mdl=fitcecoc(Spiking_data(trainIdx,:),trainlabels(trainIdx));
                display('Classifying with SVM...')
                tic
                SVM_output = predict(SVM_mdl,Spiking_data(testIdx,:));
                toc
             %# evaluate and update performance object
            cp = classperf(cp, SVM_output, testIdx);
          %SVM_acc=sum(SVM_output==testlabels)/length(testlabels);
        


          %ConfMat = confusionchart(trainlabels(testIdx), SVM_output,'fontsize',20);
          %ConfMat_all = [ConfMat_all, ConfMat];
        end

        % Save results from each iteration
         SVM_output_all = [SVM_output_all; SVM_output];
         Labels_all = [Labels_all; trainlabels(testIdx)];
         Classifier_performance = [Classifier_performance, cp];
         Accuracy = [Accuracy, cp.CorrectRate];
         
   
         
end

save('Hybrid_finger_SVM_results', 'Labels_all','SVM_output_all', 'Accuracy');

%%


Avg_Accuracy = mean(Accuracy);
 ConfMat = confusionchart(Labels_all, SVM_output_all,'Normalization','row-normalized','fontsize',30);
