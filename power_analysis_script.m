%power_analysis_for_infant_sensory_intervention_study
%this script performs all analyses for a power analysis to test whether
%enhanced sensory experience via physical intervention affects infant
%functional connectivty. Without a name for the grant, I call it the infant
%sensory intervention study. Please don't acronym it. 
%
%A few packages will be required to use the script:
%1. gramm, a data visualization package: https://github.com/piermorel/gramm
%Oscar's showM package is asked for but currently unused
%
%
%The script requires the presence of a data folder containing hard coded
%filenames -- this must be provided via request through Washington
%University in St. Louis, largely due to the nature of the data
%acquisition. A DUA will be requred to gain access to reproduce these
%analyses.

%% initial setup
%we start by declaring some variables including where the data are located
data_folder='/home/feczk001/shared/projects/FEZ_USERS/feczk001/chad_power_analysis//data/files_for_eric_predictive_power/';
%workbench command to use for running ciftiopen
wb_command='wb_command';
%paths to matlab packages we will need for loading ciftis
matlab_cifti='/home/faird/shared/code/internal/utilities/Matlab_CIFTI/';
matlab_gifti='/home/faird/shared/code/external/utilities/gifti/';
%path to gramm for data visualizations
gramm_path='/home/faird/shared/code/external/utilities/gramm';
%path to showm for matrix visualization
showm_path='/home/faird/shared/code/internal/utilities/plotting-tools/showM/';
%add all relevant paths
addpath(genpath(showm_path))
addpath(genpath(gramm_path))
addpath(genpath(matlab_cifti))
addpath(genpath(matlab_gifti))

%% loading data 

%data from control participant
control_pretreat_precision = ciftiopen(strcat(data_folder,'/sub-5000_ses-split1_ROI-precision.ptseries.nii'),wb_command);
control_posttreat_precision = ciftiopen(strcat(data_folder,'/sub-5000_ses-split2_ROI-precision.ptseries.nii'),wb_command);
control_pretreat_gordon = ciftiopen(strcat(data_folder,'/sub-5000_ses-split1_ROI-gordon.ptseries.nii'),wb_command);
control_posttreat_gordon = ciftiopen(strcat(data_folder,'/sub-5000_ses-split2_ROI-gordon.ptseries.nii'),wb_command);

control_pretreat_tmask = dlmread(strcat(data_folder,'/5000_split1.tmask'));
control_postttreat_tmask = dlmread(strcat(data_folder,'/5000_split2.tmask'));


control_pretreat_precision_ts = control_pretreat_precision.cdata;
control_posttreat_precision_ts = control_posttreat_precision.cdata;
control_pretreat_gordon_ts = control_pretreat_gordon.cdata;
control_posttreat_gordon_ts = control_posttreat_gordon.cdata;

%data from positive participant
positive_pretreat_precision = ciftiopen(strcat(data_folder,'/sub-5003_ses-pre-massage_ROI-precision.ptseries.nii'),wb_command);
positive_posttreat_precision = ciftiopen(strcat(data_folder,'/sub-5003_ses-post-massage_ROI-precision.ptseries.nii'),wb_command);
positive_pretreat_gordon = ciftiopen(strcat(data_folder,'/sub-5003_ses-pre-massage_ROI-gordon.ptseries.nii'),wb_command);
positive_posttreat_gordon = ciftiopen(strcat(data_folder,'/sub-5003_ses-post-massage_ROI-gordon.ptseries.nii'),wb_command);

positive_pretreat_tmask = dlmread(strcat(data_folder,'/5003_pre-massage.tmask'));
positive_posttreat_tmask = dlmread(strcat(data_folder,'/5003_post-massage.tmask'));

positive_pretreat_precision_ts = positive_pretreat_precision.cdata;
positive_posttreat_precision_ts = positive_posttreat_precision.cdata;
positive_pretreat_gordon_ts = positive_pretreat_gordon.cdata;
positive_posttreat_gordon_ts = positive_posttreat_gordon.cdata;

%% calculate full correlations, so we know "ground truth"
control_pretreat_precision_pconn_full = corr(control_pretreat_precision_ts(:,logical(control_pretreat_tmask))');
control_posttreat_precision_pconn_full = corr(control_posttreat_precision_ts(:,logical(control_postttreat_tmask))');
control_pretreat_gordon_pconn_full = corr(control_pretreat_gordon_ts(:,logical(control_pretreat_tmask))');
control_posttreat_gordon_pconn_full = corr(control_posttreat_gordon_ts(:,logical(control_postttreat_tmask))');

control_treatment_effect_gordon_pconn = control_posttreat_gordon_pconn_full - control_pretreat_gordon_pconn_full;
control_treatment_effect_precision_pconn = control_posttreat_precision_pconn_full - control_pretreat_precision_pconn_full;

positive_pretreat_precision_pconn_full = corr(positive_pretreat_precision_ts(:,logical(positive_pretreat_tmask))');
positive_posttreat_precision_pconn_full = corr(positive_posttreat_precision_ts(:,logical(positive_posttreat_tmask))');
positive_pretreat_gordon_pconn_full = corr(positive_pretreat_gordon_ts(:,logical(positive_pretreat_tmask))');
positive_posttreat_gordon_pconn_full = corr(positive_posttreat_gordon_ts(:,logical(positive_posttreat_tmask))');

positive_treatment_effect_gordon_pconn = positive_posttreat_gordon_pconn_full - positive_pretreat_gordon_pconn_full;
positive_treatment_effect_precision_pconn = positive_posttreat_precision_pconn_full - positive_pretreat_precision_pconn_full;

%% save out effect pconns so we can visualize them in workbench
pconn_positive_precision_template = ciftiopen(strcat(data_folder,'/sub-5003_ses-pre-massage_ROI-precision.pconn.nii'),wb_command);
pconn_gordon_template = ciftiopen(strcat(data_folder,'/sub-5003_ses-pre-massage_ROI-gordon.pconn.nii'),wb_command);
pconn_control_precision_template = ciftiopen(strcat(data_folder,'/sub-5000_ses-split1_ROI-precision.pconn.nii'),wb_command);

pconn_positive_precision_template.cdata = positive_treatment_effect_precision_pconn;
pconn_gordon_template.cdata = positive_treatment_effect_gordon_pconn;

ciftisave(pconn_positive_precision_template,strcat(data_folder,'/sub-5003_ses-treatment-effect_ROI-precision.pconn.nii'),wb_command);
ciftisave(pconn_gordon_template,strcat(data_folder,'/sub-5003_ses-treatment-effect_ROI-gordon.pconn.nii'),wb_command);

pconn_control_precision_template.cdata = control_treatment_effect_precision_pconn;
pconn_gordon_template.cdata = control_treatment_effect_gordon_pconn;

ciftisave(pconn_control_precision_template,strcat(data_folder,'/sub-5000_ses-treatment-effect_ROI-precision.pconn.nii'),wb_command);
ciftisave(pconn_gordon_template,strcat(data_folder,'/sub-5000_ses-treatment-effect_ROI-gordon.pconn.nii'),wb_command);


%% calculate effect sizes and conducted statistical tests via resampling

%declare needed variables for this stage
TR=1.2;
gordon_ROIs = [23,191];
positive_precision_ROIs = [22,329];
control_precision_ROIs = [43,324];
frames_per_min = 60/TR;
frame_lengths = frames_per_min*5:frames_per_min*5:frames_per_min*60;
nsubs_per_group = 40;
sample_size_sets = [10 20 30 40];
npermutations = 1000;
time = [0 1];
%extract indices to randomly sample for tmasks
positive_pretreat_tmask_idx = find(positive_pretreat_tmask == 1);
positive_posttreat_tmask_idx = find(positive_posttreat_tmask == 1);
control_pretreat_tmask_idx = find(control_pretreat_tmask == 1);
control_posttreat_tmask_idx = find(control_postttreat_tmask == 1);
%declare tables for storing simulated data
control_pretreat_gordon_sim_matrices = zeros(nsubs_per_group,length(frame_lengths));
positive_pretreat_gordon_sim_matrices = control_pretreat_gordon_sim_matrices;
control_pretreat_precision_sim_matrices = control_pretreat_gordon_sim_matrices;
positive_pretreat_precision_sim_matrices = control_pretreat_gordon_sim_matrices;
control_posttreat_gordon_sim_matrices = zeros(nsubs_per_group,length(frame_lengths));
positive_posttreat_gordon_sim_matrices = control_pretreat_gordon_sim_matrices;
control_posttreat_precision_sim_matrices = control_pretreat_gordon_sim_matrices;
positive_posttreat_precision_sim_matrices = control_pretreat_gordon_sim_matrices;

gordon_power_matrix = zeros(length(sample_size_sets),length(frame_lengths));
precision_power_matrix = gordon_power_matrix;
%simulate all data needed for a comparison
rng('shuffle');
for iter = 1:npermutations
    for curr_subs = 1:nsubs_per_group
        frame_length_idx = 1;
        for curr_frame_length = frame_lengths
            temp_frames = control_pretreat_tmask_idx(randperm(length(control_pretreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(control_pretreat_gordon_ts(gordon_ROIs,temp_frames)');
            control_pretreat_gordon_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
    
            temp_frames = control_posttreat_tmask_idx(randperm(length(control_posttreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(control_posttreat_gordon_ts(gordon_ROIs,temp_frames)');
            control_posttreat_gordon_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
    
            temp_frames = positive_pretreat_tmask_idx(randperm(length(positive_pretreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(positive_pretreat_gordon_ts(gordon_ROIs,temp_frames)');
            positive_pretreat_gordon_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
            
            temp_frames = positive_posttreat_tmask_idx(randperm(length(positive_posttreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(positive_posttreat_gordon_ts(gordon_ROIs,temp_frames)');
            positive_posttreat_gordon_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
      
            temp_frames = control_pretreat_tmask_idx(randperm(length(control_pretreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(control_pretreat_precision_ts(control_precision_ROIs,temp_frames)');
            control_pretreat_precision_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
    
            temp_frames = control_posttreat_tmask_idx(randperm(length(control_posttreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(control_posttreat_precision_ts(control_precision_ROIs,temp_frames)');
            control_posttreat_precision_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
    
            temp_frames = positive_pretreat_tmask_idx(randperm(length(positive_pretreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(positive_pretreat_precision_ts(positive_precision_ROIs,temp_frames)');
            positive_pretreat_precision_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
            
            temp_frames = positive_posttreat_tmask_idx(randperm(length(positive_posttreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(positive_posttreat_precision_ts(positive_precision_ROIs,temp_frames)');
            positive_posttreat_precision_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
           
            frame_length_idx = frame_length_idx + 1;
        end
    end
    %prep simulated data for comparisons
    %calculate comparisons and extract p values for power
    curr_subset_count = 1;
    for curr_subset = sample_size_sets
        for curr_frame = 1:length(frame_lengths)
            gordon_analysis_table = table([control_pretreat_gordon_sim_matrices(randperm(nsubs_per_group,curr_subset),curr_frame) ; positive_pretreat_gordon_sim_matrices(randperm(nsubs_per_group,curr_subset),curr_frame) ], [ control_posttreat_gordon_sim_matrices(randperm(nsubs_per_group,curr_subset),curr_frame) ; positive_posttreat_gordon_sim_matrices(randperm(nsubs_per_group,curr_subset),curr_frame)] , [repmat("control",curr_subset,1) ; repmat("positive",curr_subset,1)],'VariableNames',["t0","t1","Group"] );
            gordon_model = fitrm(gordon_analysis_table,'t0-t1 ~ Group','WithinDesign',time);
            gordon_anova_table = ranova(gordon_model);
            if table2array(gordon_anova_table(2,5)) <= 0.05
                gordon_power_matrix(curr_subset_count,curr_frame) = gordon_power_matrix(curr_subset_count,curr_frame) + 1;
            end
            precision_analysis_table = table([control_pretreat_precision_sim_matrices(randperm(nsubs_per_group,curr_subset),curr_frame) ; positive_pretreat_precision_sim_matrices(randperm(nsubs_per_group,curr_subset),curr_frame) ], [ control_posttreat_precision_sim_matrices(randperm(nsubs_per_group,curr_subset),curr_frame) ; positive_posttreat_precision_sim_matrices(randperm(nsubs_per_group,curr_subset),curr_frame)] , [repmat("control",curr_subset,1) ; repmat("positive",curr_subset,1)],'VariableNames',["t0","t1","Group"] );
            precision_model = fitrm(precision_analysis_table,'t0-t1 ~ Group','WithinDesign',time);
            precision_anova_table = ranova(precision_model);
            if table2array(precision_anova_table(2,5)) <= 0.05
                precision_power_matrix(curr_subset_count,curr_frame) = precision_power_matrix(curr_subset_count,curr_frame) + 1;
            end     
        end
        curr_subset_count = curr_subset_count + 1;
    end
end
gordon_power_matrix = gordon_power_matrix./npermutations;
precision_power_matrix = precision_power_matrix./npermutations;
gordon_power_80 = zeros(length(sample_size_sets),2);
precision_power_80 = gordon_power_80;
gordon_power_30 = gordon_power_80;
precision_power_30 = precision_power_80;
for iter = 1:length(sample_size_sets)
    temp_index = find(gordon_power_matrix(iter,:) >= 0.8);
    gordon_power_80(iter,1) = sample_size_sets(iter)*2;
    try
       gordon_power_80(iter,2) = frame_lengths(temp_index(1))*TR/60; 
    catch
        warning(strcat("sample size of ", num2str(sample_size_sets(iter)*2), " has low power for gordon matrices; setting to NA" ));
        gordon_power_80(iter,2) = NaN;
    end
    temp_index = find(precision_power_matrix(iter,:) >= 0.8);
    precision_power_80(iter,1) = sample_size_sets(iter)*2;
    try
        precision_power_80(iter,2) = frame_lengths(temp_index(1))*TR/60; 
    catch
        warning(strcat("sample size of ", num2str(sample_size_sets(iter)*2), " has low power for precision matrices; setting to NA" ));
        precision_power_80(iter,2) = NaN;
    end
end
figure()
power_80 = [gordon_power_80 ; precision_power_80];
power_80_labels = [repmat({'gordon'},length(sample_size_sets),1); repmat({'precision'},length(sample_size_sets),1)];
power_gramm = gramm("x",power_80(:,2),"y",power_80(:,1),"color",power_80_labels);
power_gramm.geom_line()
power_gramm.set_title("Data needed to achieve 80 percent power",'FontSize',16)
power_gramm.set_names("x","#minutes","y","#subjects","color","ROI")
power_gramm.draw()
figure()
for iter = 1:length(sample_size_sets)
    temp_index = find(gordon_power_matrix(iter,:) >= 0.3);
    gordon_power_30(iter,1) = sample_size_sets(iter)*2;
    try
       gordon_power_30(iter,2) = frame_lengths(temp_index(1))*TR/60; 
    catch
        warning(strcat("sample size of ", num2str(sample_size_sets(iter)*2), " has low power for gordon matrices; setting to NA" ));
        gordon_power_30(iter,2) = NaN;
    end
    temp_index = find(precision_power_matrix(iter,:) >= 0.3);
    precision_power_30(iter,1) = sample_size_sets(iter)*2;
    try
        precision_power_30(iter,2) = frame_lengths(temp_index(1))*TR/60; 
    catch
        warning(strcat("sample size of ", num2str(sample_size_sets(iter)*2), " has low power for precision matrices; setting to NA" ));
        precision_power_30(iter,2) = NaN;
    end
end
power_30 = [gordon_power_30 ; precision_power_30];
power_30_labels = [repmat({'gordon'},length(sample_size_sets),1); repmat({'precision'},length(sample_size_sets),1)];
power_gramm = gramm("x",power_30(:,2),"y",power_30(:,1),"color",power_30_labels);
power_gramm.geom_line()
power_gramm.set_title("Data needed to achieve 30 percent power",'FontSize',16)
power_gramm.set_names("x","#minutes","y","#subjects","color","ROI")
power_gramm.draw()