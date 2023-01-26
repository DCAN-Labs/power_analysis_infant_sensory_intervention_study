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

positive_pretreat_gordon_ALFF = dlmread(strcat(data_folder,'/5003_pre-massage_gordon_parcels.txt'));
positive_pretreat_precision_ALFF = dlmread(strcat(data_folder,'/5003_pre-massage_individ_parcels.txt'));

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

%% declare SMN labels for HB3 power analysis

SMN_precision_indices = sort([169 74 3 75 82 72 81 84 ...
    85 56 55 83 80 70 131 67 65 54 63 53 79 51 142 39 52 66 41 40 42 50 49 48 47 ...
    477 457 456 465 293 345 363 364 365 466 464 462 367 366 373 372 369 371 360 ...
    370 342 344 350 351 354 343 359 439 355 357 331 356 353 341 358 340 440 327 ...
    292 335 326 328 329 330 332 314 313 339 337]);

SMN_gordon_indices = sort([ 53 3 39 59 54 56 58 38 45 47 50 57 31 46 32 48 2 ...
    33 36 37 35 30 212 164 218 197 270 213 217 216 201 205 207 210 215 ...
    209 214 190 163 206 204 202 193 191 194 195]);

positive_gordon_SMN_meanFC = mean(positive_pretreat_gordon_pconn_full(SMN_gordon_indices,SMN_gordon_indices));
positive_precision_SMN_meanFC = mean(positive_pretreat_precision_pconn_full(SMN_precision_indices,SMN_gordon_indices));
positive_gordon_SMN_meaneffect = mean(positive_treatment_effect_gordon_pconn(SMN_gordon_indices,SMN_gordon_indices));
positive_precision_SMN_meaneffect = mean(positive_treatment_effect_precision_pconn(SMN_precision_indices,SMN_precision_indices));

positive_gordon_SMN_ALFF = positive_pretreat_gordon_ALFF(SMN_gordon_indices);
positive_precision_SMN_ALFF = positive_pretreat_precision_ALFF(SMN_precision_indices);

%% calculate effect sizes and conducted statistical tests via resampling

%declare needed variables for this stage
TR=1.2;
gordon_ROIs = [21,180];
positive_precision_ROIs = [44,338];
control_precision_ROIs = [35,328];
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
gordon_ALFF_power = zeros(length(frame_lengths),1);
precision_ALFF_power = gordon_ALFF_power;
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
            
            temp_pretreat_gordon_frames = positive_pretreat_tmask_idx(randperm(length(positive_pretreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(positive_pretreat_gordon_ts(gordon_ROIs,temp_pretreat_gordon_frames)');
            positive_pretreat_gordon_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
            
            temp_posttreat_gordon_frames = positive_posttreat_tmask_idx(randperm(length(positive_posttreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(positive_posttreat_gordon_ts(gordon_ROIs,temp_posttreat_gordon_frames)');
            positive_posttreat_gordon_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);    
            
            temp_frames = control_pretreat_tmask_idx(randperm(length(control_pretreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(control_pretreat_precision_ts(control_precision_ROIs,temp_frames)');
            control_pretreat_precision_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
    
            temp_frames = control_posttreat_tmask_idx(randperm(length(control_posttreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(control_posttreat_precision_ts(control_precision_ROIs,temp_frames)');
            control_posttreat_precision_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);

            temp_pretreat_precision_frames = positive_pretreat_tmask_idx(randperm(length(positive_pretreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(positive_pretreat_precision_ts(positive_precision_ROIs,temp_pretreat_precision_frames)');
            positive_pretreat_precision_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
            
            temp_posttreat_precision_frames = positive_posttreat_tmask_idx(randperm(length(positive_posttreat_tmask_idx),curr_frame_length),1);
            temp_corr = corr(positive_posttreat_precision_ts(positive_precision_ROIs,temp_posttreat_precision_frames)');
            positive_posttreat_precision_sim_matrices(curr_subs,frame_length_idx) = temp_corr(1,2);
            
            if curr_subs == 1
                gordon_pretreat_temp_mat = corr(positive_pretreat_gordon_ts(:,temp_pretreat_gordon_frames)');
                gordon_posttreat_temp_mat = corr(positive_posttreat_gordon_ts(:,temp_posttreat_gordon_frames)');
                precision_pretreat_temp_mat = corr(positive_pretreat_precision_ts(:,temp_pretreat_precision_frames)');
                precision_posttreat_temp_mat = corr(positive_posttreat_precision_ts(:,temp_posttreat_precision_frames)');
                gordon_effect_matrix = gordon_pretreat_temp_mat - gordon_posttreat_temp_mat;
                precision_effect_matrix = precision_pretreat_temp_mat - precision_posttreat_temp_mat;
                positive_gordon_SMN_meaneffect_perm = mean(gordon_effect_matrix(SMN_gordon_indices,SMN_gordon_indices));
                positive_precision_SMN_meaneffect_perm = mean(precision_effect_matrix(SMN_precision_indices,SMN_precision_indices));
                ALFF_FC_gordon = anova(fitlm(positive_gordon_SMN_meaneffect_perm,positive_gordon_SMN_ALFF));
                ALFF_FC_precision = anova(fitlm(positive_precision_SMN_meaneffect_perm,positive_precision_SMN_ALFF));
                if table2array(ALFF_FC_gordon(1,5)) <= 0.05
                    gordon_ALFF_power(frame_length_idx) = gordon_ALFF_power(frame_length_idx) + 1;
                end
                if table2array(ALFF_FC_precision(1,5)) <= 0.05
                    precision_ALFF_power(frame_length_idx) = precision_ALFF_power(frame_length_idx) + 1;
                end
            end
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
gordon_ALFF_power = gordon_ALFF_power./npermutations;
precision_ALFF_power = precision_ALFF_power./npermutations;

%% create data visualizations for power analysis
gordon_power_80 = zeros(length(sample_size_sets),2);
precision_power_80 = gordon_power_80;
gordon_power_30 = gordon_power_80;
precision_power_30 = precision_power_80;
power_plot_nsubs80 = zeros(length(frame_lengths)*2,2);
power_plot_nsubs80(1:length(frame_lengths),1) = gordon_power_matrix(end,:);
power_plot_nsubs80(length(frame_lengths)+1:end,1) = precision_power_matrix(end,:);
power_plot_nsubs80(1:length(frame_lengths),2) = frame_lengths.*(TR/60);
power_plot_nsubs80(length(frame_lengths)+1:end,2) = frame_lengths.*(TR/60);
figure()
power_nsubs80_labels = [repmat({'gordon'},length(frame_lengths),1); repmat({'precision'},length(frame_lengths),1)];
power_gramm = gramm("x",power_plot_nsubs80(:,2),"y",power_plot_nsubs80(:,1),"color",power_nsubs80_labels);
power_gramm.geom_line()
power_gramm.set_title("Power vs. # minutes: N = 80",'FontSize',16)
power_gramm.set_names("x","#minutes","y","%power","color","ROI")
power_gramm.draw()
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

figure()
precision_power_by_time = [precision_ALFF_power, frame_lengths'.*(1.2/60)];
precision_power_gramm = gramm("x",precision_power_by_time(:,2),"y",precision_power_by_time(:,1));
precision_power_gramm.geom_line()
precision_power_gramm.set_title("Power for ALFF-deltaFC association by time")
precision_power_gramm.set_names("x","#minutes","y","Power (%)")
precision_power_gramm.draw()

