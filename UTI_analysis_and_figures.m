%% perform UTI analysis and make figures
clearvars -except UTI_cases

if ~exist('UTI_cases') % load the data
    load('UTI_cases.mat');
end
if ~isfolder('Tables')
      mkdir('Tables')
end

%% sensitivity grouping
%Sensitive = 1; Intermediate = 2; Resistant = 3;
params.number_drugs = 8; 
params.sensitive_group = [1 2];
SIR= 1:3; params.resistant_group = SIR(~ismember(SIR, params.sensitive_group));
clear SIR
% figure colors
params.SS_color = [ 25, 32, 128]/255; params.SR_color = [ 0 0.95 0.95];
params.RR_color = [ 128 0 128]/255; params.RS_color = [ 0.95 0 0.95];
%% generate Table S1 patient demographics 
Patient_demographics_table_agecat(UTI_cases);

%% Fig 1B
mismatched_vs_mathced_recurrence_rate(UTI_cases,params)
 
%% Fig 1C
mode_of_reccurence_pie_charts_hasdiag(UTI_cases,params)

%% Fig 1D 
params.new_order = [1 2 8 3:7]; %drug order for figures
figure
set(gcf,'color','w','name','Fig. 1D','units','centimeters','Position',[1 1 20 7])
risk_reccurence_correctly_prescribed_treat_v_untreat(UTI_cases,params)

%% Fig S3
figure
set(gcf,'color','w', 'name','Fig. S3','units','centimeters','Position',[1 1 17 7.5]);
net_change_resistance(UTI_cases,params)
%% Fig 1E rate of recuurence by day 
figure
set(gcf,'color','w','name','Fig. 1E','units','centimeters','Position',[1 1 11 6]);
risk_recurrence_by_day(UTI_cases,params)
drawnow
%% Fig S6 Mode of recurrence by day for treated and untreated UTIs. Mode of recurrence by day for treated and untreated UTIs. 
figure
set(gcf,'color','w', 'name','Fig. S6 Mode of recurrence','units','centimeters','Position',[1 1 16 26]);
mode_of_recurrence_treated_untreated(UTI_cases,params)
drawnow
%% Fig 1F matrix UTIs
figure
set(gcf,'color','w', 'name','Fig. 2F Susceptibiltiy change matrix', 'units','centimeters','Position',[1 1 12 10]);
susceptibiltiy_change_matrix(UTI_cases,params)
drawnow
%% Fig 2E
figure;
set(gcf,'color','w', 'name','Fig. 2F', 'units','centimeters','Position',[1 1 18 10]);
initally_ecoli_gained_res_changed_bac_errors_remainedS(UTI_cases,params)

%% Fig 2G
figure
set(gcf,'color','w', 'name','Fig. 2G UTIs', 'units','centimeters','Position',[1 1 5 9]);
changed_res_bac_alldrugs(UTI_cases,params)
%% Fig 3B
figure;
set(gcf,'color','w', 'name','Fig. 3B Adjusted OR');
correct_treatment_adjusted_OR(UTI_cases,params)
drawnow
%% Fig 3DEF
correct_treatment_fails_reccomend_drug_v11(UTI_cases,params)
drawnow

%% Supp Fig S5 
Adjusted_risk_of_recurrence_treat_notreat_agecont(UTI_cases,params)
%% Supp Fig S8
figure;
set(gcf,'color','w', 'name','Fig. S8 Rates of resistance by species UTIs',...
'units','centimeters','Position',[1 1 25 18]);
rate_res_by_species(UTI_cases);

%% Supp Fig S10
Prev_res_OR_byYear2(UTI_cases,params)

%% Supp Fig S11
figure;
set(gcf,'color','w', 'name','Fig. S11', 'units','centimeters','Position',[1 1 25 18]);
correct_treatment_make_figs_regression_diag_bar_PURCHandRES(UTI_cases,params)
Prev_res_OR_exactly1(UTI_cases,params)
%% Supp Fig S14BC
correct_treatment_fails_reccomend_drug_v9_2testperiods(UTI_cases,params)
drawnow
%% Supp Fig 9
species_change_matrix_gained_res(UTI_cases,params)

%% Supp Fig S16 intermediate grouped with resistant
figure;
set(gcf,'color','w', 'name','Fig. S16', 'units','centimeters','Position',[1 1 25 10]);
subplot(1,2,1)
correct_treatment_adjusted_OR_S_IR(UTI_cases,params);
correct_treatment_fails_reccomend_drug_S_IR(UTI_cases,params)
drawnow

%%
species_prevalence_table(UTI_cases)
%%
if ~isfolder('Figures')
      mkdir('Figures')
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'name');
  saveas(FigHandle,['Figures/' FigName '.fig'])
end

