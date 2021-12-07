function [] = Summary_tables(UTI_cases, params)
%% Distribution of resistance to each antibiotic before and after treatment, separated by treatments with each of the antibiotics.

treatfailure= UTI_cases.treatfailure;
Names = UTI_cases.SMP_Res_drug_names;



for treat_drug = 1:params.number_drugs+1   
   index = find(UTI_cases.PCR_sameday(:,treat_drug)  & UTI_cases.treatfailure ); 
   for test_drug = 1:params.number_drugs 
   num_SR_before_after(treat_drug,test_drug*4-3) = nnz(ismember(UTI_cases.SMP_Res(index,test_drug),[1 2]));
   num_SR_before_after(treat_drug,test_drug*4-2) = nnz(ismember(UTI_cases.SMP_Res(index,test_drug),3));
   num_SR_before_after(treat_drug,test_drug*4-1) = nnz(ismember(UTI_cases.next_res(index,test_drug),[1 2]));
   num_SR_before_after(treat_drug,test_drug*4) = nnz(ismember(UTI_cases.next_res(index,test_drug),3));
   
   end
end



%% Number of days between samples for treatment failures for each antibiotic and each mode of failure.
av_window = 7;
days = 4:1:50;
SS_SR_RR_RS = zeros(length(days)-1,4, params.number_drugs);
for kk = 1:length(days)-1
treatfailure = UTI_cases.date_diff >= days(kk) & UTI_cases.date_diff < days(kk)+av_window;

for drug = 1:params.number_drugs
 
   SS_SR_RR_RS(kk, 1, drug) = nnz(UTI_cases.PCR_sameday(:,drug) & treatfailure & ismember(UTI_cases.SMP_Res(:,drug), params.sensitive_group) & ismember(UTI_cases.next_res(:,drug), params.sensitive_group));
   SS_SR_RR_RS(kk, 2, drug) = nnz(UTI_cases.PCR_sameday(:,drug) & treatfailure & ismember(UTI_cases.SMP_Res(:,drug), params.sensitive_group) & ismember(UTI_cases.next_res(:,drug), params.resistant_group));
   SS_SR_RR_RS(kk, 3, drug) = nnz(UTI_cases.PCR_sameday(:,drug) & treatfailure & ismember(UTI_cases.SMP_Res(:,drug), params.resistant_group) & ismember(UTI_cases.next_res(:,drug), params.resistant_group));
   SS_SR_RR_RS(kk, 4, drug) = nnz(UTI_cases.PCR_sameday(:,drug) & treatfailure & ismember(UTI_cases.SMP_Res(:,drug), params.resistant_group) & ismember(UTI_cases.next_res(:,drug), params.sensitive_group));
    
end
end

drug_days = table(days(1:end-1)');
drug_days.Properties.VariableNames = {'Day'};

for drug = 1:params.number_drugs
drug_days(:,drug*4-2) = table(SS_SR_RR_RS(:, 1, drug));
drug_days.Properties.VariableNames{drug*4-2} = [UTI_cases.SMP_Res_drug_names{drug} ' # S->S cases']; 
drug_days(:,drug*4-1) = table(SS_SR_RR_RS(:, 2, drug));
drug_days.Properties.VariableNames{drug*4-1} = [UTI_cases.SMP_Res_drug_names{drug} ' # S->R cases'];
drug_days(:,drug*4) = table(SS_SR_RR_RS(:, 3, drug));
drug_days.Properties.VariableNames{drug*4} = [UTI_cases.SMP_Res_drug_names{drug} ' # R->R cases'];
drug_days(:,drug*4+1) = table(SS_SR_RR_RS(:, 4, drug));
drug_days.Properties.VariableNames{drug*4+1} = [UTI_cases.SMP_Res_drug_names{drug} ' # R->S cases'];
%UTI_cases.SMP_Res_drug_names{drug}
end
filename = 'Tables/Summary_table2.xlsx';
writetable(drug_days,filename);

%% Species prevalence before and after treatment failures for each treatment antibiotic.

A = cell2mat(UTI_cases.bug_all);
[counts,bugs] = groupcounts(A);
species_prevalence_table = table(bugs,counts);
total = sum(counts);
for ii = 1:height(species_prevalence_table)
   names(ii,1) = UTI_cases.Bugs.Name( UTI_cases.Bugs.Code == species_prevalence_table.bugs(ii));
   percent(ii,1) = counts(ii)/total*100;
end
species_prevalence_table.Names = names;
species_prevalence_table.percent = round(percent);
species_prevalence_table = sortrows(species_prevalence_table,'counts','descend');
species_prevalence_table = species_prevalence_table(1:10,:);

which_bugs = species_prevalence_table.bugs;
species_prevalence_before_after = table([species_prevalence_table.Names; 'Other']);
species_prevalence_before_after.Properties.VariableNames{1} = 'Species';

for treat_drug = 1:params.number_drugs   
 index = find(UTI_cases.PCR_sameday(:,treat_drug)  & UTI_cases.treatfailure ); 
 all_before = cell2mat(UTI_cases.bug_all(index));
 all_after = cell2mat(UTI_cases.new_bug(index));
   for ii = 1:10
      prevalence_before(ii, treat_drug) = nnz(ismember(all_before, which_bugs(ii)));
      prevalence_after(ii, treat_drug) = nnz(ismember(all_after, which_bugs(ii)) );
   end
     prevalence_before(ii+1, treat_drug) = nnz(~ismember(all_before, which_bugs));
     prevalence_after(ii+1, treat_drug) = nnz(~ismember(all_after, which_bugs));
     species_prevalence_before_after(:,treat_drug*2) = table(prevalence_before(:, treat_drug));
     species_prevalence_before_after.Properties.VariableNames{treat_drug*2} = [UTI_cases.SMP_Res_drug_names{treat_drug} ' # inital infection']; 
     species_prevalence_before_after(:,treat_drug*2+1) = table(prevalence_after(:, treat_drug));
     species_prevalence_before_after.Properties.VariableNames{treat_drug*2+1} = [UTI_cases.SMP_Res_drug_names{treat_drug} ' # recurrent infection']; 

end

filename = 'Tables/Summary_table_3.xlsx';
writetable(species_prevalence_before_after,filename);

%% Patient demographics by treatment antibiotic and treatment outcome.
names ={'Female'; 'Male';'Pregnant';'Age 0-9';'Age 10-19';...
    'Age 20-29';'Age 30-39';'Age 40-49';'Age 50-59';'Age 60-69';'Age 70-79';'Age 80-89';'Age >90';...
    'Prev SMPs 0-1'; 'Prev SMPs 2-4';'Prev SMPs 5+'; 'diab'; 'CKD'; 'Dialysis'; 'Prev60daysCath'; 'anyPrevCath'};

 Demographics_table_cases = table(names);
for drug= 1:8
Demog_treated =  UTI_cases.Demog(UTI_cases.PCR_sameday(:,drug) == 1  & UTI_cases.SMP_Res(:, drug) ~= 0 ,:);
any_SRmeasurement_treated =  sum(UTI_cases.any_SRmeasurement_byyear(UTI_cases.PCR_sameday(:,drug) == 1  & UTI_cases.SMP_Res(:, drug) ~= 0,:,:),3);

predictive_numbers(1) = nnz(Demog_treated.Gender);
predictive_numbers(2) = nnz(~Demog_treated.Gender);
predictive_numbers(3) = nnz(Demog_treated.Preg);


for ii = 1:10
predictive_numbers(ii+3) = nnz(Demog_treated.Age(:,ii));    
end

prev_smp_range = [0 1 4 ];
any_SRmeasurement = max(any_SRmeasurement_treated,[],2);
predictive_numbers(14) = nnz(ismember(any_SRmeasurement	,prev_smp_range(1):prev_smp_range(2)));
predictive_numbers(15) = nnz(ismember(any_SRmeasurement	,prev_smp_range(2)+1:prev_smp_range(3)));
predictive_numbers(16) = nnz(any_SRmeasurement	> prev_smp_range(3));

predictive_numbers(17) = nnz(Demog_treated.has_daibetes);
predictive_numbers(18) = nnz(Demog_treated.has_CKD);
predictive_numbers(19) = nnz(Demog_treated.has_dialysis);
predictive_numbers(20) = nnz(Demog_treated.any_prev60days_cath);
predictive_numbers(21) = nnz(Demog_treated.any_prev_cath);


predictive_per(1) = nnz(Demog_treated.Gender)/length(Demog_treated.Gender);
predictive_per(2) = nnz(~Demog_treated.Gender)/length(Demog_treated.Gender);
predictive_per(3) = nnz(Demog_treated.Preg)/length(Demog_treated.Gender);

for ii = 1:10
predictive_per(ii+3) = nnz(Demog_treated.Age(:,ii))/length(Demog_treated.Age(:,1)); 
end

predictive_per(14) = nnz(ismember(any_SRmeasurement	,prev_smp_range(1):prev_smp_range(2)))/length(any_SRmeasurement);
predictive_per(15) = nnz(ismember(any_SRmeasurement	,prev_smp_range(2)+1:prev_smp_range(3)))/length(any_SRmeasurement);
predictive_per(16) = nnz(any_SRmeasurement	> prev_smp_range(3))/length(any_SRmeasurement);
predictive_per(17) = nnz(Demog_treated.has_daibetes)/length(Demog_treated.has_daibetes);
predictive_per(18) = nnz(Demog_treated.has_CKD)/length(Demog_treated.has_CKD);
predictive_per(19) = nnz(Demog_treated.has_dialysis)/length(Demog_treated.has_dialysis);
predictive_per(20) = nnz(Demog_treated.any_prev60days_cath)/length(Demog_treated.any_prev60days_cath);
predictive_per(21) = nnz(Demog_treated.any_prev_cath)/length(Demog_treated.any_prev_cath);

% treatment failures
Demog_treated =  UTI_cases.Demog(UTI_cases.PCR_sameday(:,drug) == 1  & UTI_cases.SMP_Res(:, drug) ~= 0 & UTI_cases.treatfailure,:);
any_SRmeasurement_treated =  sum(UTI_cases.any_SRmeasurement_byyear(UTI_cases.PCR_sameday(:,drug) == 1  & UTI_cases.SMP_Res(:, drug) ~= 0 & UTI_cases.treatfailure,:,:),3);


predictive_numbers_fails(1) = nnz(Demog_treated.Gender);
predictive_numbers_fails(2) = nnz(~Demog_treated.Gender);
predictive_numbers_fails(3) = nnz(Demog_treated.Preg);

for ii = 1:10
predictive_numbers_fails(ii+3) = nnz(Demog_treated.Age(:,ii));    
end

any_SRmeasurement = max(any_SRmeasurement_treated,[],2);
predictive_numbers_fails(14) = nnz(ismember(any_SRmeasurement	,prev_smp_range(1):prev_smp_range(2)));
predictive_numbers_fails(15) = nnz(ismember(any_SRmeasurement	,prev_smp_range(2)+1:prev_smp_range(3)));
predictive_numbers_fails(16) = nnz(any_SRmeasurement > prev_smp_range(3));
predictive_numbers_fails(17) = nnz(Demog_treated.has_daibetes);
predictive_numbers_fails(18) = nnz(Demog_treated.has_CKD);
predictive_numbers_fails(19) = nnz(Demog_treated.has_dialysis);
predictive_numbers_fails(20) = nnz(Demog_treated.any_prev60days_cath);
predictive_numbers_fails(21) = nnz(Demog_treated.any_prev_cath);

predictive_per_fails(1) = nnz(Demog_treated.Gender)/length(Demog_treated.Gender);
predictive_per_fails(2) = nnz(~Demog_treated.Gender)/length(Demog_treated.Gender);
predictive_per_fails(3) = nnz(Demog_treated.Preg)/length(Demog_treated.Gender);

for ii = 1:10
predictive_per_fails(ii+3) = nnz(Demog_treated.Age(:,ii))/length(Demog_treated.Age(:,1)); 
end

predictive_per_fails(14) = nnz(ismember(any_SRmeasurement	,prev_smp_range(1):prev_smp_range(2)))/length(any_SRmeasurement);
predictive_per_fails(15) = nnz(ismember(any_SRmeasurement	,prev_smp_range(2)+1:prev_smp_range(3)))/length(any_SRmeasurement);
predictive_per_fails(16) = nnz(any_SRmeasurement	> prev_smp_range(3))/length(any_SRmeasurement);
predictive_per_fails(17) = nnz(Demog_treated.has_daibetes)/length(Demog_treated.has_daibetes);
predictive_per_fails(18) = nnz(Demog_treated.has_CKD)/length(Demog_treated.has_CKD);
predictive_per_fails(19) = nnz(Demog_treated.has_dialysis)/length(Demog_treated.has_dialysis);
predictive_per_fails(20) = nnz(Demog_treated.any_prev60days_cath)/length(Demog_treated.any_prev60days_cath);
predictive_per_fails(21) = nnz(Demog_treated.any_prev_cath)/length(Demog_treated.any_prev_cath);


 Demographics_table_cases(:,drug*4-2) = table(predictive_numbers');
 Demographics_table_cases.Properties.VariableNames{drug*4-2} = [UTI_cases.SMP_Res_drug_names{drug} ' # cases']; 
 Demographics_table_cases(:,drug*4-1) = table(round(predictive_per*100)');
 Demographics_table_cases.Properties.VariableNames{drug*4-1} = [UTI_cases.SMP_Res_drug_names{drug} ' # %']; 
 Demographics_table_cases(:,drug*4) = table(predictive_numbers_fails');
 Demographics_table_cases.Properties.VariableNames{drug*4} = [UTI_cases.SMP_Res_drug_names{drug} ' # fail cases']; 
 Demographics_table_cases(:,drug*4+1) = table(round(predictive_per_fails*100)');
 Demographics_table_cases.Properties.VariableNames{drug*4+1} = [UTI_cases.SMP_Res_drug_names{drug} ' % fails']; 
    


end






filename = 'Tables/Summary_table_4.xlsx';
writetable(Demographics_table_cases,filename);
 

