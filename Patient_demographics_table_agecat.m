function [] = Patient_demographics_table_agecat(FTRSp)
%% only with single PRCH and diagnosis
Demog_treated =  FTRSp.Demog(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:);
treatfailure_treated =  FTRSp.treatfailure(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:);
any_SRmeasurement_treated =  sum(FTRSp.any_SRmeasurement_byyear(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:,:),3);

any_SRmeasurement_treated =  FTRSp.any_SRmeasurement_byyear(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:,1) + ...
    FTRSp.any_SRmeasurement_byyear(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:,2)...
    + FTRSp.any_SRmeasurement_byyear(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:,3)...
    + FTRSp.any_SRmeasurement_byyear(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:,4);

% any_SRmeasurement_treated =  FTRSp.any_SRmeasurement_1to365days(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:) + ...
%     FTRSp.any_SRmeasurement_356to730days(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:) ...
%     + FTRSp.any_SRmeasurement_730to1095days(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:)...
%     + FTRSp.any_SRmeasurement_1095to1460day(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:);

%
%% only with single PRCH and diagnosis and susceptibiltiy measured

for drug= 1:8
 relavant_suceptib_measured(:,drug) = FTRSp.PCR_sameday(:,drug) == 1  & FTRSp.SMP_Res(:, drug) ~= 0;
 relavant_suceptib_notmeasured(:,drug) = FTRSp.PCR_sameday(:,drug) == 1  & FTRSp.SMP_Res(:, drug) == 0;
 suceptib_matched(:,drug) =  FTRSp.PCR_sameday(:,drug) == 1  & ismember(FTRSp.SMP_Res(:, drug), [1 2]); 
 suceptib_mismatched(:,drug) =  FTRSp.PCR_sameday(:,drug) == 1  & ismember(FTRSp.SMP_Res(:, drug), 3); 
 treat_fail(:,drug) =  FTRSp.PCR_sameday(:,drug) == 1  & FTRSp.treatfailure; 
end

relavant_suceptib_measured_any = (sum(relavant_suceptib_measured,2));

% nnz(relavant_suceptib_measured)
% 
% nnz(sum(suceptib_matched,2))
% 
% 
% a = nnz(sum(suceptib_matched,2) & FTRSp.hasdiag & FTRSp.treatfailure)
% b = nnz(sum(suceptib_mismatched,2) & FTRSp.hasdiag & FTRSp.treatfailure)
% c = nnz(sum(suceptib_matched,2) & FTRSp.hasdiag & ~FTRSp.treatfailure)
% d = nnz(sum(suceptib_mismatched,2) & FTRSp.hasdiag & ~FTRSp.treatfailure)
% 
% a + b 
% c + d
% 
% a + c 

%% test train
test_date_start = datenum('2018-05-01') ;
%test_date_start = datenum('2018-04-01') ;
test_date_end = datenum('2019-07-01') ;
test_date_range = test_date_start:test_date_end;
test_data = ismember(FTRSp.SamplingDate, test_date_range);


% a = nnz(sum(suceptib_matched,2) & FTRSp.hasdiag  & test_data)
% b = nnz(sum(suceptib_mismatched,2) & FTRSp.hasdiag  & test_data)
% c = nnz(sum(suceptib_matched,2) & FTRSp.hasdiag &  ~test_data)
% d = nnz(sum(suceptib_mismatched,2) & FTRSp.hasdiag & ~test_data)
% 
% a + c
% c + d
% 
% a + c 

%%

for drug = 1:8
val1(drug) = nnz(relavant_suceptib_notmeasured(:,drug));
end

FTRSp.hasdiag = FTRSp.hasdiag & relavant_suceptib_measured_any;
Demog_treated =  FTRSp.Demog(FTRSp.hasdiag & FTRSp.PCR_sameday_any ,:);
treatfailure_treated =  FTRSp.treatfailure(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:);
any_SRmeasurement_treated =  sum(FTRSp.any_SRmeasurement_byyear(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:,:),3);

% any_SRmeasurement_treated =  FTRSp.any_SRmeasurement_1to365days(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:) + ...
%     FTRSp.any_SRmeasurement_356to730days(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:) ...
%     + FTRSp.any_SRmeasurement_730to1095days(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:)...
%     + FTRSp.any_SRmeasurement_1095to1460day(FTRSp.hasdiag & FTRSp.PCR_sameday_any,:);

%%

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
%% treatment failures
Demog_treated =  FTRSp.Demog(FTRSp.hasdiag & FTRSp.PCR_sameday_any & FTRSp.treatfailure,:);
treatfailure_treated =  FTRSp.treatfailure(FTRSp.hasdiag & FTRSp.PCR_sameday_any & FTRSp.treatfailure,:);
any_SRmeasurement_treated =  FTRSp.any_SRmeasurement(FTRSp.hasdiag & FTRSp.PCR_sameday_any & FTRSp.treatfailure,:);
%
any_SRmeasurement_treated =  sum(FTRSp.any_SRmeasurement_byyear(FTRSp.hasdiag & FTRSp.PCR_sameday_any & FTRSp.treatfailure,:,:),3);

% any_SRmeasurement_treated =  FTRSp.any_SRmeasurement_1to365days(FTRSp.hasdiag & FTRSp.PCR_sameday_any & FTRSp.treatfailure,:) + ...
%     FTRSp.any_SRmeasurement_356to730days(FTRSp.hasdiag & FTRSp.PCR_sameday_any & FTRSp.treatfailure,:) ...
%     + FTRSp.any_SRmeasurement_730to1095days(FTRSp.hasdiag & FTRSp.PCR_sameday_any & FTRSp.treatfailure,:)...
%     + FTRSp.any_SRmeasurement_1095to1460day(FTRSp.hasdiag & FTRSp.PCR_sameday_any & FTRSp.treatfailure,:);

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


names ={'Female'; 'Male';'Pregnant';'Age 0-9';'Age 10-19';...
    'Age 20-29';'Age 30-39';'Age 40-49';'Age 50-59';'Age 60-69';'Age 70-79';'Age 80-89';'Age >90';...
    'Prev SMPs 0-1'; 'Prev SMPs 2-4';'Prev SMPs 5+'; 'diab'; 'CKD'; 'Dialysis'; 'Prev60daysCath'; 'anyPrevCath'};

 Demographics_table_cases = table(names,predictive_numbers',round(predictive_per*100)'...
    ,predictive_numbers_fails',round(predictive_per_fails*100)');

filename = 'Tables/Patient_demographics_table.xlsx';
writetable(Demographics_table_cases,filename);

% sum(Demographics_table_cases.Var2(1:2))
% sum(Demographics_table_cases.Var2(5:14))
end