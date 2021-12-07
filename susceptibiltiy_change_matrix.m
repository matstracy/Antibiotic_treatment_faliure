function [] = susceptibiltiy_change_matrix(FTRSp,params)
% check correlation between drugs
%fisk of aquiring resistane to the drug X and treatment drug
%risk of aquiring resistane to the drug X and NOT treatment drug
number_drugs_corr = params.number_drugs;

dates_to_use_start([1:4 6 8]) = min(FTRSp.SamplingDate);
dates_to_use_start(5) = min(FTRSp.SamplingDate)+7*321;
dates_to_use_start(7) = min(FTRSp.SamplingDate)+7*293;
dates_to_use_end(1:7) = max(FTRSp.SamplingDate);
dates_to_use_end(8) = min(FTRSp.SamplingDate)+7*293;


%% sensitive chance of failure
clear ratio2 fail_split2 all_sensitive2 all_sensitive_purchased2 sensitive_purchased_fails2 sensitive_purchased_fails_SS2 sensitive_purchased_fails_SR2
clear ratio_noprch2 fail_split_noprch2 fail_split_adjusted2 fail_split_noprch_adjusted2 ratio_error2 number_prch ratio_noprch_error2 number_noprch2
num_gained_resistance_to_test = zeros(number_drugs_corr);
num_gained_resistance_to_treat = zeros(number_drugs_corr);
num_gained_resistance_to_test_and_treat = zeros(number_drugs_corr);
num_gained_resistance_to_test_not_treat = zeros(number_drugs_corr);
num_gained_resistance_to_treat_not_test = zeros(number_drugs_corr);
total_num_treated = zeros(number_drugs_corr);

num_lost_resistance_to_test = zeros(number_drugs_corr);
num_lost_resistance_to_treat = zeros(number_drugs_corr);
num_lost_resistance_to_test_and_treat = zeros(number_drugs_corr);
num_lost_resistance_to_test_not_treat = zeros(number_drugs_corr);
num_lost_resistance_to_treat_not_test = zeros(number_drugs_corr);

num_gained_resistance_to_treat_testremained_sens = zeros(number_drugs_corr);
num_gained_resistance_to_test_treatremained_sens = zeros(number_drugs_corr);

for drug = 1:number_drugs_corr 
    total_prch(drug) = nnz(FTRSp.PCR_sameday(:,drug));
    
     dates_index_drug =  find(FTRSp.SamplingDate >= dates_to_use_start(drug) & FTRSp.SamplingDate <= dates_to_use_end(drug));
     
         
    for drug_to_test = 1:number_drugs_corr 
   
     dates_index_test =  find(FTRSp.SamplingDate >= dates_to_use_start(drug_to_test) & FTRSp.SamplingDate <= dates_to_use_end(drug_to_test));
  
     dates_index =  intersect(dates_index_test , dates_index_drug); 
     
   all_sensitive2 = (FTRSp.SMP_Res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2) & FTRSp.hasdiag(dates_index) ;
   all_sensitive_test = (FTRSp.SMP_Res(dates_index,drug_to_test) == 1 | FTRSp.SMP_Res(dates_index,drug_to_test) == 2) & FTRSp.hasdiag(dates_index) ;
   all_resistant2 = FTRSp.SMP_Res(dates_index,drug) == 3  & FTRSp.hasdiag(dates_index) ;
   all_resistant_test = FTRSp.SMP_Res(dates_index,drug_to_test) == 3  & FTRSp.hasdiag(dates_index) ;
   
   all_currentres = ismember(FTRSp.SMP_Res(dates_index,drug),[1 2]) & ismember(FTRSp.SMP_Res(dates_index,drug_to_test),[1 2 3]);
   
     
   all_sensitive_nextres =  FTRSp.next_res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2 ;
   all_resistant_nextres =  FTRSp.next_res(dates_index,drug) == 3 ;
   
   all_sensitive_nexttestres =  FTRSp.next_res(dates_index,drug_to_test) == 1 | FTRSp.SMP_Res(dates_index,drug_to_test) == 2 ;
   all_resistant_nexttestres =  FTRSp.next_res(dates_index,drug_to_test) == 3 ;
   
   all_nextres = ismember(FTRSp.next_res(dates_index,drug),[1 2 3]) & ismember(FTRSp.next_res(dates_index,drug_to_test),[1 2 3]);
   
   gained_resistance_to_test = all_sensitive_test & all_resistant_nexttestres & FTRSp.PCR_sameday(dates_index,drug);
   remained_sensitive_to_test = all_sensitive_test & all_sensitive_nexttestres & FTRSp.PCR_sameday(dates_index,drug);
   
   gained_resistance_to_treat = all_sensitive2 & all_resistant_nextres & FTRSp.PCR_sameday(dates_index,drug);
   remained_sensitive_to_treat = all_sensitive2 & all_sensitive_nextres & FTRSp.PCR_sameday(dates_index,drug);
   
   lost_resistance_to_test = all_resistant_test & all_sensitive_nexttestres & FTRSp.PCR_sameday(dates_index,drug);
   lost_resistance_to_treat = all_resistant2 & all_sensitive_nextres & FTRSp.PCR_sameday(dates_index,drug);
   
   
   total_num_treated(drug,drug_to_test) = nnz(FTRSp.PCR_sameday(dates_index,drug) & all_nextres & all_currentres);
   num_gained_resistance_to_test(drug,drug_to_test) = nnz(gained_resistance_to_test & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_gained_resistance_to_treat(drug,drug_to_test) = nnz(gained_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_gained_resistance_to_test_and_treat(drug,drug_to_test) = nnz(gained_resistance_to_test & gained_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_gained_resistance_to_test_not_treat(drug,drug_to_test) = nnz(gained_resistance_to_test  & ~gained_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres) ;
   num_gained_resistance_to_treat_not_test(drug,drug_to_test) = nnz(~gained_resistance_to_test  & gained_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_gained_resistance_to_test_treatremained_sens(drug,drug_to_test) = nnz(gained_resistance_to_test  & remained_sensitive_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres) ;
   num_gained_resistance_to_treat_testremained_sens(drug,drug_to_test) = nnz(gained_resistance_to_treat  & remained_sensitive_to_test & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres) ;
   
   num_lost_resistance_to_test(drug,drug_to_test) = nnz(lost_resistance_to_test & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_treat(drug,drug_to_test) = nnz(lost_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_test_and_treat(drug,drug_to_test) = nnz(lost_resistance_to_test & lost_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_test_not_treat(drug,drug_to_test) = nnz(lost_resistance_to_test  & ~lost_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_treat_not_test(drug,drug_to_test) = nnz(~lost_resistance_to_test  & lost_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
  
   
end
end 

%rates
rate_gained_resistance_to_treat = num_gained_resistance_to_treat./total_num_treated*100;
rate_gained_resistance_to_test = num_gained_resistance_to_test./total_num_treated*100;
rate_gained_resistance_to_test_not_treat = num_gained_resistance_to_test_not_treat./total_num_treated*100;
rate_gained_resistance_to_treat_not_test = num_gained_resistance_to_treat_not_test./total_num_treated*100;

rate_lost_resistance_to_treat = num_lost_resistance_to_treat./total_num_treated*100;
rate_lost_resistance_to_test = num_lost_resistance_to_test./total_num_treated*100;
rate_lost_resistance_to_test_not_treat = num_lost_resistance_to_test_not_treat./total_num_treated*100;
rate_lost_resistance_to_treat_not_test = num_lost_resistance_to_treat_not_test./total_num_treated*100;

rate_gained_resistance_to_test_treatremained_sens = num_gained_resistance_to_test_treatremained_sens./total_num_treated*100;
rate_gained_resistance_to_treat_testremained_sens = num_gained_resistance_to_treat_testremained_sens./total_num_treated*100;


%%

rate_gained_resistance_to_treat(rate_gained_resistance_to_treat == Inf) = 0;
rate_gained_resistance_to_test(rate_gained_resistance_to_test == Inf) = 0;
rate_gained_resistance_to_test_not_treat(rate_gained_resistance_to_test_not_treat == Inf) = 0;
rate_gained_resistance_to_treat_not_test(rate_gained_resistance_to_treat_not_test == Inf) = 0;


%% net change in resistance
rate_lost_resistance_to_treat(rate_lost_resistance_to_treat == Inf) = 0;
rate_lost_resistance_to_test(rate_lost_resistance_to_test == Inf) = 0;
rate_lost_resistance_to_test_not_treat(rate_lost_resistance_to_test_not_treat == Inf) = 0;
rate_lost_resistance_to_treat_not_test(rate_lost_resistance_to_treat_not_test == Inf) = 0;

n = 10;
lost_res_num = 15;%-round(min(min(rot90(rate_gained_resistance_to_test-rate_lost_resistance_to_test)*n)));
gained_res_num = 55;%round(max(max(rot90(rate_gained_resistance_to_test-rate_lost_resistance_to_test)*n)));
total_num = gained_res_num + lost_res_num;

cyan = [ linspace(0, 1, gained_res_num)' ones(gained_res_num,1) ones(gained_res_num,1)];
magenta = [  ones(gained_res_num,1) linspace(0, 1, gained_res_num)' ones(gained_res_num,1)];
magenta = flipud(magenta);
cmp_grey = [ linspace(0, 1, gained_res_num)' linspace(0, 1, gained_res_num)' linspace(0, 1, gained_res_num)'];

cyan_white_magenta = [cyan; magenta(1:lost_res_num,:)];
colormap(flipud(cyan_white_magenta))


net_change_res = rate_gained_resistance_to_test-rate_lost_resistance_to_test;
net_change_res(isnan(rate_lost_resistance_to_test)) = 0;
net_change_res = net_change_res(params.new_order, params.new_order);


rate_gained_resistance_to_test_treated = rate_gained_resistance_to_test;

%% chance of failure untreated
clear ratio2 fail_split2 all_sensitive2 all_sensitive_purchased2 sensitive_purchased_fails2 sensitive_purchased_fails_SS2 sensitive_purchased_fails_SR2
clear ratio_noprch2 fail_split_noprch2 fail_split_adjusted2 fail_split_noprch_adjusted2 ratio_error2 number_prch ratio_noprch_error2 number_noprch2
num_gained_resistance_to_test = zeros(number_drugs_corr);
num_gained_resistance_to_treat = zeros(number_drugs_corr);
num_gained_resistance_to_test_and_treat = zeros(number_drugs_corr);
num_gained_resistance_to_test_not_treat = zeros(number_drugs_corr);
num_gained_resistance_to_treat_not_test = zeros(number_drugs_corr);
total_num_treated = zeros(number_drugs_corr);

num_lost_resistance_to_test = zeros(number_drugs_corr);
num_lost_resistance_to_treat = zeros(number_drugs_corr);
num_lost_resistance_to_test_and_treat = zeros(number_drugs_corr);
num_lost_resistance_to_test_not_treat = zeros(number_drugs_corr);
num_lost_resistance_to_treat_not_test = zeros(number_drugs_corr);

num_gained_resistance_to_treat_testremained_sens = zeros(number_drugs_corr);
num_gained_resistance_to_test_treatremained_sens = zeros(number_drugs_corr);

for drug = 1:number_drugs_corr 
    total_prch(drug) = nnz(FTRSp.PCR_sameday(:,10));
    
 dates_index_drug =  find(FTRSp.SamplingDate >= dates_to_use_start(drug) & FTRSp.SamplingDate <= dates_to_use_end(drug));
     
         
    for drug_to_test = 1:number_drugs_corr 
   
     dates_index_test =  find(FTRSp.SamplingDate >= dates_to_use_start(drug_to_test) & FTRSp.SamplingDate <= dates_to_use_end(drug_to_test));
  
     dates_index =  intersect(dates_index_test , dates_index_drug); 
    
   
 
   all_sensitive2 = (FTRSp.SMP_Res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2) & FTRSp.hasdiag(dates_index) ;
   all_sensitive_test = (FTRSp.SMP_Res(dates_index,drug_to_test) == 1 | FTRSp.SMP_Res(dates_index,drug_to_test) == 2) & FTRSp.hasdiag(dates_index) ;
   all_resistant2 = FTRSp.SMP_Res(dates_index,drug) == 3  & FTRSp.hasdiag(dates_index) ;
   all_resistant_test = FTRSp.SMP_Res(dates_index,drug_to_test) == 3  & FTRSp.hasdiag(dates_index) ;
   
   all_currentres = ismember(FTRSp.SMP_Res(dates_index,drug),[1 2 3]) & ismember(FTRSp.SMP_Res(dates_index,drug_to_test),[1 2 3]);
   
     
   all_sensitive_nextres =  FTRSp.next_res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2 ;
   all_resistant_nextres =  FTRSp.next_res(dates_index,drug) == 3 ;
   
   all_sensitive_nexttestres =  FTRSp.next_res(dates_index,drug_to_test) == 1 | FTRSp.SMP_Res(dates_index,drug_to_test) == 2 ;
   all_resistant_nexttestres =  FTRSp.next_res(dates_index,drug_to_test) == 3 ;
   
   all_nextres = ismember(FTRSp.next_res(dates_index,drug),[1 2 3]) & ismember(FTRSp.next_res(dates_index,drug_to_test),[1 2 3]);
   
   gained_resistance_to_test = all_sensitive_test & all_resistant_nexttestres & FTRSp.PCR_sameday(dates_index,10);
   remained_sensitive_to_test = all_sensitive_test & all_sensitive_nexttestres & FTRSp.PCR_sameday(dates_index,10);
   
   gained_resistance_to_treat = all_sensitive2 & all_resistant_nextres & FTRSp.PCR_sameday(dates_index,10);
   remained_sensitive_to_treat = all_sensitive2 & all_sensitive_nextres & FTRSp.PCR_sameday(dates_index,10);
   
   lost_resistance_to_test = all_resistant_test & all_sensitive_nexttestres & FTRSp.PCR_sameday(dates_index,10);
   lost_resistance_to_treat = all_resistant2 & all_sensitive_nextres & FTRSp.PCR_sameday(dates_index,10);
   
   
   total_num_treated(drug,drug_to_test) = nnz(FTRSp.PCR_sameday(dates_index,10) & all_nextres & all_currentres);
   num_gained_resistance_to_test(drug,drug_to_test) = nnz(gained_resistance_to_test & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_gained_resistance_to_treat(drug,drug_to_test) = nnz(gained_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_gained_resistance_to_test_and_treat(drug,drug_to_test) = nnz(gained_resistance_to_test & gained_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_gained_resistance_to_test_not_treat(drug,drug_to_test) = nnz(gained_resistance_to_test  & ~gained_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres) ;
   num_gained_resistance_to_treat_not_test(drug,drug_to_test) = nnz(~gained_resistance_to_test  & gained_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_gained_resistance_to_test_treatremained_sens(drug,drug_to_test) = nnz(gained_resistance_to_test  & remained_sensitive_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres) ;
   num_gained_resistance_to_treat_testremained_sens(drug,drug_to_test) = nnz(gained_resistance_to_treat  & remained_sensitive_to_test & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres) ;
   
   num_lost_resistance_to_test(drug,drug_to_test) = nnz(lost_resistance_to_test & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_treat(drug,drug_to_test) = nnz(lost_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_test_and_treat(drug,drug_to_test) = nnz(lost_resistance_to_test & lost_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_test_not_treat(drug,drug_to_test) = nnz(lost_resistance_to_test  & ~lost_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_treat_not_test(drug,drug_to_test) = nnz(~lost_resistance_to_test  & lost_resistance_to_treat & FTRSp.treatfailure(dates_index) & all_nextres & all_currentres);
  
   
end
end 

%rates
rate_gained_resistance_to_treat = num_gained_resistance_to_treat./total_num_treated*100;
rate_gained_resistance_to_test = num_gained_resistance_to_test./total_num_treated*100;
rate_gained_resistance_to_test_not_treat = num_gained_resistance_to_test_not_treat./total_num_treated*100;
rate_gained_resistance_to_treat_not_test = num_gained_resistance_to_treat_not_test./total_num_treated*100;

rate_lost_resistance_to_treat = num_lost_resistance_to_treat./total_num_treated*100;
rate_lost_resistance_to_test = num_lost_resistance_to_test./total_num_treated*100;
rate_lost_resistance_to_test_not_treat = num_lost_resistance_to_test_not_treat./total_num_treated*100;
rate_lost_resistance_to_treat_not_test = num_lost_resistance_to_treat_not_test./total_num_treated*100;

rate_gained_resistance_to_test_treatremained_sens = num_gained_resistance_to_test_treatremained_sens./total_num_treated*100;
rate_gained_resistance_to_treat_testremained_sens = num_gained_resistance_to_treat_testremained_sens./total_num_treated*100;

%%
rate_gained_resistance_to_treat(rate_gained_resistance_to_treat == Inf) = 0;
rate_gained_resistance_to_test(rate_gained_resistance_to_test == Inf) = 0;
rate_gained_resistance_to_test_not_treat(rate_gained_resistance_to_test_not_treat == Inf) = 0;
rate_gained_resistance_to_treat_not_test(rate_gained_resistance_to_treat_not_test == Inf) = 0;

rate_gained_resistance_to_test_untreated = rate_gained_resistance_to_test;
D = diag( rate_gained_resistance_to_test_untreated );
rate_lost_resistance_to_test_untreated = rate_lost_resistance_to_test;
D2 = diag( rate_lost_resistance_to_test_untreated );
D3 = D-D2;
%%


image([rot90(net_change_res) flipud(D3)]*n+lost_res_num )
xlabel('treatment antibiotic')
ylabel('susceptibiltiy of reccurent UTI')
colormap(flipud(cyan_white_magenta))
xticks([1:1:number_drugs_corr+1])
set(gca,'xticklabel',[FTRSp.SMP_Res_drug_names(params.new_order); 'Untreated'])
xtickangle(45);
yticks([1:1:number_drugs_corr])
set(gca,'yticklabel',flipud(FTRSp.SMP_Res_drug_names(params.new_order)))
cbh = colorbar ; %Create Colorbar
x=5;
cbh.Ticks = linspace(0, length(colormap), 8)+x ; 
ticknums = (linspace(-lost_res_num, gained_res_num, 8))/n+x/n;
cbh.TickLabels = num2cell(ticknums) ; 
set(get(cbh,'label'),'string',{'net % of treated cases which changed';...
    'resistance to focal drug'});
axis image

end
