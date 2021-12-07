function [] = risk_recurrence_by_day(UTI_cases,params)
%% sensitive chance of failure
clear ratio1 fail_split1 all_sensitive1 all_sensitive_purchased1 sensitive_purchased_fails1 sensitive_purchased_fails_SS sensitive_purchased_fails_SR
clear ratio_noprch1 fail_split_noprch1 fail_split_adjusted fail_split_noprch_adjusted1 ratio_error1 number_prch1 ratio_noprch_error1 number_noprch1
clear ratio_res1 fail_split_res1     
clear ratio_noprch1 fail_split_noprch_res1 fail_split2 fail_split_res2  

av_window = 7;
days = 4:1:50;
for kk = 1:length(days)-1
clear treatfailure
treatfailure = UTI_cases.date_diff >= days(kk) & UTI_cases.date_diff < days(kk)+av_window;

for drug = 1:params.number_drugs
 
   total_prch1(drug) = nnz(UTI_cases.PCR_sameday(:,drug));
 
   %all_sensitive = UTI_cases.SMP_Res(:,drug) == 1 | UTI_cases.SMP_Res(:,drug) == 2 ;
   all_sensitive1 = ismember(UTI_cases.SMP_Res(:,drug), params.sensitive_group) & UTI_cases.hasdiag;
   all_resistant1 = ismember(UTI_cases.SMP_Res(:,drug), params.resistant_group) & UTI_cases.hasdiag;
      
   all_sensitive_nextres1 =  ismember(UTI_cases.next_res(:,drug), params.sensitive_group) ;
   all_resistant_nextres1 =  ismember(UTI_cases.next_res(:,drug), params.resistant_group) ;
   all_nextres = ismember(UTI_cases.next_res(:,drug),[1 2 3]);
   
   all_sensitive_purchased1(drug ,kk) =  nnz((all_sensitive1 | all_resistant1) & UTI_cases.PCR_sameday(:,drug));
   sensitive_purchased_fails1(drug ,kk) = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,drug) & treatfailure);
   sensitive_purchased_nextresfails1(drug ,kk) = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_nextres);
   
   sensitive_purchased_fails_SS(drug ,kk) = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_sensitive_nextres1);
   sensitive_purchased_fails_SR(drug ,kk) = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_resistant_nextres1);
   
   ratio1(drug,kk) = sensitive_purchased_fails1(drug ,kk)/all_sensitive_purchased1(drug ,kk)*100;
   number_prch1(drug) = all_sensitive_purchased1(drug ,kk);
   ratio_error1(drug ,kk) = sqrt((ratio1(drug,kk)/100).*(1-(ratio1(drug,kk)/100))./(all_sensitive_purchased1(drug,kk)))*100;
   fail_split1(1:2,kk, drug) = [sensitive_purchased_fails_SS(drug ,kk)/sensitive_purchased_nextresfails1(drug ,kk) ; sensitive_purchased_fails_SR(drug ,kk)/sensitive_purchased_nextresfails1(drug ,kk)]*ratio1(drug,kk);
   SR_error1(drug,kk) = sqrt((fail_split1(2,kk, drug)/100).*(1-(fail_split1(2,kk, drug)/100))./(all_sensitive_purchased1(drug,kk)))*100;
   
   
   %all_resistant_purchased =  nnz(all_resistant & UTI_cases.PCR_sameday(:,drug));
   all_resistant_purchased(drug ,kk) = nnz((all_sensitive1 | all_resistant1) & UTI_cases.PCR_sameday(:,drug));
   resistant_purchased_fails(drug ,kk) = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,drug) & treatfailure);
   resistant_purchased_nextresfails(drug ,kk) = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_nextres);
   resistant_purchased_fails_RS(drug ,kk) = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_sensitive_nextres1);
   resistant_purchased_fails_RR(drug ,kk) = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_resistant_nextres1);
   
   ratio_res1(drug,kk) = resistant_purchased_fails(drug ,kk)/all_resistant_purchased(drug ,kk)*100;
   number_prch_res1(drug) = all_resistant_purchased(drug ,kk);
   ratio_error_res1(drug,kk) = sqrt((ratio_res1(drug,kk)/100).*(1-(ratio_res1(drug,kk)/100))./(all_resistant_purchased(drug,kk)))*100;
   fail_split_res1(1:2,kk, drug) = [resistant_purchased_fails_RR(drug ,kk)/resistant_purchased_nextresfails(drug ,kk) ; resistant_purchased_fails_RS(drug ,kk)/resistant_purchased_nextresfails(drug ,kk)]*ratio_res1(drug,kk);
   SR_error_res1(drug,kk) = sqrt((fail_split_res1(2,kk, drug)/100).*(1-(fail_split_res1(2,kk, drug)/100))./(all_resistant_purchased(drug,kk)))*100;
   
   
   %just_cleared(drug) = (100 -ratio(drug)) +(sensitive_purchased_fails_SS/sensitive_purchased_nextresfails)*ratio(drug);
   %fail_split_adjusted(1,drug) = ((sensitive_purchased_fails_SR/sensitive_purchased_nextresfails)*ratio(drug))/just_cleared(drug)*100;
   
   %all_sensitive_noprch =  nnz(all_sensitive & UTI_cases.PCR_sameday(:,10));
   all_sensitive_noprch1 =  nnz((all_sensitive1 | all_resistant1) & UTI_cases.PCR_sameday(:,10));
   sensitive_noprch_fails1 = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,10) & treatfailure );
   sensitive_noprch_nextresfails1 = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,10) & treatfailure & all_nextres);
   
   sensitive_noprch_fails_SS = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,10) & treatfailure & all_sensitive_nextres1);
   sensitive_noprch_fails_SR = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,10) & treatfailure & all_resistant_nextres1);
   
   ratio_noprch1(drug,kk) = sensitive_noprch_fails1/all_sensitive_noprch1*100;
   number_noprch1(drug) = all_sensitive_noprch1;
   ratio_noprch_error1(drug,kk) = sqrt((ratio_noprch1(drug,kk)/100).*(1-(ratio_noprch1(drug,kk)/100))./(all_sensitive_noprch1))*100;
   fail_split_noprch1(1:2,kk, drug) = [sensitive_noprch_fails_SS/sensitive_noprch_nextresfails1 ; sensitive_noprch_fails_SR/sensitive_noprch_nextresfails1]*ratio_noprch1(drug,kk);
   SR_noprch_error1(drug,kk) = sqrt((fail_split_noprch1(2,kk, drug)/100).*(1-(fail_split_noprch1(2,kk, drug)/100))./(all_sensitive_noprch1))*100;
   
   
   %all_resistant_noprch =  nnz(all_resistant & UTI_cases.PCR_sameday(:,10));
   all_resistant_noprch1 =  nnz((all_sensitive1 | all_resistant1) & UTI_cases.PCR_sameday(:,10));
   resistant_noprch_fails1 = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,10) & treatfailure);
   resistant_noprch_nextresfails1 = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,10) & treatfailure & all_nextres);
   resistant_noprch_fails_RS = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,10) & treatfailure & all_sensitive_nextres1);
   resistant_noprch_fails_RR = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,10) & treatfailure & all_resistant_nextres1);
   
   ratio_noprch_res1(drug,kk) = resistant_noprch_fails1/all_resistant_noprch1*100;
   number_noprch_res1(drug) = all_resistant_noprch1;
   ratio_error_res1(drug,kk) = sqrt((ratio_noprch_res1(drug,kk)/100).*(1-(ratio_noprch_res1(drug,kk)/100))./(all_resistant_noprch1))*100;
   fail_split_noprch_res1(1:2,kk, drug) = [resistant_noprch_fails_RR/resistant_noprch_nextresfails1 ; resistant_noprch_fails_RS/resistant_noprch_nextresfails1]*ratio_noprch_res1(drug,kk);
   SR_error_noprch_res1(drug,kk) = sqrt((fail_split_noprch_res1(2,kk, drug)/100).*(1-(fail_split_noprch_res1(2,kk, drug)/100))./(all_resistant_noprch1))*100;
   %just_cleared_noprch(drug) = (100 -ratio_noprch(drug)) +(sensitive_noprch_fails_SS/sensitive_noprch_nextresfails)*ratio_noprch(drug);
   %fail_split_noprch_adjusted(1,drug) = ((sensitive_noprch_fails_SR/sensitive_noprch_nextresfails)*ratio_noprch(drug))/just_cleared_noprch(drug)*100;
   
   %total_error(drug) = sqrt(P.*(1-P)./(n))
   

end

fail_split1_all(1:2,kk, 1) = [sum(sensitive_purchased_fails_SS(1:params.number_drugs ,kk))/sum(sensitive_purchased_nextresfails1(1:params.number_drugs ,kk)) ; sum(sensitive_purchased_fails_SR(1:params.number_drugs ,kk))/sum(sensitive_purchased_nextresfails1(1:params.number_drugs ,kk))]*(sum(sensitive_purchased_fails1(1:params.number_drugs ,kk))/sum(all_sensitive_purchased1(1:params.number_drugs))*100);
fail_split_res1_all(1:2,kk, 1) = [sum(resistant_purchased_fails_RR(1:params.number_drugs ,kk))/sum(resistant_purchased_nextresfails(1:params.number_drugs ,kk)) ; sum(resistant_purchased_fails_RS(1:params.number_drugs ,kk))/sum(resistant_purchased_nextresfails(1:params.number_drugs ,kk))]*(sum(resistant_purchased_fails(1:params.number_drugs ,kk))/sum(all_resistant_purchased(1:params.number_drugs))*100);

      
end

fail_split2(:,:,1) = fail_split1_all(:,:,1);
fail_split_res2(:,:,1) = fail_split_res1_all(:,:,1);
%   fail_split2(:,:,1) =  mean(fail_split1(:,:,1:params.number_drugs),3);
%   fail_split_res2(:,:,1) =  mean(fail_split_res1(:,:,1:params.number_drugs),3);
%   

%b1 = bar(1:length(days)-1,[fail_split1(:,:,drug)',fail_split_res1(:,:,drug)']./av_window, 'stacked' , 'FaceColor','flat',  'EdgeColor', 'none', 'BarWidth', 1);
b1 = bar(1:length(days)-1,[fail_split2(:,:,1)',fail_split_res2(:,:,1)']./av_window, 'stacked' , 'FaceColor','flat',  'EdgeColor', 'none', 'BarWidth', 1);

b1(1).CData(:,:) = ones(length(days)-1,1)*params.SS_color;
b1(2).CData(:,:) = ones(length(days)-1,1)*params.SR_color;
b1(3).CData(:,:) = ones(length(days)-1,1)*params.RR_color;
b1(4).CData(:,:) = ones(length(days)-1,1)*params.RS_color;

ylabel('% of patients');
set(gca,'XTick',1:8:length(days));
days_string = num2str((days(1:8:length(days)))');
set(gca,'xticklabel',days_string);
xlabel('# of days between first and second sample');
xlim([0.5 length(days)-0.5])
hold on
ylim([0 0.6]);
yl = ylim;
plot([28-days(1)+1 28-days(1)+1],yl, 'k--');
end