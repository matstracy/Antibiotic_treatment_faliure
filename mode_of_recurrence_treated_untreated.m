function [] = mode_of_recurrence_treated_untreated(UTI_cases,params)
%% sensitive chance of failure
clear ratio1 fail_split1 all_sensitive1 all_sensitive_purchased1 sensitive_purchased_fails1 sensitive_purchased_fails_SS sensitive_purchased_fails_SR
clear ratio_noprch1 fail_split_noprch1 fail_split_adjusted fail_split_noprch_adjusted1 ratio_error1 number_prch1 ratio_noprch_error1 number_noprch1
clear ratio_res1 fail_split_res1     
clear ratio_noprch1 fail_split_noprch_res1      

%days = [4 8 15 22 29 36 43 50 57]
av_window = 7;
days = 4:1:60;
for kk = 1:length(days)-1
clear treatfailure
treatfailure = UTI_cases.date_diff >= days(kk) & UTI_cases.date_diff < days(kk)+av_window;


for drug = 1:params.number_drugs
    

     
   total_prch1(drug) = nnz(UTI_cases.PCR_sameday(:,drug));
 
   %all_sensitive = UTI_cases.SMP_Res(:,drug) == 1 | UTI_cases.SMP_Res(:,drug) == 2 ;
   all_sensitive1 = (UTI_cases.SMP_Res(:,drug) == 1 | UTI_cases.SMP_Res(:,drug) == 2) & UTI_cases.hasdiag ;
   all_resistant1 = (UTI_cases.SMP_Res(:,drug) == 3)  & UTI_cases.hasdiag ;
      
   all_sensitive_nextres1 =  UTI_cases.next_res(:,drug) == 1 | UTI_cases.SMP_Res(:,drug) == 2 ;
   all_resistant_nextres1 =  UTI_cases.next_res(:,drug) == 3 ;
   all_nextres = ismember(UTI_cases.next_res(:,drug),[1 2 3]);
   
   all_sensitive_purchased1 =  nnz((all_sensitive1 | all_resistant1) & UTI_cases.PCR_sameday(:,drug));
   sensitive_purchased_fails1 = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,drug) & treatfailure);
   sensitive_purchased_nextresfails1 = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_nextres);
   
   sensitive_purchased_fails_SS = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_sensitive_nextres1);
   sensitive_purchased_fails_SR = nnz(all_sensitive1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_resistant_nextres1);
   
   ratio1(drug,kk) = sensitive_purchased_fails1/all_sensitive_purchased1*100;
   number_prch1(drug) = all_sensitive_purchased1;
   ratio_error1(drug ,kk) = sqrt((ratio1(drug,kk)/100).*(1-(ratio1(drug,kk)/100))./(all_sensitive_purchased1))*100;
   fail_split1(1:2,kk, drug) = [sensitive_purchased_fails_SS/sensitive_purchased_nextresfails1 ; sensitive_purchased_fails_SR/sensitive_purchased_nextresfails1]*ratio1(drug,kk);
   SR_error1(drug,kk) = sqrt((fail_split1(2,kk, drug)/100).*(1-(fail_split1(2,kk, drug)/100))./(all_sensitive_purchased1))*100;
   
   
   %all_resistant_purchased =  nnz(all_resistant & UTI_cases.PCR_sameday(:,drug));
   all_resistant_purchased = nnz((all_sensitive1 | all_resistant1) & UTI_cases.PCR_sameday(:,drug));
   resistant_purchased_fails = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,drug) & treatfailure);
   resistant_purchased_nextresfails = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_nextres);
   resistant_purchased_fails_RS = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_sensitive_nextres1);
   resistant_purchased_fails_RR = nnz(all_resistant1 & UTI_cases.PCR_sameday(:,drug) & treatfailure & all_resistant_nextres1);
   
   ratio_res1(drug,kk) = resistant_purchased_fails/all_resistant_purchased*100;
   number_prch_res1(drug) = all_resistant_purchased;
   ratio_error_res1(drug,kk) = sqrt((ratio_res1(drug,kk)/100).*(1-(ratio_res1(drug,kk)/100))./(all_resistant_purchased))*100;
   fail_split_res1(1:2,kk, drug) = [resistant_purchased_fails_RR/resistant_purchased_nextresfails ; resistant_purchased_fails_RS/resistant_purchased_nextresfails]*ratio_res1(drug,kk);
   SR_error_res1(drug,kk) = sqrt((fail_split_res1(2,kk, drug)/100).*(1-(fail_split_res1(2,kk, drug)/100))./(all_resistant_purchased))*100;
   
   
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

end

%%
set(0, 'DefaultFigureRenderer', 'painters');
% figure
% set(gcf,'name','Fig. S2');
% set(gcf,'units','centimeters','Position',[1 1 18 30]);
first_set = [1 2 5 6 9 10 13 14];
for drug = 1:params.number_drugs
subplot(8,2,first_set(drug)) 
tots = [fail_split1(:,:,drug)',fail_split_res1(:,:,drug)'];
tots(isnan(tots))= 0;
total_vals = sum(tots,2);
%total_vals = sum([fail_split1(:,:,drug)',fail_split_res1(:,:,drug)'],2);
b1 = bar(1:length(days)-1,[fail_split1(:,:,drug)',fail_split_res1(:,:,drug)']./total_vals, 'stacked' , 'FaceColor','flat', 'EdgeColor', 'none', 'BarWidth', 1);
b1(1).CData(:,:) = ones(length(days)-1,1)*params.SS_color;
b1(2).CData(:,:) = ones(length(days)-1,1)*params.SR_color;
b1(3).CData(:,:) = ones(length(days)-1,1)*params.RR_color;
b1(4).CData(:,:) = ones(length(days)-1,1)*params.RS_color;
ylabel('% of early recurrences');
set(gca,'XTick',1:8:length(days));
set(gca, 'XTickLabel', [])
%set(gca,'xticklabel',days_string);

xlim([0.5 length(days)-0.5])
ylim([0 1])
set(gca,'YTick',0:0.5:1);
nums_string = num2str(((0:0.5:1)*100)');
set(gca,'yticklabel',nums_string);
title(UTI_cases.SMP_Res_drug_names(drug));
end


%
for drug = 1:params.number_drugs
subplot(8,2,first_set(drug)+2) 
total_vals = sum([fail_split_noprch1(:,:,drug)',fail_split_noprch_res1(:,:,drug)'],2);
b1 = bar(1:length(days)-1,[fail_split_noprch1(:,:,drug)',fail_split_noprch_res1(:,:,drug)']./total_vals, 'stacked' , 'FaceColor','flat',  'EdgeColor', 'none', 'BarWidth', 1);
b1(1).CData(:,:) = ones(length(days)-1,1)*params.SS_color;
b1(2).CData(:,:) = ones(length(days)-1,1)*params.SR_color;
b1(3).CData(:,:) = ones(length(days)-1,1)*params.RR_color;
b1(4).CData(:,:) = ones(length(days)-1,1)*params.RS_color;
%ylabel('% of treatment failures');
set(gca,'XTick',1:8:length(days));
days_string = num2str((days(1:8:length(days)))');
set(gca,'xticklabel',days_string);
xlabel('number of days between first and second sample');
xlim([0.5 length(days)-0.5])
ylim([0 1])
set(gca,'YTick',0:0.5:1);
set(gca,'yticklabel',nums_string);
end
end