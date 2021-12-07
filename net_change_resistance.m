function []= net_change_resistance(FTRSp,params)
%% sensitive chance of failure
clear ratio2 fail_split2 all_sensitive2 all_sensitive_purchased2 sensitive_purchased_fails2 sensitive_purchased_fails_SS2 sensitive_purchased_fails_SR2
clear ratio_noprch2 fail_split_noprch2 fail_split_adjusted2 fail_split_noprch_adjusted2 ratio_error2 number_prch ratio_noprch_error2 number_noprch2

for drug = 1:params.number_drugs 
   
   total_prch(drug) = nnz(FTRSp.PCR_sameday(:,drug) & FTRSp.hasdiag);
 
  
   all_sensitive2 =  ismember(FTRSp.SMP_Res(:,drug),params.sensitive_group) & FTRSp.hasdiag ;
   all_resistant2 =  ismember(FTRSp.SMP_Res(:,drug),params.resistant_group)  & FTRSp.hasdiag ;
      
   all_sensitive_nextres =  FTRSp.next_res(:,drug) == 1 | FTRSp.SMP_Res(:,drug) == 2 ;
   all_resistant_nextres =  FTRSp.next_res(:,drug) == 3 ;
   all_nextres = ismember(FTRSp.next_res(:,drug),[1 2 3]);
   
   all_sensitive_purchased2 =  nnz((all_sensitive2 | all_resistant2) & FTRSp.PCR_sameday(:,drug));
   sensitive_purchased_fails2 = nnz(all_sensitive2 & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure);
   sensitive_purchased_nextresfails2 = nnz(all_sensitive2 & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure & all_nextres);
   
   sensitive_purchased_fails_SS = nnz(all_sensitive2 & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure & all_sensitive_nextres);
   sensitive_purchased_fails_SR = nnz(all_sensitive2 & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure & all_resistant_nextres);
   
   ratio2(drug) = sensitive_purchased_fails2/all_sensitive_purchased2*100;
   number_prch(drug) = all_sensitive_purchased2;
   ratio_error2(drug) = sqrt((ratio2(drug)/100).*(1-(ratio2(drug)/100))./(all_sensitive_purchased2))*100;
   fail_split2(1:2,drug) = [sensitive_purchased_fails_SS/sensitive_purchased_nextresfails2 ; sensitive_purchased_fails_SR/sensitive_purchased_nextresfails2]*ratio2(drug);
   SR_error2(drug) = sqrt((fail_split2(2,drug)/100).*(1-(fail_split2(2,drug)/100))./(all_sensitive_purchased2))*100;
   
   
   %all_resistant_purchased =  nnz(all_resistant & FTRSp.PCR_sameday(:,drug));
   all_resistant_purchased2 = nnz((all_sensitive2 | all_resistant2) & FTRSp.PCR_sameday(:,drug));
   resistant_purchased_fails2 = nnz(all_resistant2 & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure);
   resistant_purchased_nextresfails2 = nnz(all_resistant2 & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure & all_nextres);
   resistant_purchased_fails_RS = nnz(all_resistant2 & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure & all_sensitive_nextres);
   resistant_purchased_fails_RR = nnz(all_resistant2 & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure & all_resistant_nextres);
   
   ratio_res(drug) = resistant_purchased_fails2/all_resistant_purchased2*100;
   number_prch_res(drug) = all_resistant_purchased2;
   ratio_error_res(drug) = sqrt((ratio_res(drug)/100).*(1-(ratio_res(drug)/100))./(all_resistant_purchased2))*100;
   fail_split_res2(1:2,drug) = [resistant_purchased_fails_RR/resistant_purchased_nextresfails2 ; resistant_purchased_fails_RS/resistant_purchased_nextresfails2]*ratio_res(drug);
   SR_error_res2(drug) = sqrt((fail_split_res2(2,drug)/100).*(1-(fail_split_res2(2,drug)/100))./(all_resistant_purchased2))*100;
   
   
   %just_cleared(drug) = (100 -ratio(drug)) +(sensitive_purchased_fails_SS/sensitive_purchased_nextresfails)*ratio(drug);
   %fail_split_adjusted(1,drug) = ((sensitive_purchased_fails_SR/sensitive_purchased_nextresfails)*ratio(drug))/just_cleared(drug)*100;
   
   %all_sensitive_noprch =  nnz(all_sensitive & FTRSp.PCR_sameday(:,10));
   all_sensitive_noprch2 =  nnz((all_sensitive2 | all_resistant2) & FTRSp.PCR_sameday(:,10));
   sensitive_noprch_fails2 = nnz(all_sensitive2 & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure );
   sensitive_noprch_nextresfails2 = nnz(all_sensitive2 & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure & all_nextres);
   
   sensitive_noprch_fails_SS = nnz(all_sensitive2 & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure & all_sensitive_nextres);
   sensitive_noprch_fails_SR = nnz(all_sensitive2 & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure & all_resistant_nextres);
   
   ratio_noprch2(drug) = sensitive_noprch_fails2/all_sensitive_noprch2*100;
   number_noprch(drug) = all_sensitive_noprch2;
   ratio_noprch_error(drug) = sqrt((ratio_noprch2(drug)/100).*(1-(ratio_noprch2(drug)/100))./(all_sensitive_noprch2))*100;
   fail_split_noprch2(1:2,drug) = [sensitive_noprch_fails_SS/sensitive_noprch_nextresfails2 ; sensitive_noprch_fails_SR/sensitive_noprch_nextresfails2]*ratio_noprch2(drug);
   SR_noprch_error2(drug) = sqrt((fail_split_noprch2(2,drug)/100).*(1-(fail_split_noprch2(2,drug)/100))./(all_sensitive_noprch2))*100;
   
   
   %all_resistant_noprch =  nnz(all_resistant & FTRSp.PCR_sameday(:,10));
   all_resistant_noprch2 =  nnz((all_sensitive2 | all_resistant2) & FTRSp.PCR_sameday(:,10));
   resistant_noprch_fails2 = nnz(all_resistant2 & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure);
   resistant_noprch_nextresfails2 = nnz(all_resistant2 & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure & all_nextres);
   resistant_noprch_fails_RS = nnz(all_resistant2 & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure & all_sensitive_nextres);
   resistant_noprch_fails_RR = nnz(all_resistant2 & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure & all_resistant_nextres);
   
   ratio_noprch_res2(drug) = resistant_noprch_fails2/all_resistant_noprch2*100;
   number_noprch_res(drug) = all_resistant_noprch2;
   ratio_error_res(drug) = sqrt((ratio_noprch_res2(drug)/100).*(1-(ratio_noprch_res2(drug)/100))./(all_resistant_noprch2))*100;
   fail_split_noprch_res2(1:2,drug) = [resistant_noprch_fails_RR/resistant_noprch_nextresfails2 ; resistant_noprch_fails_RS/resistant_noprch_nextresfails2]*ratio_noprch_res2(drug);
   SR_error_noprch_res2(drug) = sqrt((fail_split_noprch_res2(2,drug)/100).*(1-(fail_split_noprch_res2(2,drug)/100))./(all_resistant_noprch2))*100;
   %just_cleared_noprch(drug) = (100 -ratio_noprch(drug)) +(sensitive_noprch_fails_SS/sensitive_noprch_nextresfails)*ratio_noprch(drug);
   %fail_split_noprch_adjusted(1,drug) = ((sensitive_noprch_fails_SR/sensitive_noprch_nextresfails)*ratio_noprch(drug))/just_cleared_noprch(drug)*100;
   
   %total_error(drug) = sqrt(P.*(1-P)./(n))
   
   clear  all_sensitive all_sensitive_purchased sensitive_purchased_fails
end

%%

width = 0.5;
gap = 1.8;
plot([0 7*gap+3], [0 0],'k')

for jj = 1:params.number_drugs
    ii = params.new_order(jj);
hold on
xpos = jj*gap-0.32;
hd3(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 fail_split_noprch2(2,ii) fail_split_noprch2(2,ii)],'w','EdgeColor','k', 'FaceColor', params.SR_color); % using MATLAB "patch" to establish the border
errorbar(xpos,fail_split_noprch2(2,ii)',SR_noprch_error2(ii), 'k', 'LineStyle','none')
%hatch(hda(ii),[45,4,1],'w')
hd4(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 -fail_split_noprch_res2(2,ii) -fail_split_noprch_res2(2,ii)],'w','EdgeColor','k', 'FaceColor', params.RS_color); % using MATLAB "patch" to establish the border
errorbar(xpos,-fail_split_noprch_res2(2,ii)',SR_error_noprch_res2(ii), 'k', 'LineStyle','none')


end

for jj = 1:params.number_drugs
    ii = params.new_order(jj);
hold on
xpos = jj*gap+0.32;
hd1(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 fail_split2(2,ii) fail_split2(2,ii)],'w','EdgeColor','k', 'FaceColor', params.SR_color); % using MATLAB "patch" to establish the border
errorbar(xpos,fail_split2(2,ii)',SR_error2(ii), 'k', 'LineStyle','none')

% temp_text = {'n = '; [num2str(round(total_prch(ii)/1000)) 'k']};
% text(xpos-width/2, fail_split(2,ii)+1, temp_text);

hd2(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 -fail_split_res2(2,ii) -fail_split_res2(2,ii)],'w','EdgeColor','k', 'FaceColor', params.RS_color); % using MATLAB "patch" to establish the border
errorbar(xpos,-fail_split_res2(2,ii)',SR_error_res2(ii), 'k', 'LineStyle','none');
hatch(hd3(jj),[45,2.5,0.8],'k')
hatch(hd4(jj),[45,2.5,0.8],'k')
end
xlim([0.5 7*gap+3])
ylim([-1.5 7.5])
ylabel('% of patients who have an infection which changes resistance')

set(gca,'XTick',gap:gap:gap*8);
set(gca,'xticklabel',FTRSp.SMP_Res_drug_names(params.new_order));
xtickangle(45);

end
