function []= risk_reccurence_correctly_prescribed_treat_v_untreat(FTRSp,params)

clear ratio fail_split all_sensitive all_sensitive_purchased sensitive_purchased_fails sensitive_purchased_fails_SS sensitive_purchased_fails_SR
clear ratio_noprch fail_split_noprch fail_split_adjusted fail_split_noprch_adjusted ratio_error number_prch ratio_noprch_error number_noprch
clear  all_sensitive all_sensitive_purchased sensitive_purchased_fails
for drug = 1:params.number_drugs  
 
   total_prch(drug) = nnz(FTRSp.PCR_sameday(:,drug) & FTRSp.hasdiag);
   all_sensitive = (ismember(FTRSp.SMP_Res(:,drug), params.sensitive_group)) & FTRSp.hasdiag ;   
   all_sensitive_nextres =  ismember(FTRSp.next_res(:,drug),params.sensitive_group);
   all_resistant_nextres =  ismember(FTRSp.next_res(:,drug),params.resistant_group);
   all_nextres = ismember(FTRSp.next_res(:,drug),[1 2 3]);
   
   all_sensitive_purchased(drug) =  nnz(all_sensitive & FTRSp.PCR_sameday(:,drug));
   sensitive_purchased_fails(drug) = nnz(all_sensitive & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure);
   sensitive_purchased_nextresfails(drug) = nnz(all_sensitive & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure & all_nextres);
   
   sensitive_purchased_fails_SS(drug) = nnz(all_sensitive & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure & all_sensitive_nextres);
   sensitive_purchased_fails_SR(drug) = nnz(all_sensitive & FTRSp.PCR_sameday(:,drug) & FTRSp.treatfailure & all_resistant_nextres);
   
   ratio(drug) = sensitive_purchased_fails(drug)/all_sensitive_purchased(drug)*100;
   number_prch(drug) = all_sensitive_purchased(drug);
   ratio_error(drug) = sqrt((ratio(drug)/100).*(1-(ratio(drug)/100))./(all_sensitive_purchased(drug)))*100;
   fail_split(1:2,drug) = [sensitive_purchased_fails_SS(drug)/sensitive_purchased_nextresfails(drug) ; sensitive_purchased_fails_SR(drug)/sensitive_purchased_nextresfails(drug)]*ratio(drug);
   
     
   all_sensitive_noprch =  nnz(all_sensitive & FTRSp.PCR_sameday(:,10));
   sensitive_noprch_fails = nnz(all_sensitive & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure );
   sensitive_noprch_nextresfails = nnz(all_sensitive & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure & all_nextres);
   sensitive_noprch_fails_SS = nnz(all_sensitive & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure & all_sensitive_nextres);
   sensitive_noprch_fails_SR = nnz(all_sensitive & FTRSp.PCR_sameday(:,10) & FTRSp.treatfailure & all_resistant_nextres);
   
   ratio_noprch(drug) = sensitive_noprch_fails/all_sensitive_noprch*100;
   number_noprch(drug) = all_sensitive_noprch;
   ratio_noprch_error(drug) = sqrt((ratio_noprch(drug)/100).*(1-(ratio_noprch(drug)/100))./(all_sensitive_noprch))*100;

   fail_split_noprch(1:2,drug) = [sensitive_noprch_fails_SS/sensitive_noprch_nextresfails ; sensitive_noprch_fails_SR/sensitive_noprch_nextresfails]*ratio_noprch(drug);
   
end


width = 0.5;
gap = 1.8;
plot([0 7*gap+3], [0 0],'k')

for jj = 1:params.number_drugs
    ii = params.new_order(jj);
    hold on
xpos = jj*gap-0.32;
hd3(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 fail_split_noprch(2,ii) fail_split_noprch(2,ii)],'w','EdgeColor','k', 'FaceColor', params.SR_color); % using MATLAB "patch" to establish the border
%hatch(hda(ii),[45,4,1],'w')
hd4(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[fail_split_noprch(2,ii) fail_split_noprch(2,ii) (fail_split_noprch(1,ii)+fail_split_noprch(2,ii)) (fail_split_noprch(1,ii)+fail_split_noprch(2,ii))],'w','EdgeColor','k', 'FaceColor', params.SS_color); % using MATLAB "patch" to establish the border
%errorbar(xpos,-fail_split_noprch_res2(2,ii)',SR_error_noprch_res2(ii), 'k', 'LineStyle','none')
errorbar(xpos,(fail_split_noprch(1,ii)+fail_split_noprch(2,ii)),ratio_noprch_error(ii), 'k', 'LineStyle','none')

end

for jj = 1:params.number_drugs
    ii = params.new_order(jj);
hold on
xpos = jj*gap+0.32;
hd1(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 fail_split(2,ii) fail_split(2,ii)],'w','EdgeColor','k', 'FaceColor', params.SR_color); % using MATLAB "patch" to establish the border
% temp_text = {'n = '; [num2str(round(number_prch(ii)/1000)) 'k']};
% text(xpos-width/2, (fail_split(1,ii)+fail_split(2,ii))+2, temp_text);
hd2(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[fail_split(2,ii)  fail_split(2,ii) (fail_split(1,ii)+fail_split(2,ii)) (fail_split(1,ii)+fail_split(2,ii))],'w','EdgeColor','k', 'FaceColor', params.SS_color); % using MATLAB "patch" to establish the border
%errorbar(xpos,-fail_split_res2(2,ii)',SR_error_res2(ii), 'k', 'LineStyle','none');
hatch(hd3(jj),[45,2.5,0.8],'k')
hatch(hd4(jj),[45,2.5,0.8],'w')
errorbar(xpos,(fail_split(1,ii)+fail_split(2,ii)),ratio_error(ii), 'k', 'LineStyle','none')
temp_text = {'n = '; [num2str(round(number_prch(ii)/1000)) 'k']};
text(xpos-width/2, (fail_split(1,ii)+fail_split(2,ii))+3, temp_text);
end
xlim([0.5 7*gap+3])
ylim([0 20])
ylabel('% of treatment failures ')
set(gca,'XTick',gap:gap:gap*8);
set(gca,'xticklabel',FTRSp.SMP_Res_drug_names(params.new_order));
xtickangle(45);

end