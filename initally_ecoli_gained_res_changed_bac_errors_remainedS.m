function [] = initally_ecoli_gained_res_changed_bac_errors_remainedS(FTRSp,params)
%all recurrent cases of initially E coli
%map on changed strain
%and numbers from sequening

%%
dates_to_use_start([1:4 6 8]) = min(FTRSp.SamplingDate);
dates_to_use_start(5) = min(FTRSp.SamplingDate)+7*321;
dates_to_use_start(7) = min(FTRSp.SamplingDate)+7*293;

dates_to_use_end(1:7) = max(FTRSp.SamplingDate);
dates_to_use_end(8) = min(FTRSp.SamplingDate)+7*293;
%%

number_drugs = 8;
%SMP.changed_bac = double([diff(SMP.bug);0]~=0);
%SMP.changed_bac(SMP.same_patient == 0) = NaN;

for drug = 1:number_drugs 
 
   dates_index =  find(FTRSp.SamplingDate >= dates_to_use_start(drug) & FTRSp.SamplingDate <= dates_to_use_end(drug));
   total_prch(drug) = nnz(FTRSp.PCR_sameday(dates_index,drug) & FTRSp.hasdiag(dates_index));
 
   all_sensitive = (FTRSp.SMP_Res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2) & FTRSp.hasdiag(dates_index) & FTRSp.bug(dates_index) == 179 & FTRSp.nBac(dates_index) == 1;
      
   all_sensitive_nextres =  FTRSp.next_res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2 ;
   all_resistant_nextres =  FTRSp.next_res(dates_index,drug) == 3 ;
   %all_resistant_nextres =  FTRSp.next_res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2 ;
   
   all_nextres = ismember(FTRSp.next_res(dates_index,drug),[3]);
   
   sensitive_purchased_fails = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres);
   %sensitive_purchased_fails = nnz(all_sensitive & FTRSp.PCR_sameday(:,drug) );
   
   sensitive_purchased_fails_stayedEcoli = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres & FTRSp.changed_bac(dates_index)==0);
   sensitive_purchased_fails_changed_kleb = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres & FTRSp.changed_bac(dates_index)==1 & FTRSp.new_bac(dates_index)==225);
   sensitive_purchased_fails_changed_prot = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres & FTRSp.changed_bac(dates_index)==1 & FTRSp.new_bac(dates_index)==280);
   sensitive_purchased_fails_changed_other = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres & FTRSp.changed_bac(dates_index)==1 & FTRSp.new_bac(dates_index)~=280 & FTRSp.new_bac(dates_index)~=225 & floor(FTRSp.new_bac(dates_index))~=179 & FTRSp.new_bac(dates_index)>0);
   
   ratio_stayed(drug) = sensitive_purchased_fails_stayedEcoli/sensitive_purchased_fails*100;
   ratio_kleb(drug) = sensitive_purchased_fails_changed_kleb/sensitive_purchased_fails*100;
   ratio_prot(drug) = sensitive_purchased_fails_changed_prot/sensitive_purchased_fails*100;
   ratio_other(drug) = sensitive_purchased_fails_changed_other/sensitive_purchased_fails*100;
   
   ratio_changed(drug) =  ratio_kleb(drug) + ratio_prot(drug) + ratio_other(drug);
   ratio_error(drug) = sqrt((ratio_changed(drug)/100).*(1-(ratio_changed(drug)/100))./(sensitive_purchased_fails))*100;
   
   
   number_prch(drug) = sensitive_purchased_fails;
%   ratio_error(drug) = sqrt((ratio(drug)/100).*(1-(ratio(drug)/100))./(sensitive_purchased_fails))*100;
   
   clear  all_sensitive all_sensitive_purchased sensitive_purchased_fails
end

%%
total = sum([ratio_stayed; ratio_kleb; ratio_prot; ratio_other]);
ratio_stayed = ratio_stayed./total*100;
ratio_kleb = ratio_kleb./total*100;
ratio_prot = ratio_prot./total*100;
ratio_other = ratio_other./total*100;
sum([ratio_stayed(1:7); ratio_kleb(1:7); ratio_prot(1:7); ratio_other(1:7)])
%sum([ratio_stayed_no_prch(1:7); ratio_kleb_no_prch(1:7); ratio_prot_no_prch(1:7); ratio_other_no_prch(1:7)])

%%

gainedR = [ratio_other; ratio_prot; ratio_kleb];

%% remained sensitve 
for drug = 1:number_drugs 
 
   dates_index =  find(FTRSp.SamplingDate >= dates_to_use_start(drug) & FTRSp.SamplingDate <= dates_to_use_end(drug));
   total_prch(drug) = nnz(FTRSp.PCR_sameday(dates_index,drug) & FTRSp.hasdiag(dates_index));
 
   all_sensitive = (FTRSp.SMP_Res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2) & FTRSp.hasdiag(dates_index) & FTRSp.bug(dates_index) == 179 & FTRSp.nBac(dates_index) == 1;
      
   %all_sensitive_nextres =  FTRSp.next_res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2 ;
   %all_resistant_nextres =  FTRSp.next_res(dates_index,drug) == 3 ;
   all_resistant_nextres =  FTRSp.next_res(dates_index,drug) == 1 | FTRSp.SMP_Res(dates_index,drug) == 2 ;
   
   all_nextres = ismember(FTRSp.next_res(dates_index,drug),[3]);
   
   sensitive_purchased_fails = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres);
   %sensitive_purchased_fails = nnz(all_sensitive & FTRSp.PCR_sameday(:,drug) );
   
   sensitive_purchased_fails_stayedEcoli = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres & FTRSp.changed_bac(dates_index)==0);
   sensitive_purchased_fails_changed_kleb = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres & FTRSp.changed_bac(dates_index)==1 & FTRSp.new_bac(dates_index)==225);
   sensitive_purchased_fails_changed_prot = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres & FTRSp.changed_bac(dates_index)==1 & FTRSp.new_bac(dates_index)==280);
   sensitive_purchased_fails_changed_other = nnz(all_sensitive & FTRSp.PCR_sameday(dates_index,drug) & FTRSp.treatfailure(dates_index) & all_resistant_nextres & FTRSp.changed_bac(dates_index)==1 & FTRSp.new_bac(dates_index)~=280 & FTRSp.new_bac(dates_index)~=225 & floor(FTRSp.new_bac(dates_index))~=179 & FTRSp.new_bac(dates_index)>0);
   
   ratio_stayed(drug) = sensitive_purchased_fails_stayedEcoli/sensitive_purchased_fails*100;
   ratio_kleb(drug) = sensitive_purchased_fails_changed_kleb/sensitive_purchased_fails*100;
   ratio_prot(drug) = sensitive_purchased_fails_changed_prot/sensitive_purchased_fails*100;
   ratio_other(drug) = sensitive_purchased_fails_changed_other/sensitive_purchased_fails*100;
   
   ratio_changed(drug) =  ratio_kleb(drug) + ratio_prot(drug) + ratio_other(drug);
   ratio_error2(drug) = sqrt((ratio_changed(drug)/100).*(1-(ratio_changed(drug)/100))./(sensitive_purchased_fails))*100;
   
   
   number_prch(drug) = sensitive_purchased_fails;
%   ratio_error(drug) = sqrt((ratio(drug)/100).*(1-(ratio(drug)/100))./(sensitive_purchased_fails))*100;
   
   clear  all_sensitive all_sensitive_purchased sensitive_purchased_fails
end
%%
total = sum([ratio_stayed; ratio_kleb; ratio_prot; ratio_other]);
ratio_stayed = ratio_stayed./total*100;
ratio_kleb = ratio_kleb./total*100;
ratio_prot = ratio_prot./total*100;
ratio_other = ratio_other./total*100;
sum([ratio_stayed(1:7); ratio_kleb(1:7); ratio_prot(1:7); ratio_other(1:7)]);
%sum([ratio_stayed_no_prch(1:7); ratio_kleb_no_prch(1:7); ratio_prot_no_prch(1:7); ratio_other_no_prch(1:7)])

%%

remainedS = [ratio_other; ratio_prot; ratio_kleb];
%reorder for fig
gainedR = gainedR(:,[1 2 8 3:7]);
remainedS =remainedS(:,[1 2 8 3:7]);

% width = 0.5;
% gap = 1;
width = 0.5;
gap = 1.8;
plot([0 7*gap+3], [0 0],'k')

new_kleb_col = [0.5, 0.5, 0.5];
new_prot_col = [0.8, 0.8, 0.8];
new_other_col = [0.2, 0.2, 0.2];
for ii = 1:8
hold on
xpos = ii*gap-0.32;
hd1(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 remainedS(1,ii) remainedS(1,ii)],'w','EdgeColor','k', 'FaceColor', new_other_col); % using MATLAB "patch" to establish the border
hd2(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[remainedS(1,ii)  remainedS(1,ii) (remainedS(1,ii)+remainedS(2,ii)) (remainedS(1,ii)+remainedS(2,ii))],'w','EdgeColor','k', 'FaceColor', new_prot_col); % using MATLAB "patch" to establish the border
hd3(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[sum(remainedS(1:2,ii)) sum(remainedS(1:2,ii)) sum(remainedS(1:3,ii)) sum(remainedS(1:3,ii))],'w','EdgeColor','k', 'FaceColor', new_kleb_col ); % using MATLAB "patch" to establish the border
errorbar(xpos,sum(remainedS(1:3,ii))',ratio_error2(ii), 'k', 'LineStyle','none')

end
for ii = 1:8
hold on
xpos = ii*gap+0.32;
hd1(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 gainedR(1,ii) gainedR(1,ii)],'w','EdgeColor','k', 'FaceColor', new_other_col); % using MATLAB "patch" to establish the border
hd2(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[gainedR(1,ii)  gainedR(1,ii) (gainedR(1,ii)+gainedR(2,ii)) (gainedR(1,ii)+gainedR(2,ii))],'w','EdgeColor','k', 'FaceColor', new_prot_col); % using MATLAB "patch" to establish the border
hd3(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[sum(gainedR(1:2,ii)) sum(gainedR(1:2,ii)) sum(gainedR(1:3,ii)) sum(gainedR(1:3,ii))],'w','EdgeColor','k', 'FaceColor', new_kleb_col ); % using MATLAB "patch" to establish the border
errorbar(xpos,sum(gainedR(1:3,ii))',ratio_error(ii), 'k', 'LineStyle','none')

end

xlim([0.5 7*gap+3])
ylim([0 80])
% xlim([0 8.5])
% ylim([0 80])
ylabel({'% of initially {\it E. coli} infections which'; ...
    'gained resistance that were caused by';...
    'reinfection different species'})
xlabel('treatment antibiotic')
set(gca,'XTick',width:gap:gap*8);
set(gca,'xticklabel', FTRSp.SMP_Res_drug_names([1 2 8 3:7]));
xtickangle(45);
legend('Other species','{\it P. mirabilis}','{\it K. pneumoniae}')
