function [] = rate_res_by_species(UTI_cases)

Names = UTI_cases.SMP_Res_drug_names;
Names{end+1} = 'All purchases';
Names{contains(Names,'Cefuroxime')} = 'Cefuroxime/Axetil';
Names{8} = 'Ofloxacin';
number_drugs =8 ;


all_bugs = cell2mat(UTI_cases.bug_all);
all_res= cell2mat(UTI_cases.RES);

%% by species
bug = cell2mat(UTI_cases.bug_all);
SMP_Res = cell2mat(UTI_cases.RES);
frac_ecoli = nnz(bug == 179)/length(bug)*100;
frac_kleb = nnz(bug == 225)/length(bug)*100;
frac_prot= nnz(bug == 280)/length(bug)*100;
frac_cit = nnz(bug == 110 )/length(bug)*100;
number_drugs = 8;
Names = UTI_cases.SMP_Res_drug_names;
pos2 = [0 0.05 0.1 0.1] ;%[left bottom width height]
for kk = [110   280  225 179 ]
for drug = 1:number_drugs
pos2(1) = pos2(1)+0.1;
subplot('Position',pos2)
yes_measurement = SMP_Res(:,drug)>0 & bug == kk;
number_resistant(drug) = nnz(ismember(SMP_Res(yes_measurement,drug),3));
number_intermediate(drug) = nnz(ismember(SMP_Res(yes_measurement,drug),[2]));
number_sensitive(drug) = nnz(ismember(SMP_Res(yes_measurement,drug),[1]));

fractionR(drug) = (number_resistant(drug)/(number_resistant(drug) +  number_intermediate(drug) + number_sensitive(drug)));
fractionI(drug) = (number_intermediate(drug)/(number_resistant(drug) + number_intermediate(drug) + number_sensitive(drug)));
fractionS(drug) = (number_sensitive(drug)/(number_resistant(drug) + number_intermediate(drug) + number_sensitive(drug)));

pie([fractionS(drug),fractionI(drug), fractionR(drug)]);
set(findobj(gca,'type','text'),'fontsize',7.5);
%labels = repmat({''},size(fractionS));
%title(Names{drug})
cmap = [150 150 150 ; 75 75 75 ; 0 0 0]./255;
colormap(gca, cmap)
end
text(0, pos2(2), UTI_cases.Bugs.Name(find(UTI_cases.Bugs.Code == kk)));
pos2 = [0 pos2(2) 0.1 0.1];
pos2(2) = pos2(2)+0.2;
end

%%
for drug = 1:number_drugs
pos2(1) = pos2(1)+0.1;
subplot('Position',pos2);
yes_measurement = SMP_Res(:,drug)>0;
number_resistant(drug) = nnz(ismember(SMP_Res(yes_measurement,drug),3));
number_intermediate(drug) = nnz(ismember(SMP_Res(yes_measurement,drug),[2]));
number_sensitive(drug) = nnz(ismember(SMP_Res(yes_measurement,drug),[1]));
%fraction(drug) = (number_resistant(drug)/(number_resistant(drug) + number_sensitive(drug)));
fractionR(drug) = (number_resistant(drug)/(number_resistant(drug) +  number_intermediate(drug) + number_sensitive(drug)));
fractionI(drug) = (number_intermediate(drug)/(number_resistant(drug) + number_intermediate(drug) + number_sensitive(drug)));
fractionS(drug) = (number_sensitive(drug)/(number_resistant(drug) + number_intermediate(drug) + number_sensitive(drug)));
pie([fractionS(drug),fractionI(drug), fractionR(drug)]);
set(findobj(gca,'type','text'),'fontsize',7.5);
%labels = repmat({''},size(fractionS));
title(Names{drug})
cmap = [150 150 150 ; 75 75 75 ; 0 0 0]./255;
colormap(gca, cmap)
end
text(0, pos2(2), 'All species');
end