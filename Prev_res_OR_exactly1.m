function [] = Prev_res_OR_exactly1(UTI_cases,params)
number_of_drugs = 7;

X_age = UTI_cases.Demog.Age;
X_age(:,6) = [];
zeroto39 = sum(X_age(:,1:4),2);
eightyplus = sum(X_age(:,end-1:end),2);
X_demog = [zeroto39 X_age(:,5:end-2) eightyplus, UTI_cases.Demog.Gender, UTI_cases.Demog.Preg];

% figure;
% set(gcf,'color','w')
% clf; set(gcf,'name','Fig. S12 OR risk of emergence of resistance given exactly one previous sample.');
subplot(2,2,3)
%% OR specifically StoR failure given any previous resistance to drug
clear y_SSR p_SSR ciu_SSR cil_SSR se_SSR 

sen_res = 3;

seall = zeros(size(X_demog,2) +2, number_of_drugs);
yall = zeros(size(X_demog,2) +2, number_of_drugs);
pall = zeros(size(X_demog,2) +2, number_of_drugs);

for drug = 1:number_of_drugs

SMP_to_use =  find(UTI_cases.any_SRmeasurement_byyear(:,drug,1)==1 & (UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);
%SMP_to_use =  find(UTI_cases.any_SRmeasurement_356to730days(:,drug)==2 & (UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);

any_previous_R = UTI_cases.num_of_previous_R_byyear(SMP_to_use ,drug,1)>0 ;
%any_previous_R = UTI_cases.num_of_previous_R_356to730days(SMP_to_use ,drug)>0 ;
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;

drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);

X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1)];

Y = treatfailure(drug_prch_temp==1);

c = my_fit(X, Y, 50);

seall(:,drug) = c.se;
yall(:,drug) = c.coef;
pall(:, drug) = c.p;

end
%

%yall(isnan(yall)) = 0;
y_SSR(2:2:number_of_drugs*2) = exp(yall(end,:));
p_SSR(2:2:number_of_drugs*2) = pall(end,:);
se_SSR(2:2:number_of_drugs*2) = seall(end,:);
cil_SSR(2:2:number_of_drugs*2) = exp(yall(end,:) -  seall(end,:));
ciu_SSR(2:2:number_of_drugs*2) = exp(yall(end,:) +  seall(end,:));


%% OR specifically StoR failure given any previous resistance to drug
clear yall pall seall
sen_res = [1 2];


for drug = 1:number_of_drugs

SMP_to_use =  find(UTI_cases.any_SRmeasurement_byyear(:,drug,1)==1 & (UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);
%SMP_to_use =  find(UTI_cases.any_SRmeasurement_356to730days(:,drug)==2 & (UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);

any_previous_R = UTI_cases.num_of_previous_R_byyear(SMP_to_use ,drug,1)>0 ;
%any_previous_R = UTI_cases.num_of_previous_R_356to730days(SMP_to_use ,drug)>0 ;
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;

drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);

X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1)];

Y = treatfailure(drug_prch_temp==1);

c = my_fit(X, Y,50); 

seall(:,drug) = c.se;
yall(:,drug) = c.coef;
pall(:, drug) = c.p;

end
%

%yall(isnan(yall)) = 0;
%y_SR = exp(yall(end,:));
%p_SR = pall(end,:);
y_SSR(1:2:number_of_drugs*2-1) = exp(yall(end,:));
%y_SSR(1:2:number_of_drugs_7*2-1) = (yall(end,:));
p_SSR(1:2:number_of_drugs*2-1) = pall(end,:);
se_SSR(1:2:number_of_drugs*2-1) = seall(end,:);

cil_SSR(1:2:number_of_drugs*2-1) = exp(yall(end,:) -  seall(end,:));
ciu_SSR(1:2:number_of_drugs*2-1) = exp(yall(end,:) +  seall(end,:));

b1 = bar(y_SSR(1:number_of_drugs*2), 'FaceColor','flat',  'BarWidth', 0.8);
hold on

b1(1).CData(1:2:number_of_drugs*2-1,:) = ones(number_of_drugs,1)*params.SS_color;
b1.CData(2:2:number_of_drugs*2,:) = ones(number_of_drugs,1)*params.SR_color;
%plot(b1.XData,cil_SSR,'kx');
%plot(b1.XData,ciu_SSR,'kx');
errneg = y_SSR - cil_SSR;
errpos =  ciu_SSR -y_SSR;

errorbar(b1.XData,y_SSR,errneg,errpos,'k', 'LineStyle','none');

hold on
for ii = 1:length(y_SSR)
    if p_SSR(ii)< 0.01 && p_SSR(ii)> 0.001
    plot(b1.XData(ii), y_SSR(ii)+1, '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.001 && p_SSR(ii)> 0.00001
    plot([b1.XData(ii)-0.1 b1.XData(ii)+0.1] , [y_SSR(ii)+1 y_SSR(ii)+1], '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.00001
        plot([b1.XData(ii)-0.2  b1.XData(ii) b1.XData(ii)+0.2] , [ y_SSR(ii)+1  y_SSR(ii)+1 y_SSR(ii)+1], '*k', 'MarkerSize',4)
    end
end
% c.se(c.p> 0.05) = 0;
% hold on
% errorbar(c.coef(2:end),c.se(2:end),'.')
xticks((1:2:number_of_drugs*2)+0.5);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(1:number_of_drugs));
xtickangle(45);

ylabel({'Adjusted odds ratio of early recurrence'; ...
    'given a prior resitant sample in the past year'}, 'FontSize' , 8);

text(0,6.5,'{\it * P<0.01}');
text(0,6,'{\it** P<0.001}');
text(0,5.5,'{\it*** P<0.00001}');

set(gca, 'YScale', 'log');
yticks(1:8)
ylim([0.4 8.5])

end
