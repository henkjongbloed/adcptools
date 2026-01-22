function plotWL_CSA(F, tak, BN)
% int = [datenum('11-Aug-2014 00:00:00'), datenum('13-Aug-2014 23:50:00')];
hours = 24*(F{1}.time - F{1}.time(1));
figure;
st = {'-', '--', '-.'};
yyaxis left    
plot(hours , F{1}.eta, 'k-')
for b = 1:length(tak)    
    hold on
    hours = 24*(F{b}.time - F{b}.time(1));
    yyaxis left
    plot(hours, F{b}.uc{1,1} + F{b}.uc{3,1}, ['b',st{b}]) % CSA velocity m/s
    yyaxis right
    plot(hours, F{b}.sc{1,1} + F{b}.sc{3,1}, ['r',st{b}]) % CSA salt (psu)
end
xlabel('t')
yyaxis left
ylabel('m or m/s')
yyaxis right
ylabel('psu')
leg = {'eta'};
j=1;
for i = 2:1:(length(tak)+1)
    leg{i} = ['u: ', BN{j}];
    leg{i+length(tak)} = ['s: ', BN{j}];
    j = j+ 1;
end
legend(leg)
title('Waterlevel, cross sectional averaged flow and salinity')
grid on


end