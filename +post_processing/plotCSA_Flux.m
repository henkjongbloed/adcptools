function plotCSA_Flux(tak, F, U)
XTL = {'F_{00}', 'F_{01}^*', 'F_{02}^*', 'F_{03}^*', 'F_{04}^*',...
    'F_{10}', 'F_{11}^*',...
    'F_{20}', 'F_{21}^*',...
    'F_{30}', 'F_{31}^*', 'F_{32}^*', 'F_{33}^*', 'F_{34}^*', 'F_{35}', 'F_{36}^*',...
    'F_{40}', 'F_{41}', 'F_{42}'};

XTLT = {'F_{00}', 'F_{01}^*',...
    'F_{10}', 'F_{11}^*',...
    'F_{20}', 'F_{21}^*',...
    'F_{30}', 'F_{31}^*',...
    'F_{40}'};

XTLAbs = {'F_{00}', 'F_{01}^*',...
    'F_{10}', 'F_{11}^*',...
    'F_{20}', 'F_{21}^*',...
    'F_{30}', 'F_{31}^*',...
    'F_{40}'};


XTLDuo = {'F_{hom}', 'F^*'};


figure;
dim = 1;
for j = 1:length(tak)
    for d = dim
        subplot(length(tak), length(dim),  sub2ind([length(tak), length(dim)], j, d))
        bar(F{j}.LatT(:,d))
        title(strcat('Lateral Decomposition at ', U{j}.BN))
        xticks(1:length(XTLT))
        xticklabels(XTLT)
        grid on
    end
end


figure;
dim = 1;
for j = 1:length(tak)
    for d = dim
        subplot(length(tak), length(dim),  sub2ind([length(tak), length(dim)], j, d))
        bar(F{j}.Lat(:,d))
        title(strcat('Lateral Decomposition at ', U{j}.BN))
        xticks(1:length(XTL))
        xticklabels(XTL)
        grid on
    end
end

figure;
dim = 1;
for j = 1:length(tak)
    for d = dim
        subplot(length(tak), length(dim),  sub2ind([length(tak), length(dim)], j, d))
        bar(F{j}.LatAbs(:,d))
        title(strcat('Abs Lateral Decomposition at ', U{j}.BN))
        xticks(1:length(XTLAbs))
        xticklabels(XTLAbs)
        grid on
    end
end

figure;
dim = 1;
LD = [];
BN = [];
for j = 1:length(tak)
    LD = [LD; F{j}.LatDuo];
    BN{j} = U{j}.BN;
end
BN = categorical(BN);
b=bar(BN, LD, 'stacked');
title(strcat('Uniform vs. laterally varying') )

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

figure;
dim = 1;
LD = [];
BN = [];
for j = 1:length(tak)
    LD = [LD; F{j}.LatDuoAbs];
    BN{j} = U{j}.BN;
end
BN = categorical(BN);
b=bar(BN, LD, 'stacked');
title(strcat('Uniform vs. laterally varying') )

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

% xticks(1:length(XTLDuo))
% xticklabels(XTLDuo)

%grid on

end