% coord_example

RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
cd([RF,'WP2\TwoJunctions\code'])
addpath(genpath(strcat(RF,'Tools\adcptools'))); %path to ADCPTools
addpath(strcat(RF,'Tools\loess-master')); %path to loess-master
% addpath(strcat(RF,'Tools\T_Tide'))
addpath(strcat(RF,'Tools\BrewerMap-master')); %path to BrewerMap (see github)
% addpath(strcat(RF,'Tools\subaxis'));

F2 = 'WP2\TwoJunctions\';
addpath(strcat(RF,F2,'\data\processedData')); %path to .mat structs of processed data
addpath(strcat(RF,'Tools\plotting\plotting'));



lim = {[0,1]; [0, 3]; [0,1]};
res = [10, 31, 21];
X = get_coords(lim, res);
h = ones(size(X.y)); h(11:20) = 2;

zb = ones(size(h));
zb(11:20) = 0;
wl = 2*ones(size(X.t));
[Wl, Zb] = ndgrid(wl, zb);
H = Wl - Zb;
D = Decomposition(X=X, H=H, wl=wl, zb=zb);

% Generate F

F = -ones(size(X.T));
F(:, 11:20, 1:10) = 1;
[DF, AF] = D.decompose_function(F);






D.plot_components(DF, "velmap")
% D.plot_components(DF, "velmap")


FF=D.get_prod_components_all(DF, DF); % Salt Flux
D.plotUSF(FF, 0:7)