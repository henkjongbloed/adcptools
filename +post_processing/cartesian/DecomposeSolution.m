% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
      set(groot, default_name,'latex');
end  

% Scaling of coordinates
T = 1;
B = 1;
t0 = 0;

Nt = 20;
Ny = 20;
Ns = 20;

t = linspace(t0, T, Nt);
y = linspace(0, B, Ny);
s = linspace(0, 1, Ns);

X = {t,y,s}; % Coordinates

[Tt, Yy, Ss] = ndgrid(t, y, s);
phi = pi/3;

F = sin(2*pi*Tt).*(Ss).^2;%0*rand(size(Tt))-1/2; %all variables are of order 1
G = sin(2*pi*Tt - phi).*(1-Ss).^2;% 0*rand(size(Tt))-1/2;% - .1*rand(size(Tt));
H = 10*ones(size(Tt)) + sin(2*pi*Tt/T + phi);




% sz = size(F);
% taH = ta(H);
% waH = wa(H);
% 
% res = twa(F, H, sz, B);
% res2 = la(res, H, sz, B);
% res3 = va(res2, H, sz, B);

%disp(squeeze(res3))
% avg_op = {@twa, @va};
% plot_component(F, avg_op, H, B)
% 
%avg_op = {@va, @twa};
%plot_component(F, avg_op, H, B)

% Why cross-terms not very very small - difference mean and trapz?

[DF, AF, DFf, AFf] = decompose_function(F, H, X);
[DG, AG, DGf, AGf] = decompose_function(G, H, X);
%avg_0(F, H, B)
C = zeros(numel(DFf), numel(DFf));
SF = zeros([numel(DFf),1]);
for i = 1:numel(DFf)
    for j = 1:numel(DFf)
        C(i,j) = avg_0(DFf{i}.*DGf{j}, H, X);
    end
    SF(i) = C(i,i);
end
%plot_components(AF, X)
plot_components(DF, X)

%plot_components(AG, X)
plot_components(DG, X)

plotC(C);
figure;
bar(SF)
disp(sum(SF))
disp(avg_0(F.*G, H, X));



% Elementary functions - trapz or mean
function taf = ta(f,X)
%taf = mean(f, 1);
%taf = sum(f, 1)/(size(f,1)-1);
taf = trapz(X{1}, f, 1);
end 

function waf = wa(f, X)
%waf = mean(f, 2);
%waf = sum(f, 2)/(size(f,2)-1);
waf = trapz(X{2}, f, 2);
end

function daf = da(f, X)
%daf = mean(f,3);
%daf = sum(f, 3)/(size(f,3)-1);
daf = trapz(X{3}, f, 3);
end

function A = cs(H, X)
A = X{2}(end)*wa(H, X); % Nt by 1 vector
end

% Thickness-weighed functions
function f = repmat3(f, sz)
% Reshape matrix to size sz , repeating mean values if needed
rep = [1 1 1];
if numel(size(f)) == 2
    rep(3) = sz(3);
end
rep([size(f)==1]) = sz([size(f)==1]);
if ~isempty(rep)
    f = repmat(f, rep);
end
end


function [twaf1, twaf2] = twa(f, H, X, sz)
% Taking thickness-weighed average over t
% Could be vectorized but for now I prefer readability
% Returns squeezed AND full array

[f, ~, sf] = handle_input_dim(f, H, sz);

if sf(2) == 1 % Laterally averaged variable
    A = cs(H, X);
    twaf0 = ta(A.*f, X)./ta(A, X);
else
    twaf0 = ta(H.*f, X)./ta(H, X); % Just a factor B difference. Yes.
end
twaf1 = twaf0(1, 1:sf(2), 1:sf(3));
twaf2 = repmat3(twaf0, sz);
end

function [laf1, laf2] = la(f, H, X, sz)
% Returns squeezed AND full array
[f, ~, sf] = handle_input_dim(f, H, sz);
if sf(1) == 1 % Time-averaged variable
    A = cs(H, X);
    laf0 = X{2}(end)*wa(ta(H, X).*f, X)./ta(A, X);
else
    laf0 = wa(H.*f, X)./wa(H, X);
end
laf1 = laf0(1:sf(1), 1, 1:sf(3));
laf2 = repmat3(laf0, sz);
end



function [vaf1, vaf2] = va(f, H, X, sz)
% Returns squeezed AND full array
[f, ~, sf] = handle_input_dim(f, H, sz);
vaf0 = da(f, X);
vaf1 = vaf0(1:sf(1), 1:sf(2), 1);
vaf2 = repmat3(vaf0, sz);
end

function [f, H, sf] = handle_input_dim(f, H, sz)
sf = size(f); % Get original size of input data
if size(f,3)==1
    sf = [sf 1];
end
f = repmat3(f, sz); 
H = repmat3(H, sz);
end

function plot_components(AF, X)
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:numel(AF)
    subplot(2,4,i)
    plot_component(AF{i}, X)
end
end

function plot_component(af, X)
sf = size(af);
if size(af,3)==1
    sf = [sf 1];
end
if sf(1) == 1 % t-independent
   if sf(2) == 1 % y-independent
       if sf(3) == 1 % f0
           bar(af)                                      %1
           xticklabels({"f0"})
       else % only sig-dependent
            plot(squeeze(af), linspace(0,1,numel(af)))    %2
            xlabel("f(sig)")
            ylabel("sig")
       end
   else % y-dependent
       if sf(3) == 1 % only y-dependent
           plot(linspace(0,X{2}(end),numel(af)), squeeze(af))     %3
            xlabel("y")
            ylabel("f(y)")
       else % y and sig dependent
           contourf(squeeze(af))                        %4
           xlabel("y")
           ylabel("sig")
       end
   end
else % t-dependent
   if sf(2) == 1 % y-independent
       if sf(3) == 1 % only t-dependent
           plot(squeeze(af))
                      xlabel("t")
           ylabel("f(t)")%5
       else % only (t,sig)-dependent
                      xlabel("t")
           ylabel("sig")                        %6
       end
   else % y-dependent
       if sf(3) == 1 % only (t,y)-dependent
           contourf(squeeze(af))
                      xlabel("t")
           ylabel("y")%7
       else % t and y and sig dependent
           bar(0)
           disp('Can only display at most 2D variables')
       end
   end
end
end

function [f,ff] = avg_comp(f, avg_op, H, X)
% Arbitrary composition of averaging operators
sf0 = size(f);
if size(f,3)==1
    sf0 = [sf0 1];
end
%f0 = f;
ff = f;
for idx = 1:numel(avg_op)
    [f, ~] = avg_op{idx}(f, H, X, sf0);
    [~, ff] = avg_op{idx}(ff, H, X, sf0);
end
end

function [f, ff] = avg_0(f, H, X)
avg_op = {@twa, @la, @va};
[f,ff] = avg_comp(f, avg_op, H, X);
end

function [DF, AF, DFf, AFf] = decompose_function(f, H, X)

AF{1} = f;
AFf{1}  = f;

[AF{2}, AFf{2}] = avg_comp(f, {@twa}, H, X);    %twa(f)
[AF{3}, AFf{3}] = avg_comp(f, {@la}, H, X);     %la(f)  
[AF{4}, AFf{4}] = avg_comp(f, {@va}, H, X);     %va(f)

[AF{5}, AFf{5}] = avg_comp(f, {@twa, @la}, H, X);   %twa(la(f))
[AF{6}, AFf{6}] = avg_comp(f, {@twa, @va}, H, X);   %twa(va(f))
[AF{7}, AFf{7}] = avg_comp(f, {@la, @va}, H, X);    %la(va(f))

[AF{8}, AFf{8}] = avg_comp(f, {@twa, @la, @va}, H, X);  %f_0


DFf{1} = AFf{8};                                DF{1} = AF{8};

DFf{2} = AFf{7} - DFf{1};   	                DF{2} = DFf{2}(:,1,1);
DFf{3} = AFf{6} - DFf{1};                       DF{3} = DFf{3}(1,:,1);
DFf{4} = AFf{5} - DFf{1};                       DF{4} = DFf{4}(1,1,:);

DFf{5} = AFf{2} - DFf{1} - DFf{3} - DFf{4};     DF{5} = DFf{5}(1,:,:);
DFf{6} = AFf{3} - DFf{1} - DFf{2} - DFf{4};     DF{6} = DFf{6}(:,1,:);
DFf{7} = AFf{4} - DFf{1} - DFf{2} - DFf{3};     DF{7} = DFf{7}(:,:,1);

DFf{8} = AFf{1} - DFf{1} - DFf{2} - DFf{3} - DFf{4} - DFf{5} - DFf{6} - DFf{7};              
DF{8} = DFf{8};
end


function plotC(C)
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(C)
colorbar
pos=get(gca,'position');
[rows,cols]=size(C);
width= 1.05*(pos(3)-pos(1))/(cols-1);
height = (pos(4)-pos(2))/(rows-1);
% %create textbox annotations
for i=1:cols
    for j=1:rows              
        annotation('textbox',[pos(1)+width*(i-1),pos(2)+height*(rows-j),width,height], ...
            'string',num2str(C(i,j)),'LineStyle','none','HorizontalAlignment','center',...
            'VerticalAlignment','middle');
    end
end
xticklabels(["$u_0$", "$\bar{u}_t^t$", "$\bar{u}_y^y$", "$\bar{u}_\sigma^\sigma$",...
    "$\widehat{u}_{y\sigma}^{y\sigma}$", "$\underline{u}_{t\sigma}^{t\sigma}$", "$[u]_{ty}^{ty}$", "$u_{t y\sigma}$"]);
yticklabels(["$s_0$", "$\bar{s}_t^t$", "$\bar{s}_y^y$", "$\bar{s}_\sigma^\sigma$",...
    "$\widehat{s}_{y\sigma}^{y\sigma}$", "$\underline{s}_{t\sigma}^{t\sigma}$", "$[s]_{ty}^{ty}$", "$s_{t y\sigma}$"]);

end