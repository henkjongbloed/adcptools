function csa_lab = get_label(di, bname)
% Make a dummy fig and generate legend to be used in csa figs
hold on
plot(1:3, 1:3, 'k')
c=get_c();
% colororder('gem')
b=bar(rand(3,3));
for idx = 1:numel(b)
    b(idx).FaceColor = c{idx};
end
legend(["Waterlevel", bname{di}{:}],'AutoUpdate','off')
fontname(gcf, "Book Antiqua")


end