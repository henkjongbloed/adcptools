function make_legends_colorbars()

maxs{1} = 18;
maxs{2} = 3; % max depth-averaged salinities -> integers
maxs{3} = 23;
maxs{4} = 6; % max depth-averaged salinities -> integers


% Four colorbars for the map plots. Residual N, S, Instantaneous N, S
fig = makefigure(10, 5);
for i = 1:4
ax=subplot(2,2,i);
    c = colorbar(ax);
    ax.Visible = 'off';
    colormap(helpers.cmaps('salmap'))
    clim([0, maxs{i}])
    c.Ticks = linspace(0, maxs{i}, 3);
    %c.LineWidth = 3
end


fontsize(fig, 10, "points")
fontname(fig, "Book Antiqua")

%%
fig = makefigure(10, 5);

leg = {["HC", "OMS", "OMN"], ["NM", "OM", "NWW"]};
for di = 1:2
    ax=subplot(1,2,di);
    %c = colorbar(ax);
    
    bar([1;2;3], eye(3,3))
    
    legend(leg{di})
    %ax.Visible = 'off';
    %c.LineWidth = 3
end
fontname(fig, "Book Antiqua")
fontsize(fig, 10, "points")

%%
fig = makefigure(10, 5);

legd = ["Waterlevel"];
plot([1,2,3], [0,1,0], 'k')
legend(legd)

fontname(fig, "Book Antiqua")
fontsize(fig, 10, "points")