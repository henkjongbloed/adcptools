function Tlim = get_Tlim(flow, bidx)

for bi = bidx
    tim{bi} = flow{1, bi}.solver.adcp.time;
end
Tlim(1)= min([tim{:}]);

M2T = flow{bidx(1)}.solver.model.periods(1,1)/(3600*24);
Tlim(2) = Tlim(1) + M2T;

end