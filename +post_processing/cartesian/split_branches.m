function ctd_new = split_branches(ctd, DS, BN)

if strcmp(DS, 'OMHA14')
    ind = {1:929, 930:2639, 2640:length(ctd.T)}; %from inspection
    for b = 1:3
        ctd_new{b}.T =ctd.T(ind{b});
        ctd_new{b}.X =ctd.X(ind{b});
        ctd_new{b}.Y =ctd.Y(ind{b});
        ctd_new{b}.Z =ctd.Z(ind{b});
        ctd_new{b}.S =ctd.S(ind{b});

        %ctd_new{b}.C =ctd.C(ind{b});
        %ctd_new{b}.Temp =ctd.Temp(ind{b});
        ctd_new{b}.Date = ctd.Day(ind{b});

        %ctd_new{b}.Omg =ctd.Omg(ind{b}); %specific to this dataset - not
        %really needed.
        ctd_new{b}.BN = BN{b};
    end
elseif strcmp(DS, 'NMOMNW15')
    ind = {1:3046, 3047:5724, 5725:length(ctd.T)}; %from inspection
    for b = 1:3
        ctd_new{b}.T =ctd.Time(ind{b});
        ctd_new{b}.X =ctd.X(ind{b});
        ctd_new{b}.Y =ctd.Y(ind{b});
        ctd_new{b}.Z =ctd.Z(ind{b});
        ctd_new{b}.S =ctd.S(ind{b});
        
        %ctd_new{b}.C =ctd.C(ind{b});
        %ctd_new{b}.Temp =ctd.T(ind{b});
        

        ctd_new{b}.Date =ctd.Date(ind{b});
        ctd_new{b}.T2 = ctd.Time2(ind{b});
        %ctd_new{b}.Omg =ctd.Omg(ind{b});
        ctd_new{b}.BN = BN{b};
    end
end



ctd_new = ctd_new';

end