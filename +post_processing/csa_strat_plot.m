function csa_strat_plot(didx, bidx, F, D, AF)
%% CSA fig



llimpos = [.3, .4];
llimneg = [-.4, -.3];


vlimbot = [.1, .2];
vlimsur = [.8, .9];

% close all
% Instantaneous Function Bottom / Surface
for di = didx
    for bi = bidx
        for vi = 1:3 %Flow, Salt, Flux
            for ri = 1
                for c = 1:size(F{di}{bi}{vi}, 2)

                    % LATERAL SHEAR VARIABLES
                    IFP{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_llim(F{di}{bi}{vi}{ri, c}, llimpos);
                    IFN{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_llim(F{di}{bi}{vi}{ri, c}, llimneg);

                    daIFP{di}{bi}{vi}{ri, c} = D{di}{bi}.da(IFP{di}{bi}{vi}{ri, c});
                    daIFN{di}{bi}{vi}{ri, c} = D{di}{bi}.da(IFN{di}{bi}{vi}{ri, c});


                    % BOTTOM AND SURFACE VARIABLES
                    IFB{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(F{di}{bi}{vi}{ri, c}, vlimbot); % U-Flow, Salinity, Flux instantaneous
                    IFS{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(F{di}{bi}{vi}{ri, c}, vlimsur); % U-Flow, Salinity, Flux instantaneous

                    % Laterally averaged bottom and surface variables.
                    laIFB{di}{bi}{vi}{ri, c} = D{di}{bi}.la(IFB{di}{bi}{vi}{ri, c});
                    % laterally averaged stratification for all components, function of (t)
                    laIFS{di}{bi}{vi}{ri, c} = D{di}{bi}.la(IFS{di}{bi}{vi}{ri, c});

                    % Tidally averaged bottom and surface variables.
                    twaIFB{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(AF{di}{bi}{vi}{ri, c}{5}, vlimbot); % TWA of the above
                    twaIFS{di}{bi}{vi}{ri, c} = D{di}{bi}.extract_inst_vlim(AF{di}{bi}{vi}{ri, c}{5}, vlimsur); % TWA of the above
                end
            end
        end
    end
end


for di = didx
    [fig, TL] = prep_fig([8,12], [4,1]);
    %     figure;
    csi_figs(di, bidx, D, AF) % (didx, bidx, D, flow, AF, tak_name)
    fontsize(fig, 10, "points")
    %     [fig, TL] = prep_fig([9,12], [1,1]);
    %      get_label(di, bname);
    [fig, TL] = prep_fig([8, 9], [3, 2]);
    
    %     figure
    strat_figs(di, bidx, D, laIFB, laIFS, daIFP, daIFN)
    fontsize(fig, 10, "points")
end
end