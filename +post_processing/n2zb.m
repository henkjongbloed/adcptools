function zb = n2zb(nvec, mesh)
% Get bed elevation from n coordinate, mesh m.

nvec(nvec>max(mesh.nb_all)) = max(mesh.nb_all);
nvec(nvec<min(mesh.nb_all)) = min(mesh.nb_all);


zb = interp1(mesh.nb_all, mesh.zb_all, nvec, 'makima');



end