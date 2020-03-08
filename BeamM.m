% This function returns the mass matrix for a 2D, two-node beam element.
% The mass matrix relates the external forces and moments on the beam
% element to the translations and rotations of the nodes. M is a 4x4 matrix.
% [M]*[w1" phi1" w2" phi2"]' =[F1 M1 F2 M2]'
% In which w1", w2" are the transverse accelerations and phi1", phi2" are
% the rotational accelerations of the nodes 1 and 2. F1, F2 are the external
% transverse forces applied to the beam element and M1, M2 are the external
% moments applied to the beam element. 
function M = BeamM(rhoA, L)
M = [156 22*L 54 -13*L;
    22*L 4*L^2 13*L -3*L^2;
    54 13*L 156 -22*L;
    -13*L -3*L^2 -22*L 4*L^2] .* rhoA.*L/420;