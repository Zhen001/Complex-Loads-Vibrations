% This function returns the stiffness matrix for a 2D, two-node beam element.
% Shear is not taken into account, only bending.
% The stiffness matrix relates the external forces and moments on the beam
% element to the translations and rotations of the nodes. K is a 4x4 matrix.
% [K]*[w1 phi1 w2 phi2]' =[F1 M1 F2 M2]'
% In which w1, w2 are the transverse translations and phi1, phi2 are the
% rotations of the nodes 1 and 2. F1, F2 are the external transverse
% forces applied to the beam element and M1, M2 are the external moments
% applied to the beam element. 
function K = BeamK(EI,L)
K = [12 6*L -12 6*L;
    6*L 4*L^2 -6*L 2*L^2;
    -12 -6*L 12 -6*L;
    6*L 2*L^2 -6*L 4*L^2] * EI/L^3;