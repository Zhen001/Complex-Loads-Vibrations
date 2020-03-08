%% Homework Assignment Complex 2018
% Mike Pesselse
% Sven Klinkhamer

clear all
close all
clc
tic
%clf
%grid on

% Water properties
rho = 1025; %sea water density

% General particulars
L = 85;     %[m] length
Bm = 12;    %[m] beam moulded
D = 8;      %[m] depth to main deck
T = 5;      %[m] draft
Cb = 0.58;  %[-] block coeficient

% Construction material
E = 210e9; % [Pa] steel

% Propulsion train
i_gb = 4.528; % Gearbox reduction ratio

%% Body plan
% Given bodyplan [ord nr,B/Bm,A/Bm*T]
bodyplan = [
    0,.48,.08
    1,.76,.21
    2,.89,.36
    3,.96,.53
    4,1,.65
    5,1,.77
    6,1,.85
    7,1,.91
    8,1,.94
    9,1,.97
    10,1,.93
    11,.98,.89
    12,.94,.82
    13,.88,.72
    14,.80,.62
    15,.68,.50
    16,.56,.39
    17,.39,.28
    18,.23,.16
    19,.09,.12
    20,0,.11];       

%% PART 1 PREPARATION
%% 2)Midship section moment of inertia
h = 0.9;    % [m] tank top height 
t = 0.01;   % [m] plate thickness 
w = 1.5;    % [m] width double wall 

% Determining location neutral axis
%   1                          2 %
% [--]                       [--]%
% [  ]                       [  ]%
% [  ]                       [  ]%
%5[  ]6                     7[  ]8%
% [  ]                       [  ]%
% [  ]                       [  ]%
% [  ]           3           [  ]%
% [-----------------------------]%
%9[        10[      ]11         ]12%
% [-----------------------------]%
%                4               %

y = [w;w;Bm;Bm;t;t;t;t;t;t;t;t];
z = [t;t;t;t;D-h;D-h;D-h;D-h;h;h;h;h];
A = zeros(1,12);    % Area plate
cg = zeros(1,12);   % Center of gravity plate (z=0)

for i=1:12
    if i==1 || i==2
        A(i) = y(i)*z(i);
        cg(i) = D;
    elseif i==3
        A(i) = y(i)*z(i);
        cg(i) = h;
    elseif i==4
        A(i) = y(i)*z(i);
        cg(i) = 0;
    elseif i==5 || i==6 || i==7 || i==8
        A(i) = y(i)*z(i);
        cg(i) = (D-h)/2+h;
    elseif i==9 || i==10 || i==11 || i==12
        A(i) = y(i)*z(i);
        cg(i) = h/2;
    end
end

na = (sum(cg.*A))/(sum(A)); % Calculate z coordinate neutral axis

%Calculation moment of inertia per plate
for i=1:12
    I(i)=(y(i).*z(i).^3)/12+A(i).*(na-cg(i)).^2;
end
%Calculate moment of inertia midship section
I_y = sum(I);

clearvars A h t cg i w y z I na
%% 3)Calculate bending stiffness distribution
ord_x = L/(length(bodyplan)-1);      % Length per ordinate  
EI_dist=I_y*bodyplan(:,2)*E;         % Bending stiffness distribution (EI)

% figure
% plot(bodyplan(:,1)*ord_x,EI_dist)
% title('Bending Stiffness Distribution')
% xlabel('Ship length [m]')
% ylabel('EI [n/m]')
% grid on
% print('plots/3','-dpng')

%% 4)Calculate mass distribution
nabla = L*Bm*T*Cb;                           % Volume ship
m_ship = nabla*rho;                          % Mass ship [kg]

% Mass ship distributed per ordinate [kg/m]
m_dist = (bodyplan(:,3)/sum(bodyplan(:,3)))*m_ship/ord_x;            
        
% figure
% plot(bodyplan(:,1)*ord_x,m_dist)                                     
% title('Mass Distribution')
% xlabel('Ship length [m]')
% ylabel('Mass [kg/m]')
% grid on
% print('plots/4','-dpng')

clearvars nabla

%% 5)Calculate distribution spring stiffness of the water
g = 9.81;       % [m/s^2]

A_w = 786.76;    % [m^2] waterlijn oppervlak uit excel via Simpson (nog berekenen in matlab)

k_water = rho*g*A_w*bodyplan(:,2);                    
% figure
% plot(bodyplan(:,1)*ord_x,k_water)
% title('Water stiffness Distribution')
% xlabel('Ship length [m]')
% ylabel('k [N/m^2]')
% grid on
% print('plots/5','-dpng')

%% PART 2 
%% 6)Calculate distribution of added mass of water
B_ord = bodyplan(:,2).*Bm;  %[m] breedte per ordinaat vector
T_ord = ones(21,1)*T;       %[m] diepgang per ordinaat vector
bdratio = B_ord./T_ord;     %[-] Beam/draft ratio vector


% sectional area coefficient KLOPT NIET
sac = (bodyplan(:,3).*(Bm))./(B_ord);    %[-] 

Cv = [.85;.76;.77;.82;.88;.99;1.09;1.18;1.27;...
    1.36;1.26;1.16;1.125;1.05;.99;.94;0.89;0.91;0.875;2;2]; 
b_ord = B_ord./2; % halve breedte per ordinaat

for n=1:8 %1: heave, 2: pitch, 3-8: bending
    if n==1
        J(1:21,n) = 1; 
    else 
        J(n) = 1.02 - 3*(1.2-(1/n))*(Bm/L);
    end
    
    Mf(:,n) = J(n).*Cv*(pi/2)*rho.*b_ord.^2;%added mass [kg/m]/ordinate/mode   
    m(:,n) = (m_dist+Mf(:,n));        % total mass per ordinate [kg/m]
end

[x,y] = size(Mf);
ord_afstand = 0:1:x-1;

% figure
% for i = 1:y
%     hold on
%     plot(ord_afstand,Mf(:,i));
%     plot(ord_afstand,m_dist,'--');
% end
% 
% title('distribution of added mass of water')
% xlabel('ordinaat [-]')
% ylabel('massa [kg]')
% legend('mode shape 1','mode shape 2','mode shape 3','mode shape 4',...
%'mode shape 5','mode shape 6','mode shape 7','modeshape 8','mass distribution')
% print('plots/6','-dpng')


%% 7) Stiffness and mass matrice of the ship         
A_ord = bodyplan(:,3)*(Bm*T);

M1 = zeros(44);
Mn = zeros(44);
K1 = zeros(44);
Kn = zeros(44);

for i = 1:21
    M_m(:,:,i) = BeamM(rho*A_ord(i),L/20);
    K_m(:,:,i) = BeamK(EI_dist(i),L/20);   
        
    n = 2*i-1;
    K1(n:n+3,n:n+3,i) =  K_m(:,:,i);
    Kn(:,:,i+1) = Kn(:,:,i) + K1(:,:,i);
        
    M1(n:n+3,n:n+3,i) = M_m(:,:,i);
    Mn(:,:,i+1) = Mn(:,:,i) + M1(:,:,i); 
end

K = Kn(1:42,1:42,end);
M = Mn(1:42,1:42,end);

% clearvars Kn Mn n M_m K_m M1 Mn K1 Kn
%% 8) Solve for the first 8 natural frequencies and mode shapes
%K distribution water (k_water) + K matrix schip (K)
Kd = zeros(42); %K distribution water optellen bij stiffness matrice ship
for i=1:21
    n=2*i-1;
    Kd(n,n) = k_water(i);
end
K = Kd + K;

%Massa matrixen ook nog bij elkaar optellen M en Mf
%hiero

%Berekenen van eigenfrequenties
omega = abs(sqrt(eig(inv(M)*K)));
for i = 1:8
    ev(:,i) = (1/(2*pi))*null(K-M*omega(i)^2);
end

% %Lengte schip opdelen in 42 stukken, halve ordinaat 
% L2 = ones(1,42); 
% for i = 1:42
%    L2(1,i) = (ord_x/2)*i;
% end
% 
% figure
% plot(L2,ev)
% title('First 8 mode shapes')
% xlabel('lengte [m]')
% ylabel('w [Hz]')
% legend('mode shape 1','mode shape 2','mode shape 3','mode shape 4',...
% 'mode shape 5','mode shape 6','mode shape 7','modeshape 8')
% %print('plots/8','-dpng')

ev2 = zeros(0.5*42,8);

for n = 1:8
    for i = 1:20
        ev2(i,n) = (1/(2*pi))*ev(2*i-1,n);
    end
end

figure
plot(bodyplan(:,1)*ord_x,ev2)
title('First 8 mode shapes')
xlabel('lengte [m]')
ylabel('w [Hz]')
legend('mode shape 1','mode shape 2','mode shape 3','mode shape 4',...
'mode shape 5','mode shape 6','mode shape 7','modeshape 8')
%print('plots/8','-dpng')

toc

%PROBEREN: INPLAATS VAN 21 MATRIX -> 42 TE MAKEN, 42 MATRIX NAAR 21 -> MBV
%ODD FUNCTION
