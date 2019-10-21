%% MAE 480 Project
% Sarah Nguyen
% Michael Angeles

%% Problem 1: Pitching moment coefficient for the fuselage
clear,clc
% Fuselage Dimensions
Sections = [1:8]';
Dx = [10;7;2;12;8.25;8.75;12;3.75];

cy_2 = (7/3)*(9+(2*4.5))/(9+4.5)
bf_2 = 9 - (2*((9+2*4.5)*(9-4.5))/(6*(9+4.5)))

cy_6 = (8.75/3)*(9+(2*7.214))/(9+7.214)
bf_6 = 9 - (2*((9+2*7.214)*(9-7.214))/(6*(9+7.214)))

cy_7 = (12/3)*(7.214+(2*4.765))/(7.214+4.765)
bf_7 = 7.214 - (2*((7.214+2*4.765)*(7.214-4.765))/(6*(7.214+4.765)))

cy_8 = (3.75/3)*(4.765+(2*4))/(4.765+4)
bf_8 = 4.765 - (2*((4.765+2*4)*(4.765-4))/(6*(4.765+4)))

% bf of trapezoid = d-2((d+2a)(d-a)/6(d+a))
bf = [4.5;5.841;9;9;9;8.14;6.07;4.39];
x1 = [14;5.11;2;0;4.125;12.464;22.591;30.82];
x1_cre = [1.1667;0.426;0.1667;0;0;0;0;0];

% Find Upwash Values (figure9b)
deda = @(x1_cre) -0.0442*x1_cre.^3 + 0.2946*x1_cre.^2 - 0.6771*x1_cre + 0.6100;
deuda1 = deda(1.166667);
deuda2 = deda(0.425926);

% Find Upwash near LE (figure9a)
deda1 = @(x1_cre) -6.4815*x1_cre.^3 + 15.5556*x1_cre.^2 - 13.0833*x1_cre ...
           + 4.7593;
deuda3 = deda1(0.166667);

% Tail Dimensions
cr = 12;
ct = 4;
lambdaT = ct/cr;
c_bar_t_integral = @(y) (0.263374*y^3-12.642*y^2+202.72*y);
c__bar_t_lower = c_bar_t_integral(2.5)*(2/204);
c_bar_t_upper = c_bar_t_integral(11.5)*(2/204);
c_bar_t_sum = c_bar_t_upper - c__bar_t_lower+(144*2.5)*(2/204);
cbar_t = c_bar_t_sum;
unknownd= 12-cbar_t;

% Find lh
lh = (cbar_t/4)+8.25+8.75+unknownd;
x1_lh = x1(5:8)/lh;
X1_lh = [0;0;0;0;x1_lh];

% Find Downwash (eq 3.10)
Deda = 0.430431;
deuda5 = (x1_lh(1)*(1-Deda))-1;
deuda6 = (x1_lh(2)*(1-Deda))-1;
deuda7 = (x1_lh(3)*(1-Deda))-1;
deuda8 = (x1_lh(4)*(1-Deda))-1;

deuda = [deuda1;deuda2;deuda3;0;deuda5;deuda6;deuda7;deuda8];

% Wing
crw = 12;
ctw = 6.61;
b = 62;
lambdaw = ctw/crw;
cbarw = (2/3)*crw*(1+lambdaw+lambdaw^2)/(1+lambdaw);
S = 605.288;

% Find Cmaf
for i = bf(1:8)
    i=bf.^2;
    for j = deuda(1:8)
        j = 1+j;
    end
    for k = Dx(1:8)
        k = k;
    end
    l = i.*j.*k;
    m = sum(l)-l(4)
    Cmaf = (pi()/(2*S*cbarw))*m;
end


% Table of upwash and downwash values
T = table(Sections,Dx,bf,x1,x1_cre,deuda,X1_lh)
Cmaf

%% Problem 2: Lift Curve Slope of the wing
% Givens
V = 60;                   % Steady-State Speed [ft/s]
h = 500;                  % Cruise Altitude [ft]
Gamma = 1.4;              % Air Density  
R = 287;                  % Ideal Gas Constant [J/kg-K]
T = 287.1594;             % Temperature [K] - found from appendix at h = 500 ft
phiTE = 15;               % [degree]
W.T_C = 0.124;            % Wing ratio
Re = 10^6;
Sweep_Angle_W = 11.5;     % [degree]
a = sqrt(Gamma*R*T)*3.28084;      % Speed of Sound [ft/s]
M = V/a;                  % Mach Number
Beta = sqrt(1-M^2);


% Function to solve for aw
ao = @(ao_theory,ao_O_ao_theory,M) (1.05/sqrt(1-M^2))*(ao_O_ao_theory)*(ao_theory);
k = @(ao) ao/(2*pi);
alpha_w = @(A,mid_slope,k) (2*pi*A)/(2+sqrt((A^2*Beta^2/k^2)*(1+(tan(mid_slope)^2/Beta^2))+4));

%Solutions using Given Function 3-13a
ao_theory_Main = fig3_13a(W.T_C);

%Solutions using Given Function 3-13b
ao_O_ao_theory_Main = fig3_13b(deg2rad(tan(phiTE/2)),Re)

main_midcs = atan(tand(Sweep_Angle_W) - ((12-6.61)/62))

ao_Main = ao(ao_theory_Main,ao_O_ao_theory_Main,M);

k_Main = k(ao_Main);
Area_Main = 605.288;    %in2
Area_Horz = 204;       %in2

b_Main = 62;
b_Horz = 23;

AR = @(b,S) b^2/S;
AR_Main = AR(b_Main, Area_Main)
AR_Horz = AR(b_Horz, Area_Horz);

aW = @(AR,k,midcs) ((2*pi)*AR)/(2+sqrt(((AR^2)*(Beta^2)/(k^2))*(1+(((tan(midcs))^2)/Beta))+4))
aW_Main = aW(AR_Main,k_Main,main_midcs)
%% Problem 3: Lift Curve Slope of the Horizontal Tail
% Givens
V = 60;                   % Steady-State Speed [ft/s]
h = 500;                  % Cruise Altitude [ft]
Gamma = 1.4;              % Air Density  
R = 287;                  % Ideal Gas Constant [J/kg-K]
T = 287.1594;             % Temperature [K] - found from appendix at h = 500 ft
H.T_C = 0.080;            % Tail ratio
Re = 10^6;
Cr = 12;
Ct = 4;
b = 23;
lambda = Ct/Cr;
S = 204;
AR = b^2/S
% NACA 00-0XX8
a = sqrt(Gamma*R*T)*3.28084;      % Speed of Sound [ft/s]
M = V/a;                  % Mach Number
Beta = sqrt(1-M^2);

% Function to solve for aw
ao = @(ao_theory,ao_O_ao_theory,M) (1.05/sqrt(1-M^2))*(ao_O_ao_theory)*(ao_theory);
k = @(ao) ao/(2*pi);
at = @(AR,k,midcs) ((2*pi)*AR)/(2+sqrt(((AR^2)*(Beta^2)/(k^2))*(1+(((tan(midcs))^2)/Beta))+4))

%Solutions using Given Function 3-13a
ao_theory_Horz = fig3_13a(H.T_C); 
% Solutions using Given Function 3-13c
phiTE = fig3_13c(H.T_C)

%Solutions using Given Function 3-13b
ao_O_ao_theory_Horiz = fig3_13b((tand(phiTE/2)),Re)

sweep_angle = 90 - atand(9/8)
tail_midcs = atan(tand(sweep_angle) - ((12-4)/23))

ao_Horz = ao(ao_theory_Horz,ao_O_ao_theory_Horiz,M)
k_Horz = k(ao_Horz)
aW_tail = at(AR,k_Horz,tail_midcs)

%% Problem 4
rho_ssl = 0.07651; %lb/ft^3

de_o_da_tail = 2*aW_tail/(pi*AR_Horz);

%C_fw = fig3_24(c_bar_t_sum,M,1000000);
C_fw = 0.045;
L = 2;
Main_Lambda = 0.124;

R_LS = @(Lam_LE_w) -1.0400*cos(Lam_LE_w)^2 + 2.0600*cos(Lam_LE_w) + 0.0400;     % Fig 3.25
R_LS_Wing = R_LS((deg2rad(lambda)));

Swet_o_S = (2*Area_Main)/Area_Main;

C_Do_w = C_fw*(1+L*0.124+100*0.124^4)*R_LS_Wing*Swet_o_S;

dq_o_q = 2.42*(sqrt(C_Do_w))/((lh/cbarw)+0.3);

neta = 1-dq_o_q;

lift_curve_slope = aW_Main+aW_tail*(1-deg2rad(1)-de_o_da_tail)*neta*(Area_Horz/Area_Main);

%% Problem 5
lt = 18.254;
lt = 22.004
bar_v1 = (Area_Horz*lt)/(Area_Main*cbarw);

x_cg_bar = (13.06)/cbarw;
x_cg_bar = (32.057-30)/cbarw
x_bar_ac_wb = 0.25;
x_bar_a = x_cg_bar-x_bar_ac_wb;


dcm_o_da = aW_Main*x_bar_a+Cmaf-aW_tail*(1-de_o_da_tail)*neta*bar_v1;

%% Problem 6
N_O =  x_bar_ac_wb - (Cmaf/lift_curve_slope) + (aW_tail/aW_Main)*(1-de_o_da_tail)*bar_v1*neta

H = N_O - x_cg_bar
