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
%lh = 13.08
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
HT_C = 0.080;            % Tail ratio
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
ao_theory_Horz = fig3_13a(HT_C); 
% Solutions using Given Function 3-13c
phiTE = fig3_13c(HT_C,'00XX-X8')

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
bar_v1 = (Area_Horz*lt)/(Area_Main*cbarw);

x_cg_bar = (13.06)/cbarw

x_bar_ac_wb = 0.25;
x_bar_a = x_cg_bar-x_bar_ac_wb;


dcm_o_da = aW_Main*x_bar_a+Cmaf-aW_tail*(1-de_o_da_tail)*neta*bar_v1;

%% Problem 6
N_O =  x_bar_ac_wb - (Cmaf/lift_curve_slope) + (aW_tail/aW_Main)*(1-de_o_da_tail)*bar_v1*neta;

H = N_O - x_cg_bar;

%% Report II
% % Move the wing and tail backwards
% lt = 22.004;                 % Distance from the Cg to the AC of the tail
% bar_v1 = (Area_Horz*lt)/(Area_Main*cbarw);
% x_cg_bar = (32.057-30)/cbarw;
% x_bar_ac_wb = 0.25;
% x_bar_a = -(x_cg_bar-x_bar_ac_wb);
% 
% 
% dcm_o_da = aW_Main*x_bar_a+Cmaf-aW_tail*(1-de_o_da_tail)*neta*bar_v1
% N_O =  x_bar_ac_wb - (Cmaf/lift_curve_slope) + (aW_tail/aW_Main)*(1-de_o_da_tail)*bar_v1*neta;
% 
% H = N_O - x_cg_bar

% Move CG location
N_O = 0.2265;       % From report I
Xcg_wrtF = 21.104;
Xcg_bar = (Xcg_wrtF-19)/cbarw;
H = N_O - Xcg_bar;
Xac_w = 0.25;
x_bar_a = (Xcg_bar - Xac_w)

lt = 30.226;                 % Distance from the Cg to the AC of the tail
bar_v1 = (Area_Horz*lt)/(Area_Main*cbarw)

%% Problem 1

cf = 1.25;      % Chord Flap
ct = cbar_t;    % Tail chord
cfc = cf/ct;    % Flap Chord Factor
HT_tc = 0.08;   % Horizontal tail ratio (t/c)

adCLadcl = fig3_35(cfc,AR_Main);       % Use AR of the wing
yo = (b_Horz/2)-1.5;            % outboard location (fig. 3.34)
yi = (b_Horz/2)-8.5;            % inboard location (fig. 3.34)
neta_o = 2*yo/b_Horz;           % Based on control surface
neta_i = 2*yi/b_Horz;
K_bo = fig3_36(neta_o,lambdaT);
K_bi = fig3_36(neta_i,lambdaT);
K_b = K_bo-K_bi;

cld_theory = fig3_37_a(cfc,HT_tc);

cld_cld_theory = fig3_37_b(cfc,ao_O_ao_theory_Horiz);
cld = cld_cld_theory*cld_theory;
Tau = (cld/ao_Horz)*adCLadcl*K_b;
cmd_e = -aW_tail*bar_v1*neta*Tau     % [/rad]

%% Problem 2

W = 12;     % Weight [lbf]
rho_500 = ((1.2017-1.225)/200)*(152.4) + 1.225;      % Linear Interpolation (Back of the text)
rho_500 = rho_500/515.379;      % Density [slug/ft^3]
v = 60;         % Velocity at steady flight [ft/s^2]
Sw = 4.1825;       % Area of the wing [ft^2]
Sw_i = 602.285;            % Area of the wing [in^2]
Cl = (2*W)/(rho_500*(v^2)*Sw)

Clmax = 1.1;        % Given
dcmdcl = -H;
detrim = -0.05;     % Given [rad]
detrim0 = detrim + (dcmdcl/cmd_e)*Cl
demax = -25;        % Negative so that it produces an upward deflection (p.226)

xcgf = N_O - (deg2rad(demax)-detrim0)*(cmd_e/Clmax)

dxcgf = N_O - xcgf
%% Problem 3
VTH = 10;           % Given [in]
HVT = 2;            % Given [in]
bv = VTH+HVT;       % Given [in]
Cvt = Ct;           % Given [in]
Cvr = 19;           % Given [in]
lambdaV = Cvt/Cvr;
r1 = 3.2808;        % Given [in]
Sv = (bv/2)*Cvr*(1+lambdaV);        % Vertical Tail Area [in^2]
Sw_i = 602.285;                     % Wing Reference Area [in^2]
Av = 2*bv/(Cvr*(1+lambdaV));        % Aspect Ratio
cbarV = (2/3)*Cvr*((1+lambdaV+lambdaV^2)/1+lambdaV);
ymac = 2*(bv/6)*(1+2*lambdaV)/(1+lambdaV);
kv = fig3_75(bv/(2*r1));
VTtc = 0.075;
ao_theory_VT = fig3_13a(VTtc);       % Given t/c of VT
phiTE_VT = fig3_13c(VTtc,'00XX-X8');
Re = 10^6;          % Given
%Solutions using Given Function 3-13b
ao_O_ao_theory_VT = fig3_13b((tand(phiTE_VT/2)),Re);

tail_LE_VT = 90 - atand(10/15);     % Found from diagram on paper
VT_c_4_angle = atand(tand(tail_LE_VT) - (Cvr - Cvt)/(2*10));    % Quarter-chord angle (b = 10)
VT_c_2_angle = atand(tand(tail_LE_VT) - (Cvr - Cvt)/(10));    % Half-chord angle (b = 10)
ao_VT = ao(ao_theory_VT,ao_O_ao_theory_VT,M);
k_VT = ao_VT/(2*pi);
AVB_AV = fig3_77(bv/(2*r1),lambdaV);
Av_eff = AVB_AV*Av;
av_VT = ((2*pi)*Av_eff)/(2+sqrt(((Av_eff^2)*(Beta^2)/(kv^2))*(1+(((tan(VT_c_4_angle))^2)/Beta))+4))
zw = 2;
dfmax = 4.166;      % Found in CAD
sidewashdynamicratio = 0.724 + (3.06*(Sv/Sw_i)/(1+cosd(VT_c_4_angle))) + (0.4*zw/dfmax) + 0.0009*Av

Cyb = -kv*av_VT*sidewashdynamicratio*(Sv/Sw_i)

quarterchordpt = cbarV/4;
% Directional Stability caused by wing sweep
% cnb = (cnb)_w + (cnb)_bw + (cnb)_vfix
% (cnb)_w = (cnb)_circular + (cnb)_vw

cnb_vw = (Cl^2)*((1/(4*pi*AR_Main)) - (tand(VT_c_4_angle)/(pi*AR_Main*(AR_Main + 4*cosd(VT_c_4_angle))))*(cosd(VT_c_4_angle) -(AR_Main/2) - (AR_Main^2)/(8*cosd(VT_c_4_angle)) + 6*x_bar_a*sind(VT_c_4_angle)/AR_Main))
cnb_w = cnb_vw
lf = 63.75;     % Length of the fuselage
Sbs = 343;         % Found by Michael on paper
lf_Sbs = (lf^2)/Sbs
h1 = 8.317;     % Found from CAD
h2 = 8.067;     % Found from CAD
sh1_h2 = sqrt(h1/h2)
h_bfmax = 1     % Found from CAD (Total height of aircraft/ max width) = (9/9)
xm = 21.104;    % Distance from nose to CG
xm_lf = xm/lf
Kn = 0.0008;    % Found from figure 3.73
Krl = 1;        % Found from figure 3.74
cnb_bw = -Kn*Krl*(Sbs/Sw_i)*(lf/10)         % b = 10
k_ = fig3_75(bv/(2*r1));        % Found from figure 3.75
Vbar_2 = (Sv*lt)/(Sw_i*b_Main);  
cnb_vfix = k_*av_VT*sidewashdynamicratio*Vbar_2;

Cnb = cnb_w + cnb_bw + cnb_vfix

%% Problem 4
neta_t = neta; % From problem 5 on report 1
neta_v = neta_t;    % Assumed based on book and Kanistras 
cf_r = 1;     % Possible to be changed *** - The width of the rudder
cvt = cbarV;   % Pulled from report 1 - cbar of the vertical tail
cfc_r = cf_r/cvt; 

% Finding tau_r = tau2
adCLadcl_r = fig3_35(cfc_r,AR_Main);
yi_r = 1.5+2;  % Using fig 3.34 and adding the 2 in to the center of the fuselage
yo_r = 7+1.5+2; % Using fig 3.34 and adding the 2 in to the center of the fuselage
neta_i_r = (yi_r)/bv; 
neta_o_r = (yo_r)/bv;
K_bi_r = fig3_36(neta_i_r,lambdaV);
K_bo_r = fig3_36(neta_o_r,lambdaV);
K_b_r = K_bo_r - K_bi_r;
cld_theory_r = fig3_37_a(cfc_r,0.075); %0.075 is the t/c of the vertical tail
cld_cld_theory_r = fig3_37_b(cfc_r,ao_O_ao_theory_VT);
cld_r = cld_cld_theory_r*cld_theory_r;
Tau_r = (cld_r/ao_VT)*adCLadcl_r*K_b_r;

% Calculate Cndr with values from Problem 3
cnd_r = -k_*neta_v*Vbar_2*av_VT*Tau_r

%% Question 5
% Find lateral stability coefficient for the entire aircraft
b_Tail = 10;
Area_Tail = 52.5;
AR = @(b,S) b^2/S;
AR_Tail = AR(b_Tail, Area_Tail);

r_gam = 0; % No dihedral angle

CLB_CL_c_2 = fig3_96(VT_c_2_angle,AR_Tail);

K_MLambda_x = M*cosd(VT_c_2_angle);
factor = AR_Tail/cosd(VT_c_2_angle);
K_MLambda = fig3_97(K_MLambda_x,factor) %Code says A_cosc2, but figure uses Aspect Ratio

lf_prime_b = 56.25/b_Tail; %Measurement to half chord point on vertical tail to tip of the aircraft (NEEDS TO MOVE TO Orginal Location)
K_f = fig3_98(lf_prime_b,AR_Tail);

CL_B_CL_A = fig3_99(AR_Tail,lambdaV); % Per degree? all the others weren't

CL_B_Gamma = fig3_100(AR_Tail, VT_c_2_angle);

K_MGamma = fig3_101(K_MLambda_x,AR_Tail);

C_lp = sqrt(2)/2; % Assumed

d= 9; % Average fuselage diameter
DCLB_Gamma = -0.0005*sqrt(AR_Tail)*(d/b)^2;
DCLB_zw = ((1.2*sqrt(AR_Tail))/57.3)*(zw/b)*((2*d)/b);

CL_B_WB = Cl*(CLB_CL_c_2*K_MLambda*K_f+CL_B_CL_A)+(r_gam*(CL_B_Gamma*K_MGamma+DCLB_Gamma))+DCLB_zw;
k_CL_B_VT = 0.075;

zv = zw; % vertical distance between the cg and the vertical tail AC
lv = lt; % between the cg and the vertical tail AC
alpha = 0;
CL_B_VT = -k_CL_B_VT*av_VT*(sidewashdynamicratio)*(Area_Tail/Area_Main)*((zv*cos(alpha)-lv*sin(alpha))/b);

CLB = CL_B_WB + CL_B_VT
%% Problem 6 

% Finding Tau_a
cf_w = 1.50;        % Chord Flap
cw = cbarw;       % Main wing chord
cfc_w = cf_w/cw;    % Flap Chord Factor
tc_w = 0.124;       % Given

adCLadcl_w = fig3_35(cfc_w,AR_Main);       % Use AR of the wing
yo_w = 62/2;            % outboard location (fig. 3.105 pg 328)
yi_w = 62/2 - 24.5;     % inboard location (fig. 3.105 pg 328)
neta_o_w = 2*yo_w/b_Main;           % Based on control surface
neta_i_w = 2*yi_w/b_Main;
K_bo_w = fig3_36(neta_o_w,lambdaw);
K_bi_w = fig3_36(neta_i_w,lambdaw);
K_b_w = K_bo_w - K_bi_w;

cld_theory_w = fig3_37_a(cfc_w,tc_w);

cld_cld_theory_w = fig3_37_b(cfc_w,ao_O_ao_theory_Main);
cld_w = cld_cld_theory_w*cld_theory_w;
Tau_a = (cld_w/ao_Main)*adCLadcl_w*K_b_w;

% Solving for the Integral 
c_yy = @(y) (y^2)/2 + (2/3)*((lambdaw-1)/b_Main)*y^3;
clda = -(2*aW_Main*Tau_a*crw)/(S*b_Main)*(c_yy(yo_w)-c_yy(yi_w))
