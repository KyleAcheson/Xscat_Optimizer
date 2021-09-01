function [factor] = f_functions_electron()
% [f_H, f_C, f_Xe] = f_functions()
% Generates interpolating functions for H, C, & Xe scattering factors. From
% International Tables for Crystallography, Volume C, section 6.1 (p. 476)
%
% Using interpolating function (6.1.1.15), p. 487:
% f(q) := c + Sum(a_i*Exp(-b_i (q/(4*pi))^2),i,1,4)
% constants are given in table 6.1.1.4
% functions valid for q < 8*pi
% error magnitudes also given in table 6.1.1.4

% Note that q is given in inverse Angstroms
% Also the original text gives the equation as a function of
% sin(theta)/lambda, but James has converted this using 
%    q = 2k sin(theta) = 4pi sin(theta)/lambda

% general function:
% a:   vector
% b:   vector
% c:   scalar
% Z_i: charge for each nucleus
% q:   scattering vector; q:= 4*pi/lambda * Sin(theta)
f = @(a,b,c,Z_i,q) (Z_i - (c + sum(a .* exp(-b *(q/(4*pi))^2))))/q^2;

% Hydrogen (SDS)
% max error: 0.000; mean error: 0.000
a_H = [ 0.493002, 0.322912, 0.140191,  0.040810];
b_H = [10.5109,  26.1257,   3.14236,  57.7997];
c_H = 0.003038;
Z_i_H = 1;

% Carbon (HF)
% max error: 0.001; mean error: 0.000
a_C = [2.26069, 1.56165,  1.05075, 0.839259];
b_C = [22.6907, 0.656665, 9.75618, 55.5949];
c_C = 0.286977;
Z_i_C = 6;

% Nitrogen (RHF)
% max error: 0.007; mean error: 0.002
a_N = [12.2126, 3.13220,  2.01250, 1.16630];
b_N = [.005700, 9.89330, 28.9975, .582600];
c_N = -11.529;
Z_i_N = 7;

% Fluorine (RHF)
% max error: .001; mean error: .000
a_F = [3.53920, 2.64120, 1.51700, 1.02430];
b_F = [10.2825, 4.29440, .261500, 26.1476];
c_F = 0.277600;
Z_i_F = 9;

% Sulfur (RHF)
% max error: .005; mean error: .002
a_S = [6.90530, 5.20340, 1.43790, 1.58630];
b_S = [1.46790, 22.2151, 0.253600, 56.1720];
c_S = 0.866900;
Z_i_S = 16;

% Iodine (RHF)
% max error: .037; mean error: 0.009
a_I = [20.1472, 18.9949, 7.51380, 2.27350];
b_I = [4.34700, .381400, 27.7660, 66.8776];
c_I = 4.07120;
Z_i_I = 53;

% Xenon (RHF)
% max error: 0.038; mean error: 0.009
a_Xe = [20.2933, 19.0298, 8.97670, 1.99000];
b_Xe = [3.92820, 0.344000, 26.4659, 64.2658];
c_Xe = 3.71180;
Z_i_Xe = 54;


factor.H  = @(q) arrayfun(@(q) f(a_H, b_H, c_H, Z_i_H, q),q);
factor.C  = @(q) arrayfun(@(q) f(a_C, b_C, c_C, Z_i_C, q),q);
factor.N  = @(q) arrayfun(@(q) f(a_N, b_N, c_N, Z_i_N, q),q);
factor.F  = @(q) arrayfun(@(q) f(a_F, b_F, c_F, Z_i_F, q),q);
factor.S  = @(q) arrayfun(@(q) f(a_S, b_S, c_S, Z_i_S, q),q);
factor.I  = @(q) arrayfun(@(q) f(a_I, b_I, c_I, Z_i_I, q),q);
factor.Xe = @(q) arrayfun(@(q) f(a_Xe,b_Xe,c_Xe,Z_i_Xe,q),q);

end
