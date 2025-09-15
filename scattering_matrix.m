clc, clear all, close all

% Scattering matrix

H = 1; C0 = 1; Bi = 1; Bj = 0.1; Bk = 0.9;

a = 1;

incoming = [a, 0, 0]';

A = [Bi*C0, Bj*C0, Bk*C0;
    -H, H, 0;
    0, -H, H];

B = [Bi*C0, Bj*C0, Bk*C0;
    H, -H, 0;
    0, H, -H];


outgoing = A^(-1)*B*incoming