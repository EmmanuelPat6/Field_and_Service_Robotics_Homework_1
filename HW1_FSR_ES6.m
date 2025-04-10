clc
clear
close all
syms I m l d theta psi

q = [theta, psi, I];
A = [I+m*(l+d)^2, m*(l+d)^2, 0];

G=null(A);
disp('Vector Fields of the Kinematic Model:')
G
if length(G)>1
    Lie_Bracket = jacobian(G(:,2),q)*G(:,1) - jacobian(G(:,1),q) * G(:,2);
    disp('Lie Bracket [g1,g2]: ')
    Lie_Bracket
end

F = [G(:,1), G(:,2), Lie_Bracket];
disp('F Matrix related to the Accessibility Distribution Δ_A : ')
F
disp('rank(F): ')
disp(rank(F))
if(rank(F)==3)
    disp('The System is Controllable and Completely Nonholonomic')
else
    disp('The System is Not Controllable, Holonomic and Integrable')
end
