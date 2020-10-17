% Function to output the dynamics matrices. Uses the lagrangian method.
% Inputs:
%   K: Kinetic Energy scalar
%   V: Potential Energy scalar
%   q: Generalized coordinates
%   dq: Time-derivative of the generalized coordinates
%   q_act: Actuated angles of the system
% Outputs:
%   D: D(q) Inertia matrix
%   C: C(q,dq) Coriollis matrix
%   G: G(q) Gravity matrix
%   B: B(q) Input Matrix?
function [D, C, G, B] = LagrangianDynamics(K, V, q, dq, q_act)

D = simplify( jacobian(jacobian(K,dq), dq) ) ;
for k=1:length(q)
    for j=1:length(q)
        C(k,j) = sym(0) ;
        for i=1:length(q)
            C(k,j) = C(k,j) + 1/2 * ( diff(D(k,j),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k)) ) * dq(i) ;
        end
    end
end
C = simplify(C) ;
G = simplify( jacobian(V,q) )' ;
B = jacobian(q_act, q)' ;