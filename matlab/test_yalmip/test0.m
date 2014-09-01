x = sdpvar(2,1);
Constraints = [sum(x) <= 1, x(1)==0, x(2) >= 0.5];
Objective = x'*x;
options = sdpsettings('verbose',1,'solver','fmincon');
sol = solvesdp(Constraints,Objective,options);
if sol.problem == 0
 solution = double(x);
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end