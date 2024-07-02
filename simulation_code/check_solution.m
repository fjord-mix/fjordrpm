function status = check_solution(p, Q)

% CHECK_SOLUTION Check for errors in the solution.
%   STATUS = CHECK_SOLUTION(P, Q) checks that the volume of the layers is
%   conserved and returns an error status S with default value 0, or 1
%   if there was an error.

status = 0;

% Check that volume is conserved
if any(abs(Q.QVg-Q.QVs+Q.QVk+Q.QVmi+Q.QVv)>1e-10)
    disp("Volume not conserved")
    status = 1;
end

end