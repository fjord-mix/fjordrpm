function status = check_solution(Q)

% CHECK_SOLUTION Check for errors in the solution.
%   STATUS = CHECK_SOLUTION(Q) checks that the volume of the layers is
%   conserved and returns an error status S with default value 0, or 1
%   if there was an error.

status = 0;

% Check that volume is conserved
maxdVdt = max(abs(Q.QVg-Q.QVs+Q.QVk+Q.QVi+Q.QVv));
if maxdVdt>1e-10
    disp(['Warning: volume possibly not conserved, max dV/dt = ',...
           num2str(maxdVdt)]);
    status = 1;
end

end