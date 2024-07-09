function status = check_solution(i, s)

% CHECK_SOLUTION Check for errors in the solution.
%   status = CHECK_SOLUTION(i, s) checks that the volume of the layers is
%   conserved and returns an error status s with default value 0, or 1
%   if there was an error.

status = 0;

% Check that volume is conserved
maxdVdt = max(abs(s.QVp(:,i)+s.QVs(:,i)+s.QVk(:,i)+s.QVi(:,i)+s.QVv(:,i)));
if maxdVdt>1e-10
    disp(['Warning: volume possibly not conserved, max dV/dt = ',...
           num2str(maxdVdt)]);
    status = 1;
end

end