function status = check_solution(i, s)

% CHECK_SOLUTION Check for errors in the solution.
%   status = CHECK_SOLUTION(i, s) checks that the volume of the layers is
%   conserved and returns an error status s with default value 0, or 1
%   if there was an error.

status = 0;

% get max volume change
maxdVdt = max(abs(sum(s.QVp(:,:,i),1)'+s.QVs(:,i)+s.QVk(:,i)+s.QVi(:,i)+s.QVv(:,i)+s.QVsurf(:,i)));

% change error status and display warning if above tolerance
if maxdVdt>1e-8
    disp(['Warning: volume possibly not conserved, max dV/dt = ',...
           num2str(maxdVdt)]);
    status = 1;
end

end