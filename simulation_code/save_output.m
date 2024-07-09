function save_output(p, t, f, a, s, path_out)

% SAVE_OUTPUT Save the outputs of the z-model simulation.
%   SAVE_OUTPUT(p, t, f, a, s, path_out) saves a file containing s, f, t,
%   p, a in the location path_out.

    fjord_output.s = s;
    fjord_output.f = f;
    fjord_output.t = t;
    fjord_output.p = p;
    fjord_output.a = a;
    % v7.3 allows large files (> 2GB), possible in very long runs
    save(path_out,'fjord_output','-v7.3'); 

end