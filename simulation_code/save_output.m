function save_zmodel_output(s, f, t, p, a, path_out)

% SAVE_ZMODEL_OUTPUTS Save the outputs of the z-model simulation.
%   SAVE_ZMODEL_OUTPUTS(S, F, T, P, A, PATH_OUT) saves a file containing S,
%   F, T, P, A in the location PATH_OUT.

    fjord_output.s = s;
    fjord_output.f = f;
    fjord_output.t = t;
    fjord_output.p = p;
    fjord_output.a = a;
    % v7.3 allows large files (> 2GB), which might happen in very long
    % runs.
    save(path_out,'fjord_output','-v7.3'); 

end