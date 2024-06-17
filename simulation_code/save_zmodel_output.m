function save_zmodel_output(s, f, t, p, a, path_out)

    fjord_output.s = s;
    fjord_output.f = f;
    fjord_output.t = t;
    fjord_output.p = p;
    fjord_output.a = a;
    save(path_out,'fjord_output','-v7.3'); % v7.3 allows large files (> 2GB), which might happen in very long runs

end