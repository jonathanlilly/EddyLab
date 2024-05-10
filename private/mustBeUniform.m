function mustBeUniform(x)
    if ~isuniform(x)
        eid = 'Spacing:notUniform';
        msg = 'Argument must be uniformly spaced.';
        error(eid,msg)
    end
end