function mustBeCompatible(x,y,z)
    if ~isequal(length(y),size(z,1))|~isequal(length(x),size(z,2))
        eid = 'Size:notCompatible';
        msg = 'Argument must have length(y) rows and length(x) columns.';
        error(eid,msg)
    end
end