function mustBeCompatible(x,y,z)
    if ~isequal(length(x),size(z,1))|~isequal(length(y),size(z,2))
        eid = 'Size:notCompatible';
        msg = 'Argument must have length(x) rows and length(y) columns.';
        error(eid,msg)
    end
end