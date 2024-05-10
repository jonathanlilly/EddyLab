function mustHaveSameSize(a,b)
    if ~isequal(size(a),size(b))
        eid = 'Size:notEqual';
        msg = 'The sizes of the two arguments must be the same.';
        error(eid,msg)
    end
end