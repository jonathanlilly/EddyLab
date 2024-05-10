function mustHaveSameRowLength(a,b)
    if ~isequal(size(a,1),size(b,1))
        eid = 'Size:notSameRowLength';
        msg = 'The number of rows in the two arguments must be the same.';
        error(eid,msg)
    end
end