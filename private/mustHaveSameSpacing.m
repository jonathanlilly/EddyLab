function mustHaveSameSpacing(x,y)
    if (x(2)-x(1))~=(y(2)-y(1))
        eid = 'Spacing:notSame';
        msg = 'Arguments x and y must have the same spacing or sampling interval.';
        error(eid,msg)
    end
end