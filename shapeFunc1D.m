function N = shapeFunc1D(var,position)
    if position == 1
        N = 1/2 * (1-var);
    elseif position == 2
        N = 1/2 * (1+var);
    end
end