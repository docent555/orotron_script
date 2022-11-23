classdef hField < handle
    properties
        Field, Zeta
    end
    methods
        function h = hField(nr, nc)
            h.Field = zeros(nr, 1);
            h.Zeta = zeros(nr, 1);
        end
    end
end
