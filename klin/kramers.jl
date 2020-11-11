######################################################
## Kramers example

module kramers

## Kramers drift field
function Drift(;A=1.0,B=1.0)
    # V(x)  = B*(x^2-A)^2/4
    # DV(x) = B*(x^2-A)*x
    function(x)
        -B*(x^2-A)*x
    end
end

end#module
