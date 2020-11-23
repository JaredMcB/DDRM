######################################################
## linear example

module linear

## Kramers drift field
function Drift(;A=1.0)
    function(x)
        -A*x
    end
end

end#module
