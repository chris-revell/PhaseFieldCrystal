



function fibrilNumber()

    sums = Float64[]
    for ut in u
        push!(sums, sum(filter(x->x>0.5,ut)))
    end

    