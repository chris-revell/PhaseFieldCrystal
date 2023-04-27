#https://en.wikipedia.org/wiki/Student%27s_t-test#Dependent_t-test_for_paired_samples

using CSV
using DataFrames
using Statistics

alldata = DataFrame(CSV.File("/Users/christopher/Dropbox (The University of Manchester)/Chris-Oliver Shared/PFC Project/defectCountsAll2.csv"))

function tTest(D1,D2)
    difference = D1.-D2
    x̄ = mean(difference)   # Mean
    s = std(difference)    # Sample standard deviation 
    n = length(difference)
    return x̄/(s/sqrt(n))    
end

@show tTest(alldata[!,:defectProportion], alldata[!,:defectProportion_1])



# (sumD/N)/(sqrt( (sumDsquared-(sumD^2)/N)/((N-1)*N) ))
# tTest(XminusYsum,XminusYsquaredSum,30)