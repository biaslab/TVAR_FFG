using LAR
using LAR.Data
using Random
import Distributions: Normal, MvNormal

Random.seed!(42)

n_samples = 100
n_dataset = 100
max_order = 10
order_set = div(n_dataset, max_order)

# auxilary function to check if time-series converge
function infinite_ar(ar_array; lim=10.0)
    for x in ar_array
        if findnext(x -> abs(x) > 10.0, x, 1) !== nothing
            return false
        end
    end
    return true
end

# time-varying AR (RW prior)
function generateTVAR(num::Int, order::Int; nvar=1.0)
    coefs = Vector{Vector{Float64}}(undef, num+3*order)
    coef_var = rand([1.0, 1e-1, 1e-2])
    coefs[1] = Data.generate_coefficients(order)
    inits = randn(order)
    data = Vector{Vector{Float64}}(undef, num+3*order)
    data[1] = inits
    for i in 2:num+3*order
        coefs[i] = rand(MvNormal(coefs[i-1], coef_var))
        data[i] = insert!(data[i-1][1:end-1], 1, rand(Normal(coefs[i]'data[i-1], sqrt(nvar)), 1)[1])
    end
    data = data[1+3*order:end]
    coefs = coefs[1+3*order:end]
    coef_var = coef_var
    return coefs, coef_var, data
end

# statitionary AR dataset
ar_data_process = []
ar_data_observations = []
ar_coefs_set = Vector{Vector{Float64}}(undef, n_dataset)
mnv = 1.0
i = 1
ar_coefs_set[i] = Data.generate_coefficients(div(i-1, order_set) + 1)
while i <= n_dataset
    global i
    tmp = generateAR(n_samples, length(ar_coefs_set[i]), nvar=1.0, coefs=ar_coefs_set[i])[2]
    if infinite_ar(tmp, lim=10.0)
        push!(ar_data_process, tmp)
        # Observations
        push!(ar_data_observations, [x[1] + sqrt(mnv)*randn() for x in tmp])
        i += 1
        if i > n_dataset
            break
        end
        ar_coefs_set[i] = Data.generate_coefficients(div(i-1, order_set) + 1)
    end
end

# time-varying AR dataset
tvar_data_process = []
tvar_data_observations = []
tvar_coefs_set = []
tvar_var_set = []
mnv = 1.0
i = 1
while i <= n_dataset
    global i
    tmp_coefs, tmp_var, tmp = generateTVAR(n_samples, div(i-1, order_set) + 1, nvar=1.0)
    if infinite_ar(tmp, lim=100.0)
        push!(tvar_data_process, tmp)
        # Observations
        push!(tvar_data_observations, [x[1] + sqrt(mnv)*randn() for x in tmp])
        push!(tvar_coefs_set, tmp_coefs)
        push!(tvar_var_set, tmp_var)
        i += 1
    end
end