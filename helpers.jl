using ForneyLab
using CSV, DataFrames
import ForneyLab: ruleSVBGaussianMeanPrecisionMGVD, ruleMGaussianMeanPrecisionGGD, 
                  @symmetrical, prod!, unsafeWeightedMean, unsafePrecision

ruleSPMultiplicationAGPN(msg_out::Message{F, Multivariate}, msg_in1::Message{PointMass, Multivariate}, msg_a::Nothing) where F<:Gaussian = ruleSPMultiplicationIn1GNP(msg_out, nothing, msg_in1)

function ruleSVBGaussianMeanPrecisionMGVD(msg_out::Message{GaussianMeanVariance, Univariate},
                                          dist_mean::Any,
                                          dist_prec::ProbabilityDistribution)
    d_out = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_out.dist)
    Message(Multivariate, GaussianMeanVariance, m=[d_out.params[:m]], v=[d_out.params[:v]] + cholinv(unsafeMean(dist_prec)))
end

function ruleMGaussianMeanPrecisionGGD(msg_out::Message{GaussianMeanVariance,Univariate},
                                       msg_mean::Message{GaussianWeightedMeanPrecision,Multivariate},
                                       dist_prec::ProbabilityDistribution) where {F1<:Gaussian, F2<:Gaussian}

    m = [msg_out.dist.params[:m]]
    v = mat(msg_out.dist.params[:v])
    
    msg_out = Message(Multivariate, GaussianMeanVariance, m=m, v=v)
    return ruleMGaussianMeanPrecisionGGD(msg_out, msg_mean, dist_prec)
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Multivariate, F1},
    y::ProbabilityDistribution{Univariate, F2},
    z::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}=ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[NaN], w=transpose([NaN]))) where {F1<:Gaussian, F2<:Gaussian}

    if dims(x) == 1
        z.params[:xi] = unsafeWeightedMean(x) + [unsafeWeightedMean(y)]
        z.params[:w] = unsafePrecision(x) + [unsafePrecision(y)]
        return z
    else 
        throw(DimensionMismatch([x, y]))
    end
end

function beautify_process(process_set)
    dump = Vector{Vector{Float64}}()
    for entry in process_set
        tmp = Vector{Float64}()
        for entity in entry
            push!(tmp, entity[1])
        end
        push!(dump, tmp)
    end
    process_set = dump
    return process_set
end

unzip(a, index) = map(d -> d[index], a)

logPDF(mx, x, vx) = log(vx) + ((mx - x)^2)/(vx)


# Filter bad data
FloatParse(x) =
try
    parse(Float64, x)
    return true
catch
    return false
end

function loadAR(filepath::String; col, delim=',')
    df = CSV.File(filepath, delim=delim) |> DataFrame
    x = []
    df = DataFrame(value=df[!, col])
    # Data
    for i in range(1, stop=size(df, 1))
        if typeof(df[i, 1]) == String && FloatParse(df[i, 1])
            xi = parse(Float64, df[i, 1])
            push!(x, xi)
        elseif typeof(df[i, 1]) == Float64
            xi = convert(Float64, df[i, 1])
            push!(x, xi)
        end
    end
    return x
end