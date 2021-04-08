function coordinates(Λ, data)
    x = data[!,"Points:0"]
    y = data[!,"Points:1"] .-  0.75inches
    z = data[!,"Points:2"] 

    vx = data[!,"U:0"] 
    vy = data[!,"U:1"] 
    vz = data[!,"U:2"] 

    r = sqrt.(y.^2.0 .+ z.^2.0)
    θ = atan.(z, y)
   
    @>> begin
        DataFrame(x = x, y = y, z = z, vx = -vx, vy = vy, vz = vz, θ = θ, r = r)
        filter(x-> (x[:r] >= Λ.r₁) & (x[:r] <= Λ.r₂))
    end
end

function convert_to_gridded_array(df)
    xs = @>> df x->x[!,:x] unique sort
    rs = @>> df x->x[!,:r] x->round.(x;digits = 5) unique sort
    θs = @>> df x->x[!,:θ] x->round.(x;digits = 4) unique sort
    
    vx = zeros(length(xs), length(rs), length(θs)) 
    vy = zeros(length(xs), length(rs), length(θs)) 
    vz = zeros(length(xs), length(rs), length(θs)) 

    map(eachrow(df)) do row
        i = argmin(abs.(row[:x] .- xs))
        j = argmin(abs.(row[:r] .- rs))
        k = argmin(abs.(row[:θ] .- θs))
        vx[i, j, k] = row[:vx] 
        vy[i, j, k] = row[:vy] 
        vz[i, j, k] = row[:vz] 
    end

    return xs, rs, θs, vx, vy, vz
end

rθtoyz(r,θ) = (r*cos(θ), r*sin(θ))
yztorθ(y,z) = (sqrt.(y.^2.0 .+ z.^2.0), atan.(z, y))

function vrθ(r, θ, vy, vz)     
    dt = 1e-6
    y, z = rθtoyz(r, θ)
    rn, θn = yztorθ(y + vy*dt, z + vz*dt)
    ((rn-r)/dt, (θn-θ)/dt)
end


function interpolated_field_functions()
    df = @>> begin
        CSV.read("Data/velocity.csv", DataFrame)
        coordinates(Λ)
    end

    xs, rs, θs, vx, vy, vz = @> convert_to_gridded_array(df)
    fvx = @> interpolate((xs, rs, θs), vx, Gridded(Linear())) extrapolate(0)
    fvy = @> interpolate((xs, rs, θs), vy, Gridded(Linear())) extrapolate(0)
    fvz = @> interpolate((xs, rs, θs), vz, Gridded(Linear())) extrapolate(0)

    function condition(θ) 
        p1(θ) = @match θ begin
            θ && if θ < -π end => θ + 2.0π
            θ && if θ > π end => θ - 2.0π
            _ => θ
        end
    
        p2(θ) = @match θ begin
            θ && if θ < θs[1] end => θs[1] 
            θ && if θ > θs[end] end => θs[end]
            _ => θ
        end
        
        @> θ p1 p2
    end 

    fx(x, r, θ) = fvx(x, r, condition(θ)) 

    function fr(x, r, θ)
        mt = condition(θ)
        vr, vθ = vrθ(r, θ, fvy(x, r, mt), fvz(x, r, mt))
        vr
    end

    function fθ(x, r, θ)
        mt = condition(θ)
        vr, vθ = vrθ(r, θ, fvy(x, r, mt), fvz(x, r, mt))
        vθ
    end

    return df, fx, fr, fθ
end

