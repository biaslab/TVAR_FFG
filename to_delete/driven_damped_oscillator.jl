# Pick sampling time
δ = 0.01

# Length of time-series
T = 3000

# Generate time range
TT = δ.*range(0, length=T);

# Pick dynamical parameters
m = 2.
k = 0.8

# Start state array
x = zeros(T,)

# Pick initial condition
x[1] = 1.
x[2] = 1.

# Generate input signal
A = 1.2
ω = 0.5
u = A.*sin.(ω .*TT)

# Pick dynamical parameters
m = 2.
c = 0.4
k = 0.8

# Start state array
x = zeros(T,)

# Pick initial condition
x[1] = 0.
x[2] = 0.

# State transition function
driven_damped_harmonic_oscillator(x_t, x_tmin1, m, c, k, δ, u_t) =  (2*m - c*δ - k*δ^2)/m*x_t + (c*δ - m)/m*x_tmin1 + δ^2*u_t

# Iterate
for t = 2:T-1
   x[t+1] = driven_damped_harmonic_oscillator(x[t], x[t-1], m, c, k, δ, u[t])
end

