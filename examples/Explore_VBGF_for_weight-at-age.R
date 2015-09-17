

# Parameters (to be found for life history database)
Linf = 100 # cm
k_L = 0.1 # growth coefficient for VBGF in length
L0 = 10
a = 0.01 # g per cm^3
b = 2 # isometric weight at age (DEFAULT = 3)
r = 3 # age at recruitment (whole number)

# Calculate length at-age
a_set = 0:30
L_a = L0 + (Linf-L0) * (1-exp( -k_L*a_set))

# Calculate parameters for VBGF in weight
w_r = a * L_a[which(a_set==r)]^b
Winf = a * Linf^b
k_W = k_L  # I'm not sure why this isn't k_W = 3*k_L, but that doesn't seem to work well
rho = exp( -k_W )
alpha_g = Winf * (1-rho)

# Calculate weight at age from VBGF for weight at age using recursive parameterization used by Delay-Difference model
W_a = rep(NA,length(a_set))
W_a[which(a_set==r)] = w_r
for( i in (which(a_set==r)+1):length(a_set)) W_a[i] = alpha_g + rho*W_a[i-1]

# Confirm that they're similar past age r
matplot( cbind( a*L_a^b, W_a), type="l", col=c("black","blue"), lty="solid", lwd=3)

