# avogadro constant (in mol^-1)
na = 6.0221415e23

# to convert particles/mum^3 to M
# litre = 1e-3 m^3 = 1e-3 (1e-6)^3 mum^3 = 1e-21 mum^3
# mol = 1/na
# 1 M = 1 mol/litre = 1/na/litre = 1/(na * 1e-21 mum^-3) = 
# [c] = M, [C]=particles/mum^3
# => c M = C particles/mum^3
# c = (toMolar) C
toMolar = na/1e21

# figure dimensions (according to PLoS One)
width1c = 3.27
height1c = width1c
width2c = 6.83
