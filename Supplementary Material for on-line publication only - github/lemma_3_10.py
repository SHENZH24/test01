# This script is designed to verify the correctness of the calculations in Lemma 3.10 of the paper,
# allowing reviewers and readers to easily check and validate the results.
# The calculations are complex, involving judgments of polynomial parity and signs across multiple parameters.
# Including the full details in the paper would take up significant space and make verification more difficult.

import sympy as sp
# Define symbols
x = sp.symbols('x')
m, l, r, a = sp.symbols('m l r a', real=True, positive=True)
m_lower_bound = 2*l + 2*r + 5
r_lower_bound = 3
l_lower_bound = 1
a_lower_bound = 0

# Define matrix: matrix_lem_3_10
matrix_lem_3_10 = sp.Matrix([
    [0, 1, 2, r-2, 2, 0],
    [1, 0, 2, r-2, 0, 0],
    [1, 1, 1, 0, 0, 0],
    [1, 1, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, l],
    [0, 0, 0, 0, 2, 0]
])

# The characteristic polynomial of matrix_lem_3_10
char_poly = matrix_lem_3_10.charpoly(x).as_expr()

# \sqrt{\left(\frac{1+\sqrt{4m-7}}{2}\right)^2+a-(m-2r-2l-5)}
value_to_substitute = sp.sqrt((1 + sp.sqrt(4 * m - 7))**2 / 4 + a - (m - 2 * r - 2 * l - 5))

# The formular of \Phi_{M_2\left( \sqrt{\left(\frac{1+\sqrt{4m-7}}{2}\right)^2+a-(m-2r-2l-5)}\right)
poly_substituted = char_poly.subs(x, value_to_substitute)

# The formulas of psi_3 and psi_4 given in paper
psi_3 = (
    (6*a**2 + 16*a*l + 16*a*r + 34*a + 8*l**2 + 24*l*r + 44*l + 8*r**2 + 48*r + 2*m + 40)*sp.sqrt(4*m - 7)
    + 66*a + 108*l + 34*m + 104*r + 88*a*l + 12*a*m + 96*a*r + 16*l*m + 136*l*r + 16*m*r + 16*a*l**2
    + 16*a**2*l + 16*a*r**2 + 16*a**2*r + 32*l*r**2 + 32*l**2*r + 34*a**2 + 4*a**3 + 40*l**2 + 56*r**2 + 48*a*l*r
)

psi_4 = sp.sqrt(2) * ((4*a + 6*l + 6*r + 14)*sp.sqrt(4*m - 7) + 4*a**2 + 12*a*l + 12*a*r + 28*a + 8*l**2 + 16*l*r +
                       30*l + 8*r**2 + 50*r + 4*m + 26) * sp.sqrt(2*a + 4*l + 4*r + sp.sqrt(4*m - 7) + 7)

# Compare
result_1 = sp.simplify(poly_substituted - 1/4 * (psi_3 - psi_4))
# Output
print(f'The characteristic polynomial of matrix_lem_3_10: {char_poly}')
if result_1 == 0:
    print("\nIs the formula of Phi_M_2(value_to_substitute) in paper correct?: correct")
else:
    print("\nIs the formula of Phi_M_2(value_to_substitute) in paper correct?: error")

"""The following verifies that the calculations in the paper are correct."""
#-----------------------------------------------------------------------------------------------------------------
# Compute psi_3^2(m,l,r,a) - psi_4^2(m,l,r,a)
diff_psi_squared = sp.expand(psi_3**2 - psi_4**2)
#The formula of psi_3^2 - psi_4^2 given in paper
diff_psi_squared_coefficients_1 = ((48*a + 64*l + 64*r + 104)*m**2 + (160*a**3 + 640*a**2*l + 640*a**2*r +
    1040*a**2 + 768*a*l**2 + 1792*a*l*r + 2560*a*l + 768*a*r**2 + 2688*a*r + 1296*a + 256*l**3 +
    1152*l**2*r + 1312*l**2 + 1152*l*r**2 + 3520*l*r + 1312*l + 256*r**3 + 1504*r**2 + 1248*r - 832)*m +
    48*a**5 + 320*a**4*l + 320*a**4*r + 520*a**4 + 768*a**3*l**2 + 1792*a**3*l*r + 2560*a**3*l +
    768*a**3*r**2 + 2688*a**3*r + 1184*a**3 + 768*a**2*l**3 + 3456*a**2*l**2*r + 3936*a**2*l**2 +
    3456*a**2*l*r**2 + 10560*a**2*l*r + 3488*a**2*l + 768*a**2*r**3 + 4512*a**2*r**2 + 3296*a**2*r -
    3224*a**2 + 256*a*l**4 + 2560*a*l**3*r + 1664*a*l**3 + 4864*a*l**2*r**2 + 12032*a*l**2*r +
    1376*a*l**2 + 2560*a*l*r**3 + 13056*a*l*r**2 + 6976*a*l*r - 9056*a*l + 256*a*r**4 + 2432*a*r**3 +
    608*a*r**2 - 14880*a*r - 12864*a + 512*l**4*r - 256*l**4 + 2048*l**3*r**2 + 3328*l**3*r - 1728*l**3 +
    2048*l**2*r**3 + 8576*l**2*r**2 + 1728*l**2*r - 6816*l**2 + 512*l*r**4 + 4352*l*r**3 + 1216*l*r**2 -
    18880*l*r - 14080*l - 2496*r**3 - 16096*r**2 - 24320*r - 8800)

diff_psi_squared_coefficients_2 = (16*m**3 +
    (240*a**2 + 640*a*l + 640*a*r + 1040*a + 384*l**2 + 896*l*r + 1280*l + 384*r**2 + 1344*r + 648)*m**2 +
    (240*a**4 + 1280*a**3*l + 1280*a**3*r + 2080*a**3 + 2304*a**2*l**2 + 5376*a**2*l*r + 7680*a**2*l +
    2304*a**2*r**2 + 8064*a**2*r + 3552*a**2 + 1536*a*l**3 + 6912*a*l**2*r + 7872*a*l**2 + 6912*a*l*r**2 +
    21120*a*l*r + 6976*a*l + 1536*a*r**3 + 9024*a*r**2 + 6592*a*r - 6448*a + 256*l**4 + 2560*l**3*r +
    1664*l**3 + 4864*l**2*r**2 + 12032*l**2*r + 1376*l**2 + 2560*l*r**3 + 13056*l*r**2 + 6976*l*r -
    9056*l + 256*r**4 + 2432*r**3 + 608*r**2 - 14880*r - 12864)*m + 16*a**6 + 128*a**5*l + 128*a**5*r +
    208*a**5 + 384*a**4*l**2 + 896*a**4*l*r + 1280*a**4*l + 384*a**4*r**2 + 1344*a**4*r + 312*a**4 +
    512*a**3*l**3 + 2304*a**3*l**2*r + 2624*a**3*l**2 + 2304*a**3*l*r**2 + 7040*a**3*l*r + 832*a**3*l +
    512*a**3*r**3 + 3008*a**3*r**2 + 704*a**3*r - 4576*a**3 + 256*a**2*l**4 + 2560*a**2*l**3*r +
    1664*a**2*l**3 + 4864*a**2*l**2*r**2 + 12032*a**2*l**2*r - 1312*a**2*l**2 + 2560*a**2*l*r**3 +
    13056*a**2*l*r**2 + 704*a**2*l*r - 18016*a**2*l + 256*a**2*r**4 + 2432*a**2*r**3 - 2080*a**2*r**2 -
    24288*a**2*r - 17400*a**2 + 1024*a*l**4*r - 512*a*l**4 + 4096*a*l**3*r**2 + 6656*a*l**3*r -
    5248*a*l**3 + 4096*a*l**2*r**3 + 17152*a*l**2*r**2 - 4608*a*l**2*r - 22816*a*l**2 + 1024*a*l*r**4 +
    8704*a*l*r**3 - 5632*a*l*r**2 - 62400*a*l*r - 37344*a*l - 6784*a*r**3 - 42720*a*r**2 - 57376*a*r -
    11776*a - 512*l**5 + 1024*l**4*r**2 - 3584*l**4 + 2048*l**3*r**3 + 6144*l**3*r**2 - 6400*l**3*r -
    10176*l**3 + 1024*l**2*r**4 + 7168*l**2*r**3 - 4480*l**2*r**2 - 38208*l**2*r - 16160*l**2 +
    1024*l*r**4 - 7424*l*r**3 - 53696*l*r**2 - 61888*l*r - 8320*l - 512*r**5 - 4608*r**4 - 24896*r**3 -
    46944*r**2 - 17280*r + 8736)
diff_psi_squared_in_paper = (diff_psi_squared_coefficients_1*sp.sqrt(4*m - 7) + diff_psi_squared_coefficients_2)
result_2 = sp.simplify(diff_psi_squared - diff_psi_squared_in_paper)
if result_2 == 0:
    print("\nIs the formula of psi_3^2 - psi_4^2 in paper correct?: correct")
else:
    print("\nIs the formula of psi_3^2 - psi_4^2 in paper correct?: error")
# - The above result will show that diff_psi_squared = diff_psi_squared_in_paper
# - Both diff_psi_squared_coefficients_1 and diff_psi_squared_coefficients_2 can be regarded as functions of m.
# - We need to prove that both are monotonically increasing with respect to m.
# - diff_psi_squared_coefficients_1 is a quadratic function of m.
# - Taking its first derivative with respect to m, the derivative need to be greater than 0.
# - diff_psi_squared_coefficients_2 is a cubic function of m.
# - If both its first and second derivatives are greater than 0,
# - it will be monotonically increasing with respect to m.
diff_psi_squared_coefficients_1_derivative_m = sp.diff(diff_psi_squared_coefficients_1, m)
#From the following equation, we know that the diff_psi_squared_coefficients_1 is monotonically increasing with respect to m when m >= 2*r + 2*l + 5.
#It is easy to see that the value of the following polynomial is greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_1_derivative_m = sp.collect(sp.simplify(diff_psi_squared_coefficients_1_derivative_m.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_1_derivative_m: {substitute_lb_into_diff_psi_squared_coefficients_1_derivative_m}')

diff_psi_squared_coefficients_2_derivative_m_1 = sp.diff(diff_psi_squared_coefficients_2, m)
diff_psi_squared_coefficients_2_derivative_m_2 = sp.diff(diff_psi_squared_coefficients_2_derivative_m_1, m)
#From the following equation, we know that the diff_psi_squared_coefficients_2 is monotonically increasing with respect to m when m >= 2*r + 2*l + 5.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_1 = sp.collect(sp.simplify(diff_psi_squared_coefficients_2_derivative_m_1.subs(m, m_lower_bound)),a)
substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_2 = sp.collect(sp.simplify(diff_psi_squared_coefficients_2_derivative_m_2.subs(m, m_lower_bound)),a)
print(f'\nThe first derivative function of diff_psi_squared_coefficients_2 with respect to m: {sp.collect(sp.expand(diff_psi_squared_coefficients_2_derivative_m_1),m)}')
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_1: {substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_1}')
print(f'The second derivative function of diff_psi_squared_coefficients_2 with respect to m: {sp.collect(sp.expand(diff_psi_squared_coefficients_2_derivative_m_2),m)}')
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_2: {substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m_2}')

#----------------------------------------------------------------------------------------------------------------
#Substitute m_lower_bound = 2*l + 2*r + 2 into psi_3^2 - psi_4^2.
substitute_m_into_diff_psi_squared = diff_psi_squared.subs(m,m_lower_bound)
#The formula of substitute_m_into_diff_psi_squared given in paper
substitute_m_into_coefficients_1 = (((256*a + 512*l + 512)*r**4 +
    (768*a**2 + 2560*a*l + 3968*a + 2048*l**2 + 7168*l + 2048)*r**3 +
    (768*a**3 + 3456*a**2*l + 5792*a**2 + 4864*a*l**2 + 18176*a*l
     + 10016*a + 2048*l**3 + 13184*l**2 + 17792*l - 4384)*r**2 +
    (320*a**4 + 1792*a**3*l + 3008*a**3 + 3456*a**2*l**2 + 13120*a**2*l
    + 8576*a**2 + 2560*a*l**3 + 17152*a*l**2 + 26816*a*l + 2112*a + 512*l**4
    + 6144*l**3 + 17920*l**2 + 7232*l - 16064)*r + 48*a**5 + 320*a**4*l + 520*a**4
    + 768*a**3*l**2 + 2880*a**3*l + 1984*a**3 + 768*a**2*l**3 + 5216*a**2*l**2
    + 8768*a**2*l + 1976*a**2 + 256*a*l**4 + 3200*a*l**3 + 10528*a*l**2 + 7296*a*l
    - 5184*a + 256*l**4 + 2432*l**3 + 4064*l**2 - 5504*l - 10360))

substitute_m_into_coefficients_2 =((256*a**2 + 1024*a*l + 3072*a + 1024*l**2 + 6656*l + 3072)*r**4 +
    (512*a**3 + 2560*a**2*l + 7040*a**2 + 4096*a*l**2 + 25600*a*l + 21504*a + 2048*l**3 + 22016*l**2
    + 43008*l + 1664)*r**3 + (384*a**4 + 2304*a**3*l + 5568*a**3 + 4864*a**2*l**2 + 28416*a**2*l
    + 26528*a**2 + 4096*a*l**3 + 44800*a*l**2 + 96896*a*l + 32544*a + 1024*l**4 + 20992*l**3 + 80256*l**2
    + 68608*l - 33632)*r**2 + (128*a**5 + 896*a**4*l + 1824*a**4 + 2304*a**3*l**2 + 12160*a**3*l
    + 11264*a**3 + 2560*a**2*l**3 + 27392*a**2*l**2 + 60992*a**2*l + 27936*a**2 + 1024*a*l**4 + 23552*a*l**3
    + 95616*a*l**2 + 104256*a*l - 512*a + 5632*l**4 + 40448*l**3 + 80256*l**2 + 7104*l - 68448)*r
    + 16*a**6 + 128*a**5*l + 208*a**5 + 384*a**4*l**2 + 1760*a**4*l + 1512*a**4 + 512*a**3*l**3
    + 5184*a**3*l**2 + 11392*a**3*l + 5824*a**3 + 256*a**2*l**4 + 6272*a**2*l**3 + 26528*a**2*l**2
    + 32288*a**2*l + 6360*a**2 + 2560*a*l**4 + 20736*a*l**3 + 47456*a*l**2 + 21440*a*l - 18016*a
    + 2560*l**4 + 13824*l**3 + 11360*l**2 - 31968*l - 37384)
substitute_m_into_diff_psi_squared_in_paper = (substitute_m_into_coefficients_1*sp.sqrt(8*l + 8*r + 13) +substitute_m_into_coefficients_2)
# The above result will show that substitute_m_into_diff_psi_squared = substitute_m_into_diff_psi_squared_in_paper.
# Both substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 can be regarded as functions of r.
# We need to show that both are monotonically increasing with respect to r.
# It is evident that substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 are quartic functions of r.
# If the first, second, and third derivatives of these two polynomials with respect to r are all greater than 0,
# then substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 will be monotonically increasing with respect to r.
result_3 = sp.simplify(substitute_m_into_diff_psi_squared - substitute_m_into_diff_psi_squared_in_paper)
if result_3 == 0:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: error")

substitute_m_into_coefficients_1_derivative_r_1 = sp.diff(substitute_m_into_coefficients_1, r)
substitute_m_into_coefficients_1_derivative_r_2 = sp.diff(substitute_m_into_coefficients_1_derivative_r_1, r)
substitute_m_into_coefficients_1_derivative_r_3 = sp.diff(substitute_m_into_coefficients_1_derivative_r_2, r)
#From the following equation, we know that the substitute_m_into_coefficients_1 is monotonically increasing with respect to r when r >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_r_1.subs(r, r_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_2 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_r_2.subs(r, r_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_3 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_r_3.subs(r, r_lower_bound)),a)
print(f'The first derivative function of substitute_m_into_coefficients_1 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_r_1),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_1: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_1}')
print(f'The second derivative function of substitute_m_into_coefficients_1 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_r_2),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_2: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_2}')
print(f'The third derivative function of substitute_m_into_coefficients_1 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_r_3),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_3: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_3}')

substitute_m_into_coefficients_2_derivative_r_1 = sp.diff(substitute_m_into_coefficients_2, r)
substitute_m_into_coefficients_2_derivative_r_2 = sp.diff(substitute_m_into_coefficients_2_derivative_r_1, r)
substitute_m_into_coefficients_2_derivative_r_3 = sp.diff(substitute_m_into_coefficients_2_derivative_r_2, r)
#From the following equation, we know that the substitute_m_into_coefficients_2 is monotonically increasing with respect to r when r >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_r_1.subs(r, r_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_2 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_r_2.subs(r, r_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_3 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_r_3.subs(r, r_lower_bound)),a)
print(f'\nThe first derivative function of substitute_m_into_coefficients_2 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_r_1),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_1: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_1}')
print(f'The second derivative function of substitute_m_into_coefficients_2 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_r_2),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_2: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_2}')
print(f'The third derivative function of substitute_m_into_coefficients_1 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_r_3),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_3: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_3}')

#--------------------------------------------------------------------------------------------------------------------
#Substitute r_lower_bound = 3 into substitute_m_into_diff_psi_squared.
substitute_m_r_into_diff_psi_squared = substitute_m_into_diff_psi_squared.subs(r,r_lower_bound)
#The formula of substitute_m_r_into_diff_psi_squared given in paper
substitute_m_r_into_coefficients_1 = ((256*a + 1792)*l**4 + (768*a**2 + 10880*a + 39296)*l**3 +
     (768*a**3 + 15584*a**2 + 105760*a + 231776)*l**2 +
     (320*a**4 + 8256*a**3 + 79232*a**2 + 320448*a + 411328)*l +
     48*a**5 + 1480*a**4 + 17920*a**3 + 100568*a**2 + 219168*a - 1240)

substitute_m_r_into_coefficients_2 = ((256*a**2 + 5632*a + 28672)*l**4 + (512*a**3 + 13952*a**2 +
    128256*a + 379392)*l**3 + (384*a**4 + 12096*a**3 + 152480*a**2 + 848096*a + 1651808)*l**2 +
    (128*a**5 + 4448*a**4 + 68608*a**3 + 540128*a**2 + 1980416*a + 2307168)*l +
    16*a**6 + 592*a**5 + 10440*a**4 + 103552*a**3 + 539736*a**2 + 1102784*a - 251656)

substitute_m_r_into_diff_psi_squared_in_paper = (substitute_m_r_into_coefficients_1*sp.sqrt(8*l + 37)+substitute_m_r_into_coefficients_2)
result_4 = sp.simplify(substitute_m_r_into_diff_psi_squared - substitute_m_r_into_diff_psi_squared_in_paper)
if result_4 == 0:
    print("\nIs the formula of substitute_m_r_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_r_into_diff_psi_squared in paper correct?: error")
# The above result will show that substitute_m_r_into_diff_psi_squared = substitute_m_r_into_diff_psi_squared_in_paper
# Both substitute_m_r_into_coefficients_1 and substitute_m_r_into_coefficients_2 can be regarded as functions of l.
# We need to show that both are monotonically increasing with respect to l.
# It is evident that substitute_m_r_into_coefficients_1 and substitute_m_r_into_coefficients_2 are quartic functions of l.
# If the first, second, and third derivatives of these two polynomials with respect to l are all greater than 0,
# then substitute_m_r_into_coefficients_1 and substitute_m_r_into_coefficients_2 will be monotonically increasing with respect to l.

substitute_m_r_into_coefficients_1_derivative_l_1 = sp.diff(substitute_m_r_into_coefficients_1, l)
substitute_m_r_into_coefficients_1_derivative_l_2 = sp.diff(substitute_m_r_into_coefficients_1_derivative_l_1, l)
substitute_m_r_into_coefficients_1_derivative_l_3 = sp.diff(substitute_m_r_into_coefficients_1_derivative_l_2, l)
#From the following equation, we know that the substitute_m_r_into_coefficients_1 is monotonically increasing with respect to l when l >= 1.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_l_1 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_1_derivative_l_1.subs(l, l_lower_bound)),a)
substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_l_2 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_1_derivative_l_2.subs(l, l_lower_bound)),a)
substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_l_3 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_1_derivative_l_3.subs(l, l_lower_bound)),a)
print(f'The first derivative function of substitute_m_r_into_coefficients_1 with respect to l: {sp.collect(sp.expand(substitute_m_r_into_coefficients_1_derivative_l_1),l)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_l_1: {substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_l_1}')
print(f'The second derivative function of substitute_m_r_into_coefficients_1 with respect to l: {sp.collect(sp.expand(substitute_m_r_into_coefficients_1_derivative_l_2),l)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_l_2: {substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_l_2}')
print(f'The third derivative function of substitute_m_r_into_coefficients_1 with respect to l: {sp.collect(sp.expand(substitute_m_r_into_coefficients_1_derivative_l_3),l)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_l_3: {substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_l_3}')

substitute_m_r_into_coefficients_2_derivative_l_1 = sp.diff(substitute_m_r_into_coefficients_2, l)
substitute_m_r_into_coefficients_2_derivative_l_2 = sp.diff(substitute_m_r_into_coefficients_2_derivative_l_1, l)
substitute_m_r_into_coefficients_2_derivative_l_3 = sp.diff(substitute_m_r_into_coefficients_2_derivative_l_2, l)
#From the following equation, we know that the substitute_m_r_into_coefficients_2 is monotonically increasing with respect to l when l >= 1.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_l_1 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_2_derivative_l_1.subs(l, l_lower_bound)),a)
substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_l_2 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_2_derivative_l_2.subs(l, l_lower_bound)),a)
substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_l_3 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_2_derivative_l_3.subs(l, l_lower_bound)),a)
print(f'\nThe first derivative function of substitute_m_r_into_coefficients_1 with respect to l: {sp.collect(sp.expand(substitute_m_r_into_coefficients_2_derivative_l_1),l)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_l_1: {substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_l_1}')
print(f'The second derivative function of substitute_m_r_into_coefficients_1 with respect to l: {sp.collect(sp.expand(substitute_m_r_into_coefficients_2_derivative_l_2),l)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_l_2: {substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_l_2}')
print(f'The third derivative function of substitute_m_r_into_coefficients_1 with respect to l: {sp.collect(sp.expand(substitute_m_r_into_coefficients_2_derivative_l_3),l)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_l_3: {substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_l_3}')


#--------------------------------------------------------------------------------------------------------------
#Substitute l_lower_bound = 1 into substitute_m_r_into_diff_psi_squared.
substitute_m_r_l_into_diff_psi_squared = substitute_m_r_into_diff_psi_squared.subs(l,l_lower_bound)
#The formula of substitute_m_r_l_into_diff_psi_squared given in paper
substitute_m_r_l_into_diff_psi_squared_in_paper = (
    16*a**6 +(144*sp.sqrt(5) + 720)*a**5 + (5400*sp.sqrt(5) + 15272)*a**4 + (80832*sp.sqrt(5) + 184768)*a**3 +
    (588456*sp.sqrt(5) + 1246552)*a**2 + (1969536*sp.sqrt(5) + 4065184)*a + 2048856*sp.sqrt(5) + 4115384
)
result_5 = sp.simplify(substitute_m_r_l_into_diff_psi_squared - substitute_m_r_l_into_diff_psi_squared_in_paper)
if result_5 == 0:
    print("\nIs the formula of substitute_m_r_l_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_r_l_into_diff_psi_squared in paper correct?: error")
