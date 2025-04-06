# This script is designed to verify the correctness of the calculations in Lemma 3.11 of the paper,
# allowing reviewers and readers to easily check and validate the results.
# The calculations are complex, involving judgments of polynomial parity and signs across multiple parameters.
# Including the full details in the paper would take up significant space and make verification more difficult.

import sympy as sp

# Define symbols
x = sp.symbols('x')
m, r, t, a = sp.symbols('m r t a', real=True)
m_lower_bound = 2*r + t + 8
r_lower_bound = 3
t_lower_bound = 3
a_lower_bound = 0

# Define matrix: matrix_lem_3_11
matrix_lem_3_11 = sp.Matrix([
    [0, 1, 2, r-2, t],
    [1, 0, 2, r-2, 0],
    [1, 1, 1, 0, 0],
    [1, 1, 0, 0, 0],
    [1, 0, 0, 0, 0]
])
# The characteristic polynomial of matrix_lem_3_11
char_poly = matrix_lem_3_11.charpoly(x).as_expr()

# Define the expression for x
value_to_substitute = sp.sqrt(((1 + sp.sqrt(4*m - 7)) / 2)**2 + a - (m - 2*r - t - 2))

# The formular of \Phi_{M_3}\left( \sqrt{\left( \frac{1 + \sqrt{4m - 7}}{2} \right)^2 + a - (m - 2r - t - 2)} \right)
poly_substituted = char_poly.subs(x, value_to_substitute)

# Define psi_5 and psi_6 expressions
psi_5 = sp.sqrt(2) * ((4*a + 4*r + 2*t)*sp.sqrt(4*m - 7) + 4*m + 4*r - 2*t + 8*a*r + 4*a*t + 4*r*t + 4*a**2 - 24) \
        * sp.sqrt(2*a + 4*r + 2*t + sp.sqrt(4*m - 7) + 1)

psi_6 = (8*a + 16*r + 4*t + 16)*sp.sqrt(4*m - 7) + 32*a + 8*m + 64*r + 12*t + 32*a*r + 8*a*t + 24*r*t + 8*a**2 + 32*r**2

#Compare
result_1 = sp.simplify(poly_substituted - 1/8 * (psi_5 - psi_6))

# Output
print(f'The characteristic polynomial of matrix_lem_3_11: {char_poly}')
if result_1 == 0:
    print("\nIs the formula of Phi_{M_2}(value_to_substitute) in paper correct?: correct")
else:
    print("\nIs the formula of Phi_{M_2}(value_to_substitute) in paper correct?: error")

"""The following verifies that the calculations in the paper are correct."""
# Compute psi_5^2(m,l,r,a) - psi_6^2(m,l,r,a)
diff_psi_squared = sp.expand(psi_5**2 - psi_6**2)
#The formula of psi_5^2 - psi_6^2 given in paper
diff_psi_squared_coefficients_1= (32*m**2 + (384*r**2 + (768*a + 448*t - 128)*r + 320*a**2 + 384*a*t - 64*a + 96*t**2 - 64*t - 640)*m
              + (512*a + 256*t - 768)*r**3 + (1152*a**2 + 1280*a*t - 896*a + 288*t**2 - 768*t - 4736)*r**2
              + (768*a**3 + 1344*a**2*t - 384*a**2 + 640*a*t**2 - 704*a*t - 6528*a + 64*t**3 - 256*t**2 - 3840*t - 2816)*r
              + 160*a**4 + 384*a**3*t - 64*a**3 + 288*a**2*t**2 - 192*a**2*t - 2144*a**2
              + 64*a*t**3 - 160*a*t**2 - 2496*a*t - 1408*a - 32*t**3 - 544*t**2 - 384*t + 1152)

diff_psi_squared_coefficients_2 = (320*a + 384*r + 192*t - 32)*m**2 + (640*a**3 + 2304*a**2*r + 1152*a**2*t - 192*a**2
              + 2304*a*r**2 + 2688*a*r*t - 768*a*r + 576*a*t**2 - 384*a*t - 4288*a + 512*r**3
              + 1280*r**2*t - 896*r**2 + 640*r*t**2 - 704*r*t - 6528*r + 64*t**3 - 160*t**2 - 2496*t - 1408)*m \
              + 64*a**5 + 384*a**4*r + 192*a**4*t - 32*a**4 + 768*a**3*r**2 + 896*a**3*r*t - 256*a**3*r \
              + 192*a**3*t**2 - 128*a**3*t - 2176*a**3 + 512*a**2*r**3 + 1280*a**2*r**2*t - 896*a**2*r**2 \
              + 640*a**2*r*t**2 - 704*a**2*r*t - 9216*a**2*r + 64*a**2*t**3 - 160*a**2*t**2 - 3840*a**2*t - 1184*a**2 \
              + 512*a*r**3*t - 1536*a*r**3 + 576*a*r**2*t**2 - 1536*a*r**2*t - 12160*a*r**2 + 128*a*r*t**3 - 512*a*r*t**2 \
              - 10816*a*r*t - 4736*a*r - 64*a*t**3 - 1760*a*t**2 - 320*a*t + 6784*a - 1024*r**4 + 128*r**3*t**2 \
              - 1280*r**3*t - 4864*r**3 + 64*r**2*t**3 - 544*r**2*t**2 - 7168*r**2*t - 4480*r**2 - 64*r*t**3 - 2304*r*t**2 \
              - 1280*r*t + 10496*r - 96*t**3 + 416*t**2 + 4736*t + 2944
diff_psi_squared_in_paper = (diff_psi_squared_coefficients_1*sp.sqrt(4*m - 7) + diff_psi_squared_coefficients_2)
result_2 = sp.simplify(diff_psi_squared - diff_psi_squared_in_paper)
if result_2 == 0:
    print("\nIs the formula of psi_5^2 - psi_6^2 in paper correct?: correct")
else:
    print("\nIs the formula of psi_5^2 - psi_6^2 in paper correct?: error")
#The above result will show that diff_psi_squared = diff_psi_squared_in_paper
# - Both diff_psi_squared_coefficients_1 and diff_psi_squared_coefficients_2 can be regarded as functions of m.
# - We need to prove that both are monotonically increasing with respect to m.
# - diff_psi_squared_coefficients_1 and diff_psi_squared_coefficients_2 are both quadratic functions of m.
# - Taking their first derivative with respect to m, the derivatives need to be greater than 0.
diff_psi_squared_coefficients_1_derivative_m = sp.diff(diff_psi_squared_coefficients_1, m)
#From the following equation, we know that the diff_psi_squared_coefficients_1 is monotonically increasing with respect to m when m >= 2*r + t + 8.
#It is easy to see that the value of the following polynomial is greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_1_derivative_m = sp.collect(sp.simplify(diff_psi_squared_coefficients_1_derivative_m.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_1_derivative_m: {substitute_lb_into_diff_psi_squared_coefficients_1_derivative_m}')

diff_psi_squared_coefficients_2_derivative_m = sp.diff(diff_psi_squared_coefficients_2, m)
#From the following equation, we know that the diff_psi_squared_coefficients_2 is monotonically increasing with respect to m when m >= 2*r + t + 8.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m = sp.collect(sp.simplify(diff_psi_squared_coefficients_2_derivative_m.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m: {substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m}')

#----------------------------------------------------------------------------------------------------------------
#Substitute m_lower_bound = 2*r + t + 8 into psi_5^2 - psi_6^2.
substitute_m_into_diff_psi_squared = diff_psi_squared.subs(m,m_lower_bound)
#The formula of substitute_m_into_diff_psi_squared given in paper
substitute_m_into_coefficients_1 = (((512*a + 256*t)*r**3 + (1152*a**2 + 1280*a*t + 640*a + 288*t**2
              + 512*t - 1792)*r**2 + (768*a**3 + 1344*a**2*t + 256*a**2 + 640*a*t**2 + 832*a*t
              - 512*a + 64*t**3 + 384*t**2 - 384*t - 4096)*r + 160*a**4 + 384*a**3*t - 64*a**3
              + 288*a**2*t**2 + 128*a**2*t + 416*a**2 + 64*a*t**3 + 224*a*t**2 + 512*a*t - 1920*a
              + 64*t**3 + 192*t**2 - 1024*t - 1920))

substitute_m_into_coefficients_2 = (((512*a**2 + 512*a*t + 3072*a + 128*t**2 + 1792*t - 1024)*r**3)
              + (768*a**3 + 1280*a**2*t + 3712*a**2 + 576*a*t**2 + 6144*a*t + 6016*a
              + 64*t**3 + 2016*t**2 + 3072*t - 12544)*r**2  + (384*a**4 + 896*a**3*t + 1024*a**3
              + 640*a**2*t**2 + 3904*a**2*t + 8832*a**2 + 128*a*t**3 + 3328*a*t**2 + 10432*a*t
              - 9216*a + 704*t**3 + 2944*t**2 - 6272*t - 20992)*r + 64*a**5 + 192*a**4*t - 32*a**4
              + 192*a**3*t**2 + 512*a**3*t + 2944*a**3 + 64*a**2*t**3 + 992*a**2*t**2 + 5184*a**2*t
              - 2720*a**2 + 512*a*t**3 + 2784*a*t**2 - 2560*a*t - 7040*a + 64*t**4 + 448*t**3 - 320*t**2
              - 4864*t - 10368)

substitute_m_into_diff_psi_squared_in_paper = (substitute_m_into_coefficients_1*sp.sqrt(8*r + 4*t + 25) +substitute_m_into_coefficients_2)
# The above result will show that substitute_m_into_diff_psi_squared = substitute_m_into_diff_psi_squared_in_paper
# Both substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 can be regarded as functions of r.
# We need to show that both are monotonically increasing with respect to r.
# It is evident that substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 are cubic functions of r.
# If the first and second derivatives of these two polynomials with respect to r are all greater than 0,
# then substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 will be monotonically increasing with respect to r.
result_3 = sp.simplify(substitute_m_into_diff_psi_squared - substitute_m_into_diff_psi_squared_in_paper)
if result_3 == 0:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: error")

substitute_m_into_coefficients_1_derivative_r_1 = sp.diff(substitute_m_into_coefficients_1, r)
substitute_m_into_coefficients_1_derivative_r_2 = sp.diff(substitute_m_into_coefficients_1_derivative_r_1, r)
#From the following equation, we know that the substitute_m_into_coefficients_1 is monotonically increasing with respect to r when r >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_r_1.subs(r, r_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_2 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_r_2.subs(r, r_lower_bound)),a)
print(f'The first derivative function of substitute_m_into_coefficients_1 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_r_1),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_1: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_1}')
print(f'The second derivative function of substitute_m_into_coefficients_1 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_r_2),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_2: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_r_2}')

substitute_m_into_coefficients_2_derivative_r_1 = sp.diff(substitute_m_into_coefficients_2, r)
substitute_m_into_coefficients_2_derivative_r_2 = sp.diff(substitute_m_into_coefficients_2_derivative_r_1, r)
#From the following equation, we know that the substitute_m_into_coefficients_2 is monotonically increasing with respect to r when r >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_r_1.subs(r, r_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_2 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_r_2.subs(r, r_lower_bound)),a)
print(f'The first derivative function of substitute_m_into_coefficients_2 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_r_1),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_1: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_1}')
print(f'The second derivative function of substitute_m_into_coefficients_2 with respect to r: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_r_2),r)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_2: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_r_2}')

#--------------------------------------------------------------------------------------------------------------------
#Substitute r_lower_bound = 3 into substitute_m_into_diff_psi_squared.
substitute_m_r_into_diff_psi_squared = substitute_m_into_diff_psi_squared.subs(r,r_lower_bound)
#The formula of substitute_m_r_into_diff_psi_squared given in paper
substitute_m_r_into_coefficients_1 = ((64*a + 256)*t**3 + (288*a**2 + 2144*a + 3936)*t**2
              + (384*a**3 + 4160*a**2 + 14528*a + 9344)*t + 160*a**4 + 2240*a**3
              + 11552*a**2 + 16128*a - 30336)

substitute_m_r_into_coefficients_2 = (64*t**4 + (64*a**2 + 896*a + 3136)*t**3
              + (192*a**3 + 2912*a**2 + 17952*a + 30112)*t**2
              + (192*a**4 + 3200*a**3 + 28416*a**2 + 97856*a + 52352)*t
              + 64*a**5 + 1120*a**4 + 12928*a**3 + 71008*a**2 + 102400*a - 213888)

substitute_m_r_into_diff_psi_squared_in_paper = (substitute_m_r_into_coefficients_1*sp.sqrt(4*t + 49)+substitute_m_r_into_coefficients_2)
result_4 = sp.simplify(substitute_m_r_into_diff_psi_squared - substitute_m_r_into_diff_psi_squared_in_paper)
if result_4 == 0:
    print("\nIs the formula of substitute_m_r_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_r_into_diff_psi_squared in paper correct?: error")
# The above result will show that substitute_m_r_into_diff_psi_squared = substitute_m_r_into_diff_psi_squared_in_paper
# Both substitute_m_r_into_coefficients_1 and substitute_m_r_into_coefficients_2 can be regarded as functions of t.
# We need to show that both are monotonically increasing with respect to t.
# It is evident that substitute_m_r_into_coefficients_1 is a cubic function of t,
# and substitute_m_r_into_coefficients_2 is a quartic function of t.
# If the first and second derivatives of substitute_m_r_into_coefficients_1 with respect to t,
# and the first, second and third derivatives of substitute_m_r_into_coefficients_1 with respect to t
# are all greater than 0,
# then substitute_m_r_into_coefficients_1 and substitute_m_r_into_coefficients_2 will be monotonically increasing with respect to t.
substitute_m_r_into_coefficients_1_derivative_t_1 = sp.diff(substitute_m_r_into_coefficients_1, t)
substitute_m_r_into_coefficients_1_derivative_t_2 = sp.diff(substitute_m_r_into_coefficients_1_derivative_t_1, t)
#From the following equation, we know that the substitute_m_r_into_coefficients_1 is monotonically increasing with respect to t when t >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_t_1 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_1_derivative_t_1.subs(t, t_lower_bound)),a)
substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_t_2 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_1_derivative_t_2.subs(t, t_lower_bound)),a)
print(f'The first derivative function of substitute_m_r_into_coefficients_1 with respect to t: {sp.collect(sp.expand(substitute_m_r_into_coefficients_1_derivative_t_1),t)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_t_1: {substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_t_1}')
print(f'The second derivative function of substitute_m_r_into_coefficients_1 with respect to t: {sp.collect(sp.expand(substitute_m_r_into_coefficients_1_derivative_t_2),t)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_t_2: {substitute_lb_into_substitute_m_r_into_coefficients_1_derivative_t_2}')

substitute_m_r_into_coefficients_2_derivative_t_1 = sp.diff(substitute_m_r_into_coefficients_2, t)
substitute_m_r_into_coefficients_2_derivative_t_2 = sp.diff(substitute_m_r_into_coefficients_2_derivative_t_1, t)
substitute_m_r_into_coefficients_2_derivative_t_3 = sp.diff(substitute_m_r_into_coefficients_2_derivative_t_2, t)
#From the following equation, we know that the substitute_m_r_into_coefficients_2 is monotonically increasing with respect to t when t >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_t_1 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_2_derivative_t_1.subs(t, t_lower_bound)),a)
substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_t_2 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_2_derivative_t_2.subs(t, t_lower_bound)),a)
substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_t_3 = sp.collect(sp.simplify(substitute_m_r_into_coefficients_2_derivative_t_3.subs(t, t_lower_bound)),a)
print(f'\nThe first derivative function of substitute_m_r_into_coefficients_2 with respect to t: {sp.collect(sp.expand(substitute_m_r_into_coefficients_2_derivative_t_1),t)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_t_1: {substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_t_1}')
print(f'The second derivative function of substitute_m_r_into_coefficients_2 with respect to t: {sp.collect(sp.expand(substitute_m_r_into_coefficients_2_derivative_t_2),t)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_t_2: {substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_t_2}')
print(f'The third derivative function of substitute_m_r_into_coefficients_2 with respect to t: {sp.collect(sp.expand(substitute_m_r_into_coefficients_2_derivative_t_3),t)}')
print(f'The value of substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_t_3: {substitute_lb_into_substitute_m_r_into_coefficients_2_derivative_t_3}')

#--------------------------------------------------------------------------------------------------------------
#Substitute t_lower_bound = 3 into substitute_m_r_into_diff_psi_squared.
substitute_m_r_t_into_diff_psi_squared = substitute_m_r_into_diff_psi_squared.subs(t,t_lower_bound)
#The formula of substitute_m_r_t_into_diff_psi_squared given in paper
substitute_m_r_t_into_diff_psi_squared_in_paper = ((581728 * a + sp.sqrt(61) * (160 * a ** 4 + 3392 * a ** 3
     + 26624 * a ** 2 + 80736 * a + 40032) + 184192 * a ** 2 + 24256 * a ** 3 + 1696 * a ** 4 + 64 * a ** 5
     + 304032))
result_5 = sp.simplify(substitute_m_r_t_into_diff_psi_squared - substitute_m_r_t_into_diff_psi_squared_in_paper)
if result_5 == 0:
    print("\nIs the formula of substitute_m_r_t_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_r_t_into_diff_psi_squared in paper correct?: error")


