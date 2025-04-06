# This script is designed to verify the correctness of the calculations in Cases 1 and 2 of Lemma 3.13 of the paper,
# allowing reviewers and readers to easily check and validate the results.
# The calculations are complex, involving judgments of polynomial parity and signs across multiple parameters.
# Including the full details in the paper would take up significant space and make verification more difficult.

import sympy as sp

# Define symbols
x = sp.symbols('x')
m, r, t, a = sp.symbols('m r t a', real=True)
m_lower_bound = 2*r + t + 7
r_lower_bound = 3
t_lower_bound = 2
a_lower_bound = 0

# Define matrix: matrix_lem_3_13_case_1
matrix_lem_3_13_case_1 = sp.Matrix([
    [0, 1, r, t],
    [1, 0, r, 0],
    [1, 1, 0, 0],
    [1, 0, 0, 0]
])
# The characteristic polynomial of matrix_lem_3_13_case_1
char_poly = matrix_lem_3_13_case_1.charpoly(x).as_expr()

# Define the expression for x
value_to_substitute = sp.sqrt(((1 + sp.sqrt(4*m - 7)) / 2)**2 + a - (m - 2*r - t - 1))

# The formular of \Phi_{M_4}(\sqrt{\left( \frac{1+\sqrt{4m-7}}{2}\right)^2+a-(m-2r-t-1)})
poly_substituted = char_poly.subs(x, value_to_substitute)

# Define psi_7 and psi_8 expressions
psi_7 = (a + r + t/2 - 1)*sp.sqrt(4*m - 7) + m - 2*a - 3*r - 3*t/2 + 2*a*r + a*t + r*t + a**2 - 1

psi_8 = 2*r*sp.sqrt(a + 2*r + t + (sp.sqrt(4*m - 7) - 1)/2)

#Compare
result_1 = sp.simplify(poly_substituted - (psi_7 - psi_8))

# Output
print(f'The characteristic polynomial of matrix_lem_3_13_case_1: {char_poly}')
if result_1 == 0:
    print("\nIs the formula of Phi_{M_3}(value_to_substitute) in paper correct?: correct")
else:
    print("\nIs the formula of Phi_{M_3}(value_to_substitute) in paper correct?: error")

"""The following verifies that the calculations in the paper are correct."""
# Compute psi_7^2(m,l,r,a) - psi_8^2(m,l,r,a)
diff_psi_squared = sp.expand(psi_7**2 - psi_8**2)
#The formula of psi_7^2 - psi_8^2 given in paper
diff_psi_squared_coefficients_1 = (((2*a + 2*r + t - 2)*m + 2*a**3 + 6*a**2*r + 3*a**2*t - 6*a**2 + 4*a*r**2 + 6*a*r*t - 14*a*r + a*t**2
    - 7*a*t + 2*a + 2*r**2*t - 8*r**2 + r*t**2 - 8*r*t + 4*r - 3/2*t**2 + 2*t + 2))

diff_psi_squared_coefficients_2 = (m**2 + (6*a**2 + 12*a*r + 6*a*t - 12*a + 4*r**2 + 6*r*t - 14*r + t**2 - 7*t + 2)*m
    + a**4 + 4*a**3*r + 2*a**3*t - 4*a**3 + 4*a**2*r**2 + 6*a**2*r*t - 14*a**2*r + a**2*t**2 - 7*a**2*t - 5*a**2
    + 4*a*r**2*t - 16*a*r**2 + 2*a*r*t**2 - 16*a*r*t - 6*a*r - 3*a*t**2 - 3*a*t + 18*a - 8*r**3 + r**2*t**2
    - 10*r**2*t + 4*r**2 - 3*r*t**2 + 20*r + 1/2*t**2 + 10*t - 6)

diff_psi_squared_in_paper = (diff_psi_squared_coefficients_1*sp.sqrt(4*m - 7) + diff_psi_squared_coefficients_2)
result_2 = sp.simplify(diff_psi_squared - diff_psi_squared_in_paper)
if result_2 == 0:
    print("\nIs the formula of psi_7^2 - psi_8^2 in paper correct?: correct")
else:
    print("\nIs the formula of psi_7^2 - psi_8^2 in paper correct?: error")
#The above result will show that diff_psi_squared = diff_psi_squared_in_paper
# - Both diff_psi_squared_coefficients_1 and diff_psi_squared_coefficients_2 can be regarded as functions of m.
# - We need to prove that both are monotonically increasing with respect to m.
# - diff_psi_squared_coefficients_1 is a linear function of m
# - diff_psi_squared_coefficients_2 is a quadratic function of m.
# - Taking the first derivative of diff_psi_squared_coefficients_2 with respect to m, the derivative need to be greater than 0.


#Easily, we know that the diff_psi_squared_coefficients_1 is monotonically increasing with respect to m when m >= 2*r + t + 7.
#It is easy to see that the value of the following polynomial is greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_1 = sp.collect(sp.simplify(diff_psi_squared_coefficients_1.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_1: {substitute_lb_into_diff_psi_squared_coefficients_1}')

diff_psi_squared_coefficients_2_derivative_m = sp.diff(diff_psi_squared_coefficients_2, m)
#From the following equation, we know that the diff_psi_squared_coefficients_2 is monotonically increasing with respect to m when m >= 2*r + t + 7.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m = sp.collect(sp.simplify(diff_psi_squared_coefficients_2_derivative_m.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m: {substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m}')

#----------------------------------------------------------------------------------------------------------------
#Substitute m_lower_bound = 2*r + t + 7 into psi_7^2 - psi_8^2.
substitute_m_into_diff_psi_squared = diff_psi_squared.subs(m,m_lower_bound)
#The formula of substitute_m_into_diff_psi_squared given in paper

substitute_m_into_coefficients_1 = ( ((a + r - 1/2)*t**2 + (3*a**2 + 6*a*r - 5*a + 2*r**2 - 4*r + 7)*t
    + 2*a**3 + 6*a**2*r - 6*a**2 + 4*a*r**2 - 10*a*r + 16*a - 4*r**2 + 14*r - 12))

substitute_m_into_coefficients_2 =  ((t**3 + (a**2 + 2*a*r + 3*a + r**2 + 5*r + 3/2)*t**2 + (2*a**3
    + 6*a**2*r - a**2 + 4*a*r**2 + 8*a*r + 27*a + 6*r**2 + 18*r - 23)*t + a**4 + 4*a**3*r - 4*a**3
    + 4*a**2*r**2 - 2*a**2*r + 37*a**2 + 8*a*r**2 + 54*a*r) - 66*a + 8*r**2 - 46*r + 57)

substitute_m_into_diff_psi_squared_in_paper = (substitute_m_into_coefficients_1*sp.sqrt(8*r + 4*t + 21) +substitute_m_into_coefficients_2)
# The above result will show that substitute_m_into_diff_psi_squared = substitute_m_into_diff_psi_squared_in_paper
# Both substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 can be regarded as functions of t.
# We need to show that both are monotonically increasing with respect to t.
# It is evident that substitute_m_into_coefficients_1 is a quadratic function of t,
# and substitute_m_into_coefficients_2 is a cubic functions of t.
# If the first derivative of substitute_m_into_coefficients_1 with respect to t,
# and the first and second derivatives of substitute_m_into_coefficients_2 are all greater than 0,
# then substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 will be monotonically increasing with respect to t.
result_3 = sp.simplify(substitute_m_into_diff_psi_squared - substitute_m_into_diff_psi_squared_in_paper)
if result_3 == 0:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: error")

substitute_m_into_coefficients_1_derivative_t_1 = sp.diff(substitute_m_into_coefficients_1, t)
#From the following equation, we know that the substitute_m_into_coefficients_1 is monotonically increasing with respect to t when t >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_1_derivative_t_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_t_1.subs(t, t_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_t_1: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_t_1}')

substitute_m_into_coefficients_2_derivative_t_1 = sp.diff(substitute_m_into_coefficients_2, t)
substitute_m_into_coefficients_2_derivative_t_2 = sp.diff(substitute_m_into_coefficients_2_derivative_t_1, t)
#From the following equation, we know that the substitute_m_into_coefficients_2 is monotonically increasing with respect to t when t >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_2_derivative_t_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_t_1.subs(t, t_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_2_derivative_t_2 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_t_2.subs(t, t_lower_bound)),a)
print(f'\nThe first derivative function of substitute_m_into_coefficients_2 with respect to t: {sp.collect(sp.expand(substitute_m_into_coefficients_1_derivative_t_1),t)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_t_1: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_t_1}')
print(f'The second derivative function of substitute_m_into_coefficients_2 with respect to t: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_t_2),t)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_t_2: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_t_2}')

#--------------------------------------------------------------------------------------------------------------------
#Substitute t_lower_bound = 2 into substitute_m_into_diff_psi_squared.
substitute_m_t_into_diff_psi_squared = substitute_m_into_diff_psi_squared.subs(t,t_lower_bound)
#The formula of substitute_m_t_into_diff_psi_squared given in paper
substitute_m_t_into_coefficients_1 = (4*a*r**2 + (6*a**2 + 2*a + 10)*r + 2*a**3 + 10*a)

substitute_m_t_into_coefficients_2 = (4*a**2 + 16*a + 24)*r**2 + (4*a**3 + 10*a**2 + 78*a + 10)*r  + a**4 + 39*a**2 + 25

substitute_m_t_into_diff_psi_squared_in_paper = (substitute_m_t_into_coefficients_1*sp.sqrt(8*r + 29)+substitute_m_t_into_coefficients_2)
result_4 = sp.simplify(substitute_m_t_into_diff_psi_squared - substitute_m_t_into_diff_psi_squared_in_paper)
if result_4 == 0:
    print("\nIs the formula of substitute_m_t_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_t_into_diff_psi_squared in paper correct?: error")
# The above result will show that substitute_m_t_into_diff_psi_squared = substitute_m_t_into_diff_psi_squared_in_paper
# Both substitute_m_t_into_coefficients_1 and substitute_m_t_into_coefficients_2 can be regarded as functions of r.
# We need to show that both are monotonically increasing with respect to r.
# It is evident that substitute_m_t_into_coefficients_1 and substitute_m_t_into_coefficients_2 are quadratic functions of r.
# If the first derivatives of these two polynomials with respect to r are all greater than 0,
# then substitute_m_t_into_coefficients_1 and substitute_m_t_into_coefficients_2 will be monotonically increasing with respect to r.
substitute_m_t_into_coefficients_1_derivative_r_1 = sp.diff(substitute_m_t_into_coefficients_1, r)
#From the following equation, we know that the substitute_m_t_into_coefficients_1 is monotonically increasing with respect to r when r >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_t_into_coefficients_1_derivative_r_1 = sp.collect(sp.simplify(substitute_m_t_into_coefficients_1_derivative_r_1.subs(r, r_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_t_into_coefficients_1_derivative_r_1: {substitute_lb_into_substitute_m_t_into_coefficients_1_derivative_r_1}')

substitute_m_t_into_coefficients_2_derivative_r_1 = sp.diff(substitute_m_t_into_coefficients_2, r)
#From the following equation, we know that the substitute_m_t_into_coefficients_2 is monotonically increasing with respect to r when r >= 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_r_1 = sp.collect(sp.simplify(substitute_m_t_into_coefficients_2_derivative_r_1.subs(r, r_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_r_1: {substitute_lb_into_substitute_m_t_into_coefficients_2_derivative_r_1}')

#--------------------------------------------------------------------------------------------------------------
#Substitute r_lower_bound = 3 into substitute_m_t_into_diff_psi_squared.
substitute_m_t_r_into_diff_psi_squared = substitute_m_t_into_diff_psi_squared.subs(r,r_lower_bound)
#The formula of substitute_m_t_r_into_diff_psi_squared given in paper
substitute_m_t_r_into_diff_psi_squared_in_paper =(
        378*a + sp.sqrt(53)*(2*a**3 + 18*a**2 + 52*a + 30)
        + 105*a**2 + 12*a**3 + a**4 + 271
)
result_5 = sp.simplify(substitute_m_t_r_into_diff_psi_squared - substitute_m_t_r_into_diff_psi_squared_in_paper)
if result_5 == 0:
    print("\nIs the formula of substitute_m_t_r_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_t_r_into_diff_psi_squared in paper correct?: error")