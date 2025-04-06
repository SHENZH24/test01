# This script is designed to verify the correctness of the calculations in Cases 3 of Lemma 3.13 of the paper,
# allowing reviewers and readers to easily check and validate the results.
# The calculations are complex, involving judgments of polynomial parity and signs across multiple parameters.
# Including the full details in the paper would take up significant space and make verification more difficult.

import sympy as sp

# Define symbols
x = sp.symbols('x')
m, k_1, t_1, a = sp.symbols('m k_1 t_1 a', real=True)
m_lower_bound = 2*t_1 + 2*k_1 + 3
k_1_lower_bound = 2
t_1_lower_bound = 1
a_lower_bound = 0

# Define matrix: matrix_lem_3_13_case_3
matrix_lem_3_13_case_3 = sp.Matrix([
    [0, 2, 0, t_1],
    [2, 1, k_1-1, 0],
    [0, 2, 0, 0],
    [2, 0, 0, 0]
])
# The characteristic polynomial of matrix_lem_3_13_case_3
char_poly = matrix_lem_3_13_case_3.charpoly(x).as_expr()

# Define the expression for x
value_to_substitute = sp.sqrt(((1 + sp.sqrt(4*m - 7)) / 2)**2 + a - (m - 2*t_1 - 2*k_1 - 3))

# The formular of \Phi_{M_5}\left( \sqrt{\left( \frac{1+\sqrt{4m-7}}{2}\right)^2+a-(m-2t_1-2k_1-3)}\right)
poly_substituted = char_poly.subs(x, value_to_substitute)

# Define psi_9 and psi_10 expressions
psi_9 = (
    (2*a + 2*k_1 + 2*t_1 + 1)*sp.sqrt(4*m - 7)
    + 2*a - 2*k_1 + 2*m - 10*t_1 + 4*a*k_1 + 4*a*t_1 + 8*k_1*t_1 + 2*a**2 - 5
)

psi_10 = (
    sp.sqrt(2)*(a + 2*k_1 + (3 + sp.sqrt(4*m - 7))/2)
    * sp.sqrt(2*a + 4*k_1 + 4*t_1 + sp.sqrt(4*m - 7) + 3)
)


#Compare
result_1 = sp.simplify(poly_substituted - 1/2 * (psi_9 - psi_10))

# Output
print(f'The characteristic polynomial of matrix_lem_3_11_case_3: {char_poly}')
if result_1 == 0:
    print("\nIs the formula of Phi_{M_4}(value_to_substitute) in paper correct?: correct")
else:
    print("\nIs the formula of Phi_{M_4}(value_to_substitute) in paper correct?: error")

"""The following verifies that the calculations in the paper are correct."""
# Compute psi_9^2(m,l,r,a) - psi_10^2(m,l,r,a)
diff_psi_squared = sp.expand(psi_9**2 - psi_10**2)
#The formula of psi_9^2 - psi_10^2 given in paper

diff_psi_squared_coefficients_1 = (((8*a + 8*k_1 + 8*t_1 + 2)*m
     + 8*a**3 + 24*a**2*k_1 + 24*a**2*t_1 + 6*a**2
     + 16*a*k_1**2 + 64*a*k_1*t_1 - 16*a*k_1
     + 16*a*t_1**2 - 32*a*t_1 - 34*a
     + 32*k_1**2*t_1 - 32*k_1**2 + 32*k_1*t_1**2
     - 48*k_1*t_1 - 60*k_1 - 40*t_1**2 - 52*t_1 - 20))

diff_psi_squared_coefficients_2 = (4*m**2
    + (24*a**2 + 48*a*k_1 + 48*a*t_1 + 12*a
       + 16*k_1**2 + 64*k_1*t_1 - 16*k_1 + 16*t_1**2 - 32*t_1 - 34)*m
    + 4*a**4 + 16*a**3*k_1 + 16*a**3*t_1 + 4*a**3
    + 16*a**2*k_1**2 + 64*a**2*k_1*t_1 - 16*a**2*k_1
    + 16*a**2*t_1**2 - 32*a**2*t_1 - 62*a**2
    + 64*a*k_1**2*t_1 - 64*a*k_1**2 + 64*a*k_1*t_1**2
    - 96*a*k_1*t_1 - 176*a*k_1 - 80*a*t_1**2 - 160*a*t_1 - 54*a
    - 32*k_1**3 + 64*k_1**2*t_1**2 - 64*k_1**2*t_1
    - 96*k_1**2 - 160*k_1*t_1**2 - 144*k_1*t_1
    - 20*k_1 + 72*t_1**2 + 68*t_1 + 36)


diff_psi_squared_in_paper = (diff_psi_squared_coefficients_1*sp.sqrt(4*m - 7) + diff_psi_squared_coefficients_2)
result_2 = sp.simplify(diff_psi_squared - diff_psi_squared_in_paper)
if result_2 == 0:
    print("\nIs the formula of psi_9^2 - psi_10^2 in paper correct?: correct")
else:
    print("\nIs the formula of psi_9^2 - psi_10^2 in paper correct?: error")
#The above result will show that diff_psi_squared = diff_psi_squared_in_paper
# - Both diff_psi_squared_coefficients_1 and diff_psi_squared_coefficients_2 can be regarded as functions of m.
# - We need to prove that both are monotonically increasing with respect to m.
# - diff_psi_squared_coefficients_1 is a linear function of m
# - diff_psi_squared_coefficients_2 is a quadratic function of m.
# - Taking the first derivative of diff_psi_squared_coefficients_2 with respect to m, the derivative need to be greater than 0.

#Easily, we know that the diff_psi_squared_coefficients_1 is monotonically increasing with respect to m when m >= 2*t_1 + 2*k_1 + 3.
#It is easy to see that the value of the following polynomial is greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_1 = sp.collect(sp.simplify(diff_psi_squared_coefficients_1.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_1: {substitute_lb_into_diff_psi_squared_coefficients_1}')

diff_psi_squared_coefficients_2_derivative_m = sp.diff(diff_psi_squared_coefficients_2, m)
#From the following equation, we know that the diff_psi_squared_coefficients_2 is monotonically increasing with respect to m when m >= 2*t_1 + 2*k_1 + 3.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m = sp.collect(sp.simplify(diff_psi_squared_coefficients_2_derivative_m.subs(m, m_lower_bound)),a)
print(f'The value of substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m: {substitute_lb_into_diff_psi_squared_coefficients_2_derivative_m}')

#----------------------------------------------------------------------------------------------------------------
#Substitute m_lower_bound = 2*t_1 + 2*k_1 + 3 into psi_9^2 - psi_10^2.
substitute_m_into_diff_psi_squared = diff_psi_squared.subs(m,m_lower_bound)
#The formula of substitute_m_into_diff_psi_squared given in paper
substitute_m_into_coefficients_1 =  ((16*a + 32*k_1 - 24)*t_1**2
     + (24*a**2 + 64*a*k_1 - 16*a + 32*k_1**2 - 16*k_1 - 24)*t_1
     + 8*a**3 + 24*a**2*k_1 + 6*a**2 + 16*a*k_1**2 - 10*a
     - 16*k_1**2 - 32*k_1 - 14)

substitute_m_into_coefficients_2  = (32*t_1**3
    + (16*a**2 + 64*a*k_1 + 16*a + 64*k_1**2 + 72)*t_1**2
    + (16*a**3 + 64*a**2*k_1 + 16*a**2 + 64*a*k_1**2
       + 96*a*k_1 + 8*a + 96*k_1**2 - 16*k_1 - 48)*t_1
    + 4*a**4 + 16*a**3*k_1 + 4*a**3
    + 16*a**2*k_1**2 + 32*a**2*k_1 + 10*a**2
    + 32*a*k_1**2 - 8*a*k_1 - 18*a - 64*k_1**2 - 88*k_1 - 30)

substitute_m_into_diff_psi_squared_in_paper = (substitute_m_into_coefficients_1*sp.sqrt(8*k_1 + 8*t_1 + 5) +substitute_m_into_coefficients_2)
# The above result will show that substitute_m_into_diff_psi_squared = substitute_m_into_diff_psi_squared_in_paper
# Both substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 can be regarded as functions of t_1.
# We need to show that both are monotonically increasing with respect to t_1.
# It is evident that substitute_m_into_coefficients_1 is a quadratic function of t_1,
# and substitute_m_into_coefficients_2 is a cubic functions of t_1.
# If the first derivative of substitute_m_into_coefficients_1 with respect to t_1,
# and the first and second derivatives of substitute_m_into_coefficients_2 are all greater than 0,
# then substitute_m_into_coefficients_1 and substitute_m_into_coefficients_2 will be monotonically increasing with respect to t_1.
result_3 = sp.simplify(substitute_m_into_diff_psi_squared - substitute_m_into_diff_psi_squared_in_paper)
if result_3 == 0:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_into_diff_psi_squared in paper correct?: error")

substitute_m_into_coefficients_1_derivative_t1_1 = sp.diff(substitute_m_into_coefficients_1, t_1)
#From the following equation, we know that the substitute_m_into_coefficients_1 is monotonically increasing with respect to t_1 when t_1 >= 1.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_1_derivative_t1_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_1_derivative_t1_1.subs(t_1, t_1_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_1_derivative_t1_1: {substitute_lb_into_substitute_m_into_coefficients_1_derivative_t1_1}')

substitute_m_into_coefficients_2_derivative_t1_1 = sp.diff(substitute_m_into_coefficients_2, t_1)
substitute_m_into_coefficients_2_derivative_t1_2 = sp.diff(substitute_m_into_coefficients_2_derivative_t1_1, t_1)
#From the following equation, we know that the substitute_m_into_coefficients_2 is monotonically increasing with respect to t_1 when t_1 >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_into_coefficients_2_derivative_t1_1 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_t1_1.subs(t_1, t_1_lower_bound)),a)
substitute_lb_into_substitute_m_into_coefficients_2_derivative_t1_2 = sp.collect(sp.simplify(substitute_m_into_coefficients_2_derivative_t1_2.subs(t_1, t_1_lower_bound)),a)
print(f'\nThe first derivative function of substitute_m_into_coefficients_2 with respect to t_1: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_t1_1),t_1)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_t1_1: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_t1_1}')
print(f'The second derivative function of substitute_m_into_coefficients_2 with respect to t_1: {sp.collect(sp.expand(substitute_m_into_coefficients_2_derivative_t1_2),t_1)}')
print(f'The value of substitute_lb_into_substitute_m_into_coefficients_2_derivative_t1_2: {substitute_lb_into_substitute_m_into_coefficients_2_derivative_t1_2}')

#--------------------------------------------------------------------------------------------------------------------
#Substitute t_1_lower_bound = 1 into substitute_m_into_diff_psi_squared.
substitute_m_t1_into_diff_psi_squared = substitute_m_into_diff_psi_squared.subs(t_1,t_1_lower_bound)
#The formula of substitute_m_t1_into_diff_psi_squared given in paper
substitute_m_t1_into_coefficients_1 = (((16*a + 16)*k_1**2 
     + (24*a**2 + 64*a - 16)*k_1 
     + 8*a**3 + 30*a**2 - 10*a - 62) )

substitute_m_t1_into_coefficients_2 = ((16*a**2 + 96*a + 96)*k_1**2 
    + (16*a**3 + 96*a**2 + 152*a - 104)*k_1 
    + 4*a**4 + 20*a**3 + 42*a**2 + 6*a + 26)

substitute_m_t1_into_diff_psi_squared_in_paper = (substitute_m_t1_into_coefficients_1*sp.sqrt(8*k_1 + 13)+substitute_m_t1_into_coefficients_2)
result_4 = sp.simplify(substitute_m_t1_into_diff_psi_squared - substitute_m_t1_into_diff_psi_squared_in_paper)
if result_4 == 0:
    print("\nIs the formula of substitute_m_t1_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_t1_into_diff_psi_squared in paper correct?: error")
# The above result will show that substitute_m_t1_into_diff_psi_squared = substitute_m_t1_into_diff_psi_squared_in_paper
# Both substitute_m_t1_into_coefficients_1 and substitute_m_t1_into_coefficients_2 can be regarded as functions of k_1.
# We need to show that both are monotonically increasing with respect to k_1.
# It is evident that substitute_m_t1_into_coefficients_1 and substitute_m_t1_into_coefficients_2 are quadratic functions of k_1.
# If the first derivatives of these two polynomials with respect to k_1 are all greater than 0,
# then substitute_m_t_into_coefficients_1 and substitute_m_t_into_coefficients_2 will be monotonically increasing with respect to k_1.
substitute_m_t1_into_coefficients_1_derivative_k1_1 = sp.diff(substitute_m_t1_into_coefficients_1, k_1)
#From the following equation, we know that the substitute_m_t1_into_coefficients_1 is monotonically increasing with respect to k_1 when k_1 >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_t1_into_coefficients_1_derivative_k1_1 = sp.collect(sp.simplify(substitute_m_t1_into_coefficients_1_derivative_k1_1.subs(k_1, k_1_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_t1_into_coefficients_1_derivative_k1_1: {substitute_lb_into_substitute_m_t1_into_coefficients_1_derivative_k1_1}')

substitute_m_t1_into_coefficients_2_derivative_k1_1 = sp.diff(substitute_m_t1_into_coefficients_2, k_1)
#From the following equation, we know that the substitute_m_t1_into_coefficients_2 is monotonically increasing with respect to k_1 when k_1 >= 2.
#It is easy to see that the value of the following polynomials are greater than 0.
substitute_lb_into_substitute_m_t1_into_coefficients_2_derivative_k1_1 = sp.collect(sp.simplify(substitute_m_t1_into_coefficients_2_derivative_k1_1.subs(k_1, k_1_lower_bound)),a)
print(f'The value of substitute_lb_into_substitute_m_t1_into_coefficients_2_derivative_k1_1: {substitute_lb_into_substitute_m_t1_into_coefficients_2_derivative_k1_1}')

#--------------------------------------------------------------------------------------------------------------
#Substitute k_1_lower_bound = 2 into substitute_m_t1_into_diff_psi_squared.
substitute_m_t1_k1_into_diff_psi_squared = substitute_m_t1_into_diff_psi_squared.subs(k_1,k_1_lower_bound)
#The formula of substitute_m_t1_k1_into_diff_psi_squared given in paper
substitute_m_t1_k1_into_diff_psi_squared_in_paper =(
    694*a + 182*sp.sqrt(29)*a - 30*sp.sqrt(29) + 78*sp.sqrt(29)*a**2 + 8*sp.sqrt(29)*a**3
    + 298*a**2 + 52*a**3 + 4*a**4 + 202)
result_5 = sp.simplify(substitute_m_t1_k1_into_diff_psi_squared - substitute_m_t1_k1_into_diff_psi_squared_in_paper)
if result_5 == 0:
    print("\nIs the formula of substitute_m_t1_k1_into_diff_psi_squared in paper correct?: correct")
else:
    print("\nIs the formula of substitute_m_t1_k1_into_diff_psi_squared in paper correct?: error")