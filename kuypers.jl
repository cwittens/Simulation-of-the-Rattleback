# equations for omega taken from wackelstein_gleichungen_omega.ipynb
# wi_dot_mathematica is from mathematica rattleback_omega.nb
# The algebraic optimized equations made the code more than 20% faster

#using @inline reduces number of allocations by 3 for EACH rhs call
@inline function w1_dot(alpha, beta, omega1, omega2, omega3, x1, x2, x3,
        k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
        return (A * B * h * m * omega1 * omega2 * x1 + A * B * m * omega1 * omega2 * x1 * x3 -
                A * C * F * omega1 * omega3 - A * C * m * omega1 * omega3 * x1 * x2 + A * F * h * m * omega1 * omega2 * x2 +
                A * F * m * omega1 * omega2 * x2 * x3 - A * F * m * omega1 * omega3 * x1^2 - A * F * m * omega1 * omega3 * x2^2 +
                A * h^3 * m^2 * omega1 * omega2 * x1 + 3 * A * h^2 * m^2 * omega1 * omega2 * x1 * x3 -
                A * h^2 * m^2 * omega1 * omega3 * x1 * x2 + A * h * m^2 * omega1 * omega2 * x1^3 +
                A * h * m^2 * omega1 * omega2 * x1 * x2^2 + 3 * A * h * m^2 * omega1 * omega2 * x1 * x3^2 -
                2 * A * h * m^2 * omega1 * omega3 * x1 * x2 * x3 + A * m^2 * omega1 * omega2 * x1^3 * x3 +
                A * m^2 * omega1 * omega2 * x1 * x2^2 * x3 + A * m^2 * omega1 * omega2 * x1 * x3^3 -
                A * m^2 * omega1 * omega3 * x1^3 * x2 - A * m^2 * omega1 * omega3 * x1 * x2^3 -
                A * m^2 * omega1 * omega3 * x1 * x2 * x3^2 + B^2 * C * omega2 * omega3 - B^2 * h * m * omega1 * omega2 * x1 -
                B^2 * m * omega1 * omega2 * x1 * x3 + B^2 * m * omega2 * omega3 * x1^2 + B^2 * m * omega2 * omega3 * x2^2 -
                B * C^2 * omega2 * omega3 - B * C * F * omega1 * omega3 - B * C * c_1 * omega1 -
                B * C * g * h * m * sin(alpha) + B * C * g * m * x2 * cos(alpha) * cos(beta) -
                B * C * g * m * x3 * sin(alpha) + B * C * h^2 * m * omega2 * omega3 - B * C * h * k2 * m +
                2 * B * C * h * m * omega2 * omega3 * x3 - B * C * k2 * m * x3 + B * C * k3 * m * x2 - B * C * m * omega2 * omega3 * x2^2 +
                B * C * m * omega2 * omega3 * x3^2 + B * F * h * m * omega1^2 * x1 - B * F * h * m * omega1 * omega2 * x2 -
                B * F * h * m * omega2^2 * x1 + B * F * m * omega1^2 * x1 * x3 - B * F * m * omega1 * omega2 * x2 * x3 -
                B * F * m * omega1 * omega3 * x1^2 - B * F * m * omega1 * omega3 * x2^2 - B * F * m * omega2^2 * x1 * x3 -
                B * c_1 * m * omega1 * x1^2 - B * c_1 * m * omega1 * x2^2 - B * c_3 * h * m * omega3 * x1 -
                B * c_3 * m * omega3 * x1 * x3 + B * g * h * m^2 * x1 * x2 * sin(beta) * cos(alpha) -
                B * g * h * m^2 * x2^2 * sin(alpha) + B * g * m^2 * x1^2 * x2 * cos(alpha) * cos(beta) +
                B * g * m^2 * x1 * x2 * x3 * sin(beta) * cos(alpha) + B * g * m^2 * x2^3 * cos(alpha) * cos(beta) -
                B * g * m^2 * x2^2 * x3 * sin(alpha) - B * h^3 * m^2 * omega1 * omega2 * x1 -
                3 * B * h^2 * m^2 * omega1 * omega2 * x1 * x3 + B * h^2 * m^2 * omega2 * omega3 * x1^2 -
                B * h * k1 * m^2 * x1 * x2 - B * h * k2 * m^2 * x2^2 - B * h * m^2 * omega1 * omega2 * x1^3 -
                B * h * m^2 * omega1 * omega2 * x1 * x2^2 - 3 * B * h * m^2 * omega1 * omega2 * x1 * x3^2 +
                2 * B * h * m^2 * omega2 * omega3 * x1^2 * x3 - B * k1 * m^2 * x1 * x2 * x3 - B * k2 * m^2 * x2^2 * x3 +
                B * k3 * m^2 * x1^2 * x2 + B * k3 * m^2 * x2^3 - B * m^2 * omega1 * omega2 * x1^3 * x3 -
                B * m^2 * omega1 * omega2 * x1 * x2^2 * x3 - B * m^2 * omega1 * omega2 * x1 * x3^3 +
                B * m^2 * omega2 * omega3 * x1^4 + B * m^2 * omega2 * omega3 * x1^2 * x2^2 +
                B * m^2 * omega2 * omega3 * x1^2 * x3^2 + C^2 * F * omega1 * omega3 - C^2 * h^2 * m * omega2 * omega3 -
                2 * C^2 * h * m * omega2 * omega3 * x3 + C^2 * m * omega1 * omega3 * x1 * x2 - C^2 * m * omega2 * omega3 * x1^2 -
                C^2 * m * omega2 * omega3 * x3^2 + C * F^2 * omega2 * omega3 - C * F * c_2 * omega2 -
                C * F * g * h * m * sin(beta) * cos(alpha) - C * F * g * m * x1 * cos(alpha) * cos(beta) -
                C * F * g * m * x3 * sin(beta) * cos(alpha) - C * F * h^2 * m * omega1 * omega3 + C * F * h * k1 * m -
                2 * C * F * h * m * omega1 * omega3 * x3 + C * F * k1 * m * x3 - C * F * k3 * m * x1 + C * F * m * omega1 * omega3 * x2^2 -
                C * F * m * omega1 * omega3 * x3^2 + C * F * m * omega2 * omega3 * x1 * x2 - C * c_1 * h^2 * m * omega1 -
                2 * C * c_1 * h * m * omega1 * x3 - C * c_1 * m * omega1 * x1^2 - C * c_1 * m * omega1 * x3^2 -
                C * c_2 * m * omega2 * x1 * x2 - C * g * h^3 * m^2 * sin(alpha) + C * g * h^2 * m^2 * x2 * cos(alpha) * cos(beta) -
                3 * C * g * h^2 * m^2 * x3 * sin(alpha) - C * g * h * m^2 * x1^2 * sin(alpha) -
                C * g * h * m^2 * x1 * x2 * sin(beta) * cos(alpha) + 2 * C * g * h * m^2 * x2 * x3 * cos(alpha) * cos(beta) -
                3 * C * g * h * m^2 * x3^2 * sin(alpha) - C * g * m^2 * x1^2 * x3 * sin(alpha) -
                C * g * m^2 * x1 * x2 * x3 * sin(beta) * cos(alpha) + C * g * m^2 * x2 * x3^2 * cos(alpha) * cos(beta) -
                C * g * m^2 * x3^3 * sin(alpha) - C * h^3 * k2 * m^2 - 3 * C * h^2 * k2 * m^2 * x3 + C * h^2 * k3 * m^2 * x2 +
                C * h^2 * m^2 * omega1 * omega3 * x1 * x2 - C * h^2 * m^2 * omega2 * omega3 * x1^2 + C * h * k1 * m^2 * x1 * x2 -
                C * h * k2 * m^2 * x1^2 - 3 * C * h * k2 * m^2 * x3^2 + 2 * C * h * k3 * m^2 * x2 * x3 +
                2 * C * h * m^2 * omega1 * omega3 * x1 * x2 * x3 - 2 * C * h * m^2 * omega2 * omega3 * x1^2 * x3 +
                C * k1 * m^2 * x1 * x2 * x3 - C * k2 * m^2 * x1^2 * x3 - C * k2 * m^2 * x3^3 + C * k3 * m^2 * x2 * x3^2 +
                C * m^2 * omega1 * omega3 * x1^3 * x2 + C * m^2 * omega1 * omega3 * x1 * x2^3 +
                C * m^2 * omega1 * omega3 * x1 * x2 * x3^2 - C * m^2 * omega2 * omega3 * x1^4 -
                C * m^2 * omega2 * omega3 * x1^2 * x2^2 - C * m^2 * omega2 * omega3 * x1^2 * x3^2 +
                F^2 * h * m * omega1^2 * x2 - F^2 * h * m * omega2^2 * x2 + F^2 * m * omega1^2 * x2 * x3 -
                F^2 * m * omega2^2 * x2 * x3 + F^2 * m * omega2 * omega3 * x1^2 + F^2 * m * omega2 * omega3 * x2^2 -
                F * c_2 * m * omega2 * x1^2 - F * c_2 * m * omega2 * x2^2 - F * c_3 * h * m * omega3 * x2 -
                F * c_3 * m * omega3 * x2 * x3 - F * g * h * m^2 * x1^2 * sin(beta) * cos(alpha) +
                F * g * h * m^2 * x1 * x2 * sin(alpha) - F * g * m^2 * x1^3 * cos(alpha) * cos(beta) -
                F * g * m^2 * x1^2 * x3 * sin(beta) * cos(alpha) - F * g * m^2 * x1 * x2^2 * cos(alpha) * cos(beta) +
                F * g * m^2 * x1 * x2 * x3 * sin(alpha) + F * h^3 * m^2 * omega1^2 * x1 - F * h^3 * m^2 * omega2^2 * x1 +
                3 * F * h^2 * m^2 * omega1^2 * x1 * x3 - F * h^2 * m^2 * omega1 * omega3 * x1^2 -
                3 * F * h^2 * m^2 * omega2^2 * x1 * x3 + F * h^2 * m^2 * omega2 * omega3 * x1 * x2 + F * h * k1 * m^2 * x1^2 +
                F * h * k2 * m^2 * x1 * x2 + F * h * m^2 * omega1^2 * x1^3 + F * h * m^2 * omega1^2 * x1 * x2^2 +
                3 * F * h * m^2 * omega1^2 * x1 * x3^2 - 2 * F * h * m^2 * omega1 * omega3 * x1^2 * x3 -
                F * h * m^2 * omega2^2 * x1^3 - F * h * m^2 * omega2^2 * x1 * x2^2 - 3 * F * h * m^2 * omega2^2 * x1 * x3^2 +
                2 * F * h * m^2 * omega2 * omega3 * x1 * x2 * x3 + F * k1 * m^2 * x1^2 * x3 + F * k2 * m^2 * x1 * x2 * x3 -
                F * k3 * m^2 * x1^3 - F * k3 * m^2 * x1 * x2^2 + F * m^2 * omega1^2 * x1^3 * x3 +
                F * m^2 * omega1^2 * x1 * x2^2 * x3 + F * m^2 * omega1^2 * x1 * x3^3 - F * m^2 * omega1 * omega3 * x1^4 -
                F * m^2 * omega1 * omega3 * x1^2 * x2^2 - F * m^2 * omega1 * omega3 * x1^2 * x3^2 -
                F * m^2 * omega2^2 * x1^3 * x3 - F * m^2 * omega2^2 * x1 * x2^2 * x3 - F * m^2 * omega2^2 * x1 * x3^3 +
                F * m^2 * omega2 * omega3 * x1^3 * x2 + F * m^2 * omega2 * omega3 * x1 * x2^3 +
                F * m^2 * omega2 * omega3 * x1 * x2 * x3^2 - c_1 * h^2 * m^2 * omega1 * x1^2 -
                2 * c_1 * h * m^2 * omega1 * x1^2 * x3 - c_1 * m^2 * omega1 * x1^4 - c_1 * m^2 * omega1 * x1^2 * x2^2 -
                c_1 * m^2 * omega1 * x1^2 * x3^2 - c_2 * h^2 * m^2 * omega2 * x1 * x2 - 2 * c_2 * h * m^2 * omega2 * x1 * x2 * x3 -
                c_2 * m^2 * omega2 * x1^3 * x2 - c_2 * m^2 * omega2 * x1 * x2^3 - c_2 * m^2 * omega2 * x1 * x2 * x3^2 -
                c_3 * h^3 * m^2 * omega3 * x1 - 3 * c_3 * h^2 * m^2 * omega3 * x1 * x3 - c_3 * h * m^2 * omega3 * x1^3 -
                c_3 * h * m^2 * omega3 * x1 * x2^2 - 3 * c_3 * h * m^2 * omega3 * x1 * x3^2 - c_3 * m^2 * omega3 * x1^3 * x3 -
                c_3 * m^2 * omega3 * x1 * x2^2 * x3 - c_3 * m^2 * omega3 * x1 * x3^3) / (A * B * C + A * B * m * x1^2 +
                A * B * m * x2^2 + A * C * h^2 * m + 2 * A * C * h * m * x3 + A * C * m * x1^2 + A * C * m * x3^2 +
                A * h^2 * m^2 * x1^2 + 2 * A * h * m^2 * x1^2 * x3 + A * m^2 * x1^4 + A * m^2 * x1^2 * x2^2 +
                A * m^2 * x1^2 * x3^2 + B * C * h^2 * m + 2 * B * C * h * m * x3 + B * C * m * x2^2 + B * C * m * x3^2 +
                B * h^2 * m^2 * x2^2 + 2 * B * h * m^2 * x2^2 * x3 + B * m^2 * x1^2 * x2^2 + B * m^2 * x2^4 +
                B * m^2 * x2^2 * x3^2 - C * F^2 - 2 * C * F * m * x1 * x2 + C * h^4 * m^2 + 4 * C * h^3 * m^2 * x3 +
                C * h^2 * m^2 * x1^2 + C * h^2 * m^2 * x2^2 + 6 * C * h^2 * m^2 * x3^2 + 2 * C * h * m^2 * x1^2 * x3 +
                2 * C * h * m^2 * x2^2 * x3 + 4 * C * h * m^2 * x3^3 + C * m^2 * x1^2 * x3^2 + C * m^2 * x2^2 * x3^2 +
                C * m^2 * x3^4 - F^2 * m * x1^2 - F^2 * m * x2^2 - 2 * F * h^2 * m^2 * x1 * x2 -
                4 * F * h * m^2 * x1 * x2 * x3 - 2 * F * m^2 * x1^3 * x2 - 2 * F * m^2 * x1 * x2^3 - 2 * F * m^2 * x1 * x2 * x3^2)
    end

# equations for omega_mathematica taken from wackelstein_gleichungen_omega.nb 
# this function is way shorter, but for some inital conditions it didnt seems stable.
# this was not further investigated, as it was only 1.2x faster
#=
function w1_dot_mathematica(alpha, beta, omega1, omega2, omega3, x1, x2, x3,
    k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
    return (-((-((-(F * x2) - x1 * (B + m * (x1^2 + x2^2 + (h + x3)^2))) *
                (-((C + m * (x1^2 + x2^2)) * (c_1 * omega1 + F * omega1 * omega3 + (-B + C) * omega2 * omega3 - k3 * m * x2 + k2 * m * (h + x3) -
                                            g * m * x2 * cos(alpha) * cos(beta) + g * m * (h + x3) * sin(alpha))) +
                m * x1 * (h + x3) * ((A - B) * omega1 * omega2 + F * (omega1 - omega2) * (omega1 + omega2) - c_3 * omega3 + k2 * m * x1 - k1 * m * x2 +
                                    g * m * x1 * sin(alpha) + g * m * x2 * cos(alpha) * sin(beta)))) +
            ((F + m * x1 * x2) * (C + m * (x1^2 + x2^2)) + m^2 * x1 * x2 * (h + x3)^2) *
            (x2 * (c_1 * omega1 + F * omega1 * omega3 + (-B + C) * omega2 * omega3 - k3 * m * x2 + k2 * m * (h + x3) - g * m * x2 * cos(alpha) * cos(beta) +
                    g * m * (h + x3) * sin(alpha)) + x1 * (-(c_2 * omega2) + (-A + C) * omega1 * omega3 + F * omega2 * omega3 - k3 * m * x1 + k1 * m * (h + x3) -
                                                            g * m * cos(alpha) * (x1 * cos(beta) + (h + x3) * sin(beta))))) /
            (((F + m * x1 * x2) * (C + m * (x1^2 + x2^2)) + m^2 * x1 * x2 * (h + x3)^2) * (F * x1 + x2 * (A + m * (x1^2 + x2^2 + (h + x3)^2))) -
            (m^2 * x1^2 * (h + x3)^2 - (C + m * (x1^2 + x2^2)) * (A + m * (x2^2 + (h + x3)^2))) *
            (-(F * x2) - x1 * (B + m * (x1^2 + x2^2 + (h + x3)^2))))))
end
=#

@inline  function w2_dot(alpha, beta, omega1, omega2, omega3, x1, x2, x3,
        k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
        return (-A^2 * C * omega1 * omega3 + A^2 * h * m * omega1 * omega2 * x2 +
                A^2 * m * omega1 * omega2 * x2 * x3 - A^2 * m * omega1 * omega3 * x1^2 -
                A^2 * m * omega1 * omega3 * x2^2 - A * B * h * m * omega1 * omega2 * x2 -
                A * B * m * omega1 * omega2 * x2 * x3 + A * C^2 * omega1 * omega3 +
                A * C * F * omega2 * omega3 - A * C * c_2 * omega2 - A * C * g * h * m * sin(beta) * cos(alpha) -
                A * C * g * m * x1 * cos(alpha) * cos(beta) - A * C * g * m * x3 * sin(beta) * cos(alpha) -
                A * C * h^2 * m * omega1 * omega3 + A * C * h * k1 * m - 2 * A * C * h * m * omega1 * omega3 * x3 +
                A * C * k1 * m * x3 - A * C * k3 * m * x1 + A * C * m * omega1 * omega3 * x1^2 -
                A * C * m * omega1 * omega3 * x3^2 + A * F * h * m * omega1^2 * x2 + A * F * h * m * omega1 * omega2 * x1 -
                A * F * h * m * omega2^2 * x2 + A * F * m * omega1^2 * x2 * x3 + A * F * m * omega1 * omega2 * x1 * x3 -
                A * F * m * omega2^2 * x2 * x3 + A * F * m * omega2 * omega3 * x1^2 + A * F * m * omega2 * omega3 * x2^2 -
                A * c_2 * m * omega2 * x1^2 - A * c_2 * m * omega2 * x2^2 - A * c_3 * h * m * omega3 * x2 -
                A * c_3 * m * omega3 * x2 * x3 - A * g * h * m^2 * x1^2 * sin(beta) * cos(alpha) +
                A * g * h * m^2 * x1 * x2 * sin(alpha) - A * g * m^2 * x1^3 * cos(alpha) * cos(beta) -
                A * g * m^2 * x1^2 * x3 * sin(beta) * cos(alpha) - A * g * m^2 * x1 * x2^2 * cos(alpha) * cos(beta) +
                A * g * m^2 * x1 * x2 * x3 * sin(alpha) + A * h^3 * m^2 * omega1 * omega2 * x2 +
                3 * A * h^2 * m^2 * omega1 * omega2 * x2 * x3 - A * h^2 * m^2 * omega1 * omega3 * x2^2 +
                A * h * k1 * m^2 * x1^2 + A * h * k2 * m^2 * x1 * x2 + A * h * m^2 * omega1 * omega2 * x1^2 * x2 +
                A * h * m^2 * omega1 * omega2 * x2^3 + 3 * A * h * m^2 * omega1 * omega2 * x2 * x3^2 -
                2 * A * h * m^2 * omega1 * omega3 * x2^2 * x3 + A * k1 * m^2 * x1^2 * x3 + A * k2 * m^2 * x1 * x2 * x3 -
                A * k3 * m^2 * x1^3 - A * k3 * m^2 * x1 * x2^2 + A * m^2 * omega1 * omega2 * x1^2 * x2 * x3 +
                A * m^2 * omega1 * omega2 * x2^3 * x3 + A * m^2 * omega1 * omega2 * x2 * x3^3 -
                A * m^2 * omega1 * omega3 * x1^2 * x2^2 - A * m^2 * omega1 * omega3 * x2^4 -
                A * m^2 * omega1 * omega3 * x2^2 * x3^2 + B * C * F * omega2 * omega3 +
                B * C * m * omega2 * omega3 * x1 * x2 - B * F * h * m * omega1 * omega2 * x1 -
                B * F * m * omega1 * omega2 * x1 * x3 + B * F * m * omega2 * omega3 * x1^2 +
                B * F * m * omega2 * omega3 * x2^2 - B * h^3 * m^2 * omega1 * omega2 * x2 -
                3 * B * h^2 * m^2 * omega1 * omega2 * x2 * x3 + B * h^2 * m^2 * omega2 * omega3 * x1 * x2 -
                B * h * m^2 * omega1 * omega2 * x1^2 * x2 - B * h * m^2 * omega1 * omega2 * x2^3 -
                3 * B * h * m^2 * omega1 * omega2 * x2 * x3^2 + 2 * B * h * m^2 * omega2 * omega3 * x1 * x2 * x3 -
                B * m^2 * omega1 * omega2 * x1^2 * x2 * x3 - B * m^2 * omega1 * omega2 * x2^3 * x3 -
                B * m^2 * omega1 * omega2 * x2 * x3^3 + B * m^2 * omega2 * omega3 * x1^3 * x2 +
                B * m^2 * omega2 * omega3 * x1 * x2^3 + B * m^2 * omega2 * omega3 * x1 * x2 * x3^2 -
                C^2 * F * omega2 * omega3 + C^2 * h^2 * m * omega1 * omega3 + 2 * C^2 * h * m * omega1 * omega3 * x3 +
                C^2 * m * omega1 * omega3 * x2^2 + C^2 * m * omega1 * omega3 * x3^2 -
                C^2 * m * omega2 * omega3 * x1 * x2 - C * F^2 * omega1 * omega3 - C * F * c_1 * omega1 -
                C * F * g * h * m * sin(alpha) + C * F * g * m * x2 * cos(alpha) * cos(beta) -
                C * F * g * m * x3 * sin(alpha) + C * F * h^2 * m * omega2 * omega3 - C * F * h * k2 * m +
                2 * C * F * h * m * omega2 * omega3 * x3 - C * F * k2 * m * x3 + C * F * k3 * m * x2 -
                C * F * m * omega1 * omega3 * x1 * x2 - C * F * m * omega2 * omega3 * x1^2 +
                C * F * m * omega2 * omega3 * x3^2 - C * c_1 * m * omega1 * x1 * x2 - C * c_2 * h^2 * m * omega2 -
                2 * C * c_2 * h * m * omega2 * x3 - C * c_2 * m * omega2 * x2^2 - C * c_2 * m * omega2 * x3^2 -
                C * g * h^3 * m^2 * sin(beta) * cos(alpha) - C * g * h^2 * m^2 * x1 * cos(alpha) * cos(beta) -
                3 * C * g * h^2 * m^2 * x3 * sin(beta) * cos(alpha) - C * g * h * m^2 * x1 * x2 * sin(alpha) -
                2 * C * g * h * m^2 * x1 * x3 * cos(alpha) * cos(beta) - C * g * h * m^2 * x2^2 * sin(beta) * cos(alpha) -
                3 * C * g * h * m^2 * x3^2 * sin(beta) * cos(alpha) - C * g * m^2 * x1 * x2 * x3 * sin(alpha) -
                C * g * m^2 * x1 * x3^2 * cos(alpha) * cos(beta) - C * g * m^2 * x2^2 * x3 * sin(beta) * cos(alpha) -
                C * g * m^2 * x3^3 * sin(beta) * cos(alpha) + C * h^3 * k1 * m^2 + 3 * C * h^2 * k1 * m^2 * x3 -
                C * h^2 * k3 * m^2 * x1 + C * h^2 * m^2 * omega1 * omega3 * x2^2 -
                C * h^2 * m^2 * omega2 * omega3 * x1 * x2 + C * h * k1 * m^2 * x2^2 + 3 * C * h * k1 * m^2 * x3^2 -
                C * h * k2 * m^2 * x1 * x2 - 2 * C * h * k3 * m^2 * x1 * x3 + 2 * C * h * m^2 * omega1 * omega3 * x2^2 * x3 -
                2 * C * h * m^2 * omega2 * omega3 * x1 * x2 * x3 + C * k1 * m^2 * x2^2 * x3 + C * k1 * m^2 * x3^3 -
                C * k2 * m^2 * x1 * x2 * x3 - C * k3 * m^2 * x1 * x3^2 + C * m^2 * omega1 * omega3 * x1^2 * x2^2 +
                C * m^2 * omega1 * omega3 * x2^4 + C * m^2 * omega1 * omega3 * x2^2 * x3^2 -
                C * m^2 * omega2 * omega3 * x1^3 * x2 - C * m^2 * omega2 * omega3 * x1 * x2^3 -
                C * m^2 * omega2 * omega3 * x1 * x2 * x3^2 + F^2 * h * m * omega1^2 * x1 -
                F^2 * h * m * omega2^2 * x1 + F^2 * m * omega1^2 * x1 * x3 - F^2 * m * omega1 * omega3 * x1^2 -
                F^2 * m * omega1 * omega3 * x2^2 - F^2 * m * omega2^2 * x1 * x3 - F * c_1 * m * omega1 * x1^2 -
                F * c_1 * m * omega1 * x2^2 - F * c_3 * h * m * omega3 * x1 - F * c_3 * m * omega3 * x1 * x3 +
                F * g * h * m^2 * x1 * x2 * sin(beta) * cos(alpha) - F * g * h * m^2 * x2^2 * sin(alpha) +
                F * g * m^2 * x1^2 * x2 * cos(alpha) * cos(beta) + F * g * m^2 * x1 * x2 * x3 * sin(beta) * cos(alpha) +
                F * g * m^2 * x2^3 * cos(alpha) * cos(beta) - F * g * m^2 * x2^2 * x3 * sin(alpha) +
                F * h^3 * m^2 * omega1^2 * x2 - F * h^3 * m^2 * omega2^2 * x2 + 3 * F * h^2 * m^2 * omega1^2 * x2 * x3 -
                F * h^2 * m^2 * omega1 * omega3 * x1 * x2 - 3 * F * h^2 * m^2 * omega2^2 * x2 * x3 +
                F * h^2 * m^2 * omega2 * omega3 * x2^2 - F * h * k1 * m^2 * x1 * x2 - F * h * k2 * m^2 * x2^2 +
                F * h * m^2 * omega1^2 * x1^2 * x2 + F * h * m^2 * omega1^2 * x2^3 + 3 * F * h * m^2 * omega1^2 * x2 * x3^2 -
                2 * F * h * m^2 * omega1 * omega3 * x1 * x2 * x3 - F * h * m^2 * omega2^2 * x1^2 * x2 -
                F * h * m^2 * omega2^2 * x2^3 - 3 * F * h * m^2 * omega2^2 * x2 * x3^2 + 2 * F * h * m^2 * omega2 * omega3 * x2^2 * x3 -
                F * k1 * m^2 * x1 * x2 * x3 - F * k2 * m^2 * x2^2 * x3 + F * k3 * m^2 * x1^2 * x2 + F * k3 * m^2 * x2^3 +
                F * m^2 * omega1^2 * x1^2 * x2 * x3 + F * m^2 * omega1^2 * x2^3 * x3 + F * m^2 * omega1^2 * x2 * x3^3 -
                F * m^2 * omega1 * omega3 * x1^3 * x2 - F * m^2 * omega1 * omega3 * x1 * x2^3 -
                F * m^2 * omega1 * omega3 * x1 * x2 * x3^2 - F * m^2 * omega2^2 * x1^2 * x2 * x3 - F * m^2 * omega2^2 * x2^3 * x3 -
                F * m^2 * omega2^2 * x2 * x3^3 + F * m^2 * omega2 * omega3 * x1^2 * x2^2 + F * m^2 * omega2 * omega3 * x2^4 +
                F * m^2 * omega2 * omega3 * x2^2 * x3^2 - c_1 * h^2 * m^2 * omega1 * x1 * x2 -
                2 * c_1 * h * m^2 * omega1 * x1 * x2 * x3 - c_1 * m^2 * omega1 * x1^3 * x2 - c_1 * m^2 * omega1 * x1 * x2^3 -
                c_1 * m^2 * omega1 * x1 * x2 * x3^2 - c_2 * h^2 * m^2 * omega2 * x2^2 -
                2 * c_2 * h * m^2 * omega2 * x2^2 * x3 - c_2 * m^2 * omega2 * x1^2 * x2^2 -
                c_2 * m^2 * omega2 * x2^4 - c_2 * m^2 * omega2 * x2^2 * x3^2 - c_3 * h^3 * m^2 * omega3 * x2 -
                3 * c_3 * h^2 * m^2 * omega3 * x2 * x3 - c_3 * h * m^2 * omega3 * x1^2 * x2 -
                c_3 * h * m^2 * omega3 * x2^3 - 3 * c_3 * h * m^2 * omega3 * x2 * x3^2 -
                c_3 * m^2 * omega3 * x1^2 * x2 * x3 - c_3 * m^2 * omega3 * x2^3 * x3 - c_3 * m^2 * omega3 * x2 * x3^3) /
            (A * B * C + A * B * m * x1^2 + A * B * m * x2^2 + A * C * h^2 * m + 2 * A * C * h * m * x3 +
                A * C * m * x1^2 + A * C * m * x3^2 + A * h^2 * m^2 * x1^2 + 2 * A * h * m^2 * x1^2 * x3 +
                A * m^2 * x1^4 + A * m^2 * x1^2 * x2^2 + A * m^2 * x1^2 * x3^2 + B * C * h^2 * m +
                2 * B * C * h * m * x3 + B * C * m * x2^2 + B * C * m * x3^2 + B * h^2 * m^2 * x2^2 +
                2 * B * h * m^2 * x2^2 * x3 + B * m^2 * x1^2 * x2^2 + B * m^2 * x2^4 + B * m^2 * x2^2 * x3^2 -
                C * F^2 - 2 * C * F * m * x1 * x2 + C * h^4 * m^2 + 4 * C * h^3 * m^2 * x3 + C * h^2 * m^2 * x1^2 +
                C * h^2 * m^2 * x2^2 + 6 * C * h^2 * m^2 * x3^2 + 2 * C * h * m^2 * x1^2 * x3 + 2 * C * h * m^2 * x2^2 * x3 +
                4 * C * h * m^2 * x3^3 + C * m^2 * x1^2 * x3^2 + C * m^2 * x2^2 * x3^2 + C * m^2 * x3^4 - F^2 * m * x1^2 -
                F^2 * m * x2^2 - 2 * F * h^2 * m^2 * x1 * x2 - 4 * F * h * m^2 * x1 * x2 * x3 -
                2 * F * m^2 * x1^3 * x2 - 2 * F * m^2 * x1 * x2^3 - 2 * F * m^2 * x1 * x2 * x3^2)
    end

# this function is shorter, but for (almost all?) inital conditions it didnt seems stable.
# this was not further investigated
#=
function w2_dot_again(alpha, beta, omega1, omega2, omega3, x1, x2, x3,
    k1, k2, k3, g, m, h, A, B, C, F, c1, c2, c3)
    return (-((A*C*h*k1*m - C*F*h*k2*m + C*h^3*k1*m^2 - C*c1*F*omega1 - A*C*c2*omega2 - C*c2*h^2*m*omega2 - 
    A^2*C*omega1*omega3 + A*C^2*omega1*omega3 - C*F^2*omega1*omega3 - A*C*h^2*m*omega1*omega3 + C^2*h^2*m*omega1*omega3 + 
    A*C*F*omega2*omega3 + B*C*F*omega2*omega3 - C^2*F*omega2*omega3 + C*F*h^2*m*omega2*omega3 - A*C*k3*m*x1 - C*h^2*k3*m^2*x1 + 
    F^2*h*m*omega1^2*x1 + A*F*h*m*omega1*omega2*x1 - B*F*h*m*omega1*omega2*x1 - F^2*h*m*omega2^2*x1 - c3*F*h*m*omega3*x1 + 
    A*h*k1*m^2*x1^2 - c1*F*m*omega1*x1^2 - A*c2*m*omega2*x1^2 - A^2*m*omega1*omega3*x1^2 + A*C*m*omega1*omega3*x1^2 - 
    F^2*m*omega1*omega3*x1^2 + A*F*m*omega2*omega3*x1^2 + B*F*m*omega2*omega3*x1^2 - C*F*m*omega2*omega3*x1^2 - A*k3*m^2*x1^3 + 
    C*F*k3*m*x2 + A*F*h*m*omega1^2*x2 + F*h^3*m^2*omega1^2*x2 + A^2*h*m*omega1*omega2*x2 - A*B*h*m*omega1*omega2*x2 + 
    A*h^3*m^2*omega1*omega2*x2 - B*h^3*m^2*omega1*omega2*x2 - A*F*h*m*omega2^2*x2 - F*h^3*m^2*omega2^2*x2 - 
    A*c3*h*m*omega3*x2 - c3*h^3*m^2*omega3*x2 - F*h*k1*m^2*x1*x2 + A*h*k2*m^2*x1*x2 - C*h*k2*m^2*x1*x2 - C*c1*m*omega1*x1*x2 - 
    c1*h^2*m^2*omega1*x1*x2 - C*F*m*omega1*omega3*x1*x2 - F*h^2*m^2*omega1*omega3*x1*x2 + B*C*m*omega2*omega3*x1*x2 - 
    C^2*m*omega2*omega3*x1*x2 + B*h^2*m^2*omega2*omega3*x1*x2 - C*h^2*m^2*omega2*omega3*x1*x2 + F*k3*m^2*x1^2*x2 + 
    F*h*m^2*omega1^2*x1^2*x2 + A*h*m^2*omega1*omega2*x1^2*x2 - B*h*m^2*omega1*omega2*x1^2*x2 - F*h*m^2*omega2^2*x1^2*x2 - 
    c3*h*m^2*omega3*x1^2*x2 - c1*m^2*omega1*x1^3*x2 - F*m^2*omega1*omega3*x1^3*x2 + B*m^2*omega2*omega3*x1^3*x2 - 
    C*m^2*omega2*omega3*x1^3*x2 + C*h*k1*m^2*x2^2 - F*h*k2*m^2*x2^2 - c1*F*m*omega1*x2^2 - A*c2*m*omega2*x2^2 - 
    C*c2*m*omega2*x2^2 - c2*h^2*m^2*omega2*x2^2 - A^2*m*omega1*omega3*x2^2 + C^2*m*omega1*omega3*x2^2 - 
    F^2*m*omega1*omega3*x2^2 - A*h^2*m^2*omega1*omega3*x2^2 + C*h^2*m^2*omega1*omega3*x2^2 + A*F*m*omega2*omega3*x2^2 + 
    B*F*m*omega2*omega3*x2^2 + F*h^2*m^2*omega2*omega3*x2^2 - A*k3*m^2*x1*x2^2 - c2*m^2*omega2*x1^2*x2^2 - 
    A*m^2*omega1*omega3*x1^2*x2^2 + C*m^2*omega1*omega3*x1^2*x2^2 + F*m^2*omega2*omega3*x1^2*x2^2 + F*k3*m^2*x2^3 + 
    F*h*m^2*omega1^2*x2^3 + A*h*m^2*omega1*omega2*x2^3 - B*h*m^2*omega1*omega2*x2^3 - F*h*m^2*omega2^2*x2^3 - 
    c3*h*m^2*omega3*x2^3 - c1*m^2*omega1*x1*x2^3 - F*m^2*omega1*omega3*x1*x2^3 + B*m^2*omega2*omega3*x1*x2^3 - 
    C*m^2*omega2*omega3*x1*x2^3 - c2*m^2*omega2*x2^4 - A*m^2*omega1*omega3*x2^4 + C*m^2*omega1*omega3*x2^4 + 
    F*m^2*omega2*omega3*x2^4 + A*C*k1*m*x3 - C*F*k2*m*x3 + 3*C*h^2*k1*m^2*x3 - 2*C*c2*h*m*omega2*x3 - 2*A*C*h*m*omega1*omega3*x3 + 
    2*C^2*h*m*omega1*omega3*x3 + 2*C*F*h*m*omega2*omega3*x3 - 2*C*h*k3*m^2*x1*x3 + F^2*m*omega1^2*x1*x3 + 
    A*F*m*omega1*omega2*x1*x3 - B*F*m*omega1*omega2*x1*x3 - F^2*m*omega2^2*x1*x3 - c3*F*m*omega3*x1*x3 + A*k1*m^2*x1^2*x3 + 
    A*F*m*omega1^2*x2*x3 + 3*F*h^2*m^2*omega1^2*x2*x3 + A^2*m*omega1*omega2*x2*x3 - A*B*m*omega1*omega2*x2*x3 + 
    3*A*h^2*m^2*omega1*omega2*x2*x3 - 3*B*h^2*m^2*omega1*omega2*x2*x3 - A*F*m*omega2^2*x2*x3 - 3*F*h^2*m^2*omega2^2*x2*x3 - 
    A*c3*m*omega3*x2*x3 - 3*c3*h^2*m^2*omega3*x2*x3 - F*k1*m^2*x1*x2*x3 + A*k2*m^2*x1*x2*x3 - C*k2*m^2*x1*x2*x3 - 
    2*c1*h*m^2*omega1*x1*x2*x3 - 2*F*h*m^2*omega1*omega3*x1*x2*x3 + 2*B*h*m^2*omega2*omega3*x1*x2*x3 - 
    2*C*h*m^2*omega2*omega3*x1*x2*x3 + F*m^2*omega1^2*x1^2*x2*x3 + A*m^2*omega1*omega2*x1^2*x2*x3 - 
    B*m^2*omega1*omega2*x1^2*x2*x3 - F*m^2*omega2^2*x1^2*x2*x3 - c3*m^2*omega3*x1^2*x2*x3 + C*k1*m^2*x2^2*x3 - 
    F*k2*m^2*x2^2*x3 - 2*c2*h*m^2*omega2*x2^2*x3 - 2*A*h*m^2*omega1*omega3*x2^2*x3 + 2*C*h*m^2*omega1*omega3*x2^2*x3 + 
    2*F*h*m^2*omega2*omega3*x2^2*x3 + F*m^2*omega1^2*x2^3*x3 + A*m^2*omega1*omega2*x2^3*x3 - B*m^2*omega1*omega2*x2^3*x3 - 
    F*m^2*omega2^2*x2^3*x3 - c3*m^2*omega3*x2^3*x3 + 3*C*h*k1*m^2*x3^2 - C*c2*m*omega2*x3^2 - A*C*m*omega1*omega3*x3^2 + 
    C^2*m*omega1*omega3*x3^2 + C*F*m*omega2*omega3*x3^2 - C*k3*m^2*x1*x3^2 + 3*F*h*m^2*omega1^2*x2*x3^2 + 
    3*A*h*m^2*omega1*omega2*x2*x3^2 - 3*B*h*m^2*omega1*omega2*x2*x3^2 - 3*F*h*m^2*omega2^2*x2*x3^2 - 
    3*c3*h*m^2*omega3*x2*x3^2 - c1*m^2*omega1*x1*x2*x3^2 - F*m^2*omega1*omega3*x1*x2*x3^2 + B*m^2*omega2*omega3*x1*x2*x3^2 - 
    C*m^2*omega2*omega3*x1*x2*x3^2 - c2*m^2*omega2*x2^2*x3^2 - A*m^2*omega1*omega3*x2^2*x3^2 + 
    C*m^2*omega1*omega3*x2^2*x3^2 + F*m^2*omega2*omega3*x2^2*x3^2 + C*k1*m^2*x3^3 + F*m^2*omega1^2*x2*x3^3 + 
    A*m^2*omega1*omega2*x2*x3^3 - B*m^2*omega1*omega2*x2*x3^3 - F*m^2*omega2^2*x2*x3^3 - c3*m^2*omega3*x2*x3^3 - 
    A*C*g*m*x1*cos(alpha)*cos(beta) - C*g*h^2*m^2*x1*cos(alpha)*cos(beta) - A*g*m^2*x1^3*cos(alpha)*cos(beta) + 
    C*F*g*m*x2*cos(alpha)*cos(beta) + F*g*m^2*x1^2*x2*cos(alpha)*cos(beta) - A*g*m^2*x1*x2^2*cos(alpha)*cos(beta) + 
    F*g*m^2*x2^3*cos(alpha)*cos(beta) - 2*C*g*h*m^2*x1*x3*cos(alpha)*cos(beta) - C*g*m^2*x1*x3^2*cos(alpha)*cos(beta) - 
    C*F*g*h*m*sin(alpha) + A*g*h*m^2*x1*x2*sin(alpha) - C*g*h*m^2*x1*x2*sin(alpha) - F*g*h*m^2*x2^2*sin(alpha) - 
    C*F*g*m*x3*sin(alpha) + A*g*m^2*x1*x2*x3*sin(alpha) - C*g*m^2*x1*x2*x3*sin(alpha) - F*g*m^2*x2^2*x3*sin(alpha) - 
    A*C*g*h*m*cos(alpha)*sin(beta) - C*g*h^3*m^2*cos(alpha)*sin(beta) - A*g*h*m^2*x1^2*cos(alpha)*sin(beta) + 
    F*g*h*m^2*x1*x2*cos(alpha)*sin(beta) - C*g*h*m^2*x2^2*cos(alpha)*sin(beta) - A*C*g*m*x3*cos(alpha)*sin(beta) - 
    3*C*g*h^2*m^2*x3*cos(alpha)*sin(beta) - A*g*m^2*x1^2*x3*cos(alpha)*sin(beta) + F*g*m^2*x1*x2*x3*cos(alpha)*sin(beta) - 
    C*g*m^2*x2^2*x3*cos(alpha)*sin(beta) - 3*C*g*h*m^2*x3^2*cos(alpha)*sin(beta) - C*g*m^2*x3^3*cos(alpha)*sin(beta))/
(-(A*B*C) + C*F^2 - A*C*h^2*m - B*C*h^2*m - C*h^4*m^2 - A*B*m*x1^2 - A*C*m*x1^2 + F^2*m*x1^2 - A*h^2*m^2*x1^2 - 
    C*h^2*m^2*x1^2 - A*m^2*x1^4 + 2*C*F*m*x1*x2 + 2*F*h^2*m^2*x1*x2 + 2*F*m^2*x1^3*x2 - A*B*m*x2^2 - B*C*m*x2^2 + 
    F^2*m*x2^2 - B*h^2*m^2*x2^2 - C*h^2*m^2*x2^2 - A*m^2*x1^2*x2^2 - B*m^2*x1^2*x2^2 + 2*F*m^2*x1*x2^3 - B*m^2*x2^4 - 
    2*A*C*h*m*x3 - 2*B*C*h*m*x3 - 4*C*h^3*m^2*x3 - 2*A*h*m^2*x1^2*x3 - 2*C*h*m^2*x1^2*x3 + 4*F*h*m^2*x1*x2*x3 - 
    2*B*h*m^2*x2^2*x3 - 2*C*h*m^2*x2^2*x3 - A*C*m*x3^2 - B*C*m*x3^2 - 6*C*h^2*m^2*x3^2 - A*m^2*x1^2*x3^2 - 
    C*m^2*x1^2*x3^2 + 2*F*m^2*x1*x2*x3^2 - B*m^2*x2^2*x3^2 - C*m^2*x2^2*x3^2 - 4*C*h*m^2*x3^3 - C*m^2*x3^4)))
end
=#

@inline function w3_dot(omega1_dot, omega2_dot, alpha, beta, omega1, omega2, omega3, x1, x2, x3,
        k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
        return (A * omega1 * omega2 - B * omega1 * omega2 + F * omega1^2 - F * omega2^2 -
                c_3 * omega3 + g * m * x1 * sin(alpha) + g * m * x2 * sin(beta) * cos(alpha) +
                h * m * omega1_dot * x1 + h * m * omega2_dot * x2 - k1 * m * x2 + k2 * m * x1 +
                m * omega1_dot * x1 * x3 + m * omega2_dot * x2 * x3) / (C + m * x1^2 + m * x2^2)
    end

@inline function w3_dot_mathematica(omega1_dot, omega2_dot, alpha, beta, omega1, omega2, omega3, x1, x2, x3,
    k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
    return (((A - B) * omega1 * omega2 + F * (omega1 - omega2) * (omega1 + omega2) - c_3 * omega3 +
            m * (k2 * x1 - k1 * x2 + (omega1_dot * x1 + omega2_dot * x2) * (h + x3)) + g * m * x1 * sin(alpha) + g * m * x2 * cos(alpha) * sin(beta)) /
            (C + m * (x1^2 + x2^2)))
end

function dydt_torque(y, parameter, t)
    alpha, beta, gamma, omega1, omega2, omega3 = y
    (g, m, h, A, B, C, F, a, b, c, c_1, c_2, c_3) = parameter

    # calculate the derivatives of the angles alpha, beta, gamma
    alpha_dot = omega1 * cos(beta) + omega3 * sin(beta)
    beta_dot = (omega1 * sin(beta) - omega3 * cos(beta)) * tan(alpha) + omega2
    gamma_dot = (-omega1 * sin(beta) + omega3 * cos(beta)) / cos(alpha)

    # calculate k_i 
    mu1 = cos(alpha) * sin(beta)
    mu2 = -sin(alpha)
    mu3 = -cos(alpha) * cos(beta)
    p̃ = sqrt((a * mu1)^2 + (b * mu2)^2 + (c * mu3)^2)

    x1 = a^2 * mu1 / p̃
    x2 = b^2 * mu2 / p̃
    x3 = c^2 * mu3 / p̃

    # formels from CAS/rattleback_equations_x_dot.ipynb
    x1_dot = a^2 * (alpha_dot * b^2 * mu2 * sin(beta) + b^2 * beta_dot * mu3 * cos(alpha)^2 + b^2 * beta_dot * cos(alpha) * cos(beta) + beta_dot * c^2 * cos(alpha)^3 * cos(beta)) / p̃^3
    x2_dot = b^2 * (4 * a^2 * alpha_dot * cos(2 * beta) - 4 * a^2 * alpha_dot + a^2 * beta_dot * cos(2 * alpha - 2 * beta) - a^2 * beta_dot * cos(2 * alpha + 2 * beta) - 4 * alpha_dot * c^2 * cos(2 * beta) - 4 * alpha_dot * c^2 - beta_dot * c^2 * cos(2 * alpha - 2 * beta) + beta_dot * c^2 * cos(2 * alpha + 2 * beta)) * cos(alpha) / (8 * p̃^3)
    x3_dot = c^2 * (a^2 * beta_dot * mu1 * cos(alpha)^2 + alpha_dot * b^2 * sin(alpha) * cos(beta) - b^2 * beta_dot * mu1 * cos(alpha)^2 + b^2 * beta_dot * mu1) / p̃^3

    v1 = omega3 * x2 - omega2 * (x3 + h)
    v2 = omega1 * (x3 + h) - omega3 * x1
    v3 = omega2 * x1 - omega1 * x2

    k1 = omega2 * (v3 - x3_dot) - omega3 * (v2 - x2_dot)
    k2 = omega3 * (v1 - x1_dot) - omega1 * (v3 - x3_dot)
    k3 = omega1 * (v2 - x2_dot) - omega2 * (v1 - x1_dot)



    omega1_dot = w1_dot(alpha, beta, omega1, omega2, omega3, x1, x2, x3,
        k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
    omega2_dot = w2_dot(alpha, beta, omega1, omega2, omega3, x1, x2, x3,
        k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
    omega3_dot = w3_dot_mathematica(omega1_dot, omega2_dot, alpha, beta, omega1, omega2, omega3,
        x1, x2, x3, k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
    return [alpha_dot, beta_dot, gamma_dot, omega1_dot, omega2_dot, omega3_dot]
end

# define the ODE for Kuypers
function dydt_torque!(dy, y, parameter, t)
    alpha, beta, gamma, omega1, omega2, omega3 = y
    (g, m, h, A, B, C, F, a, b, c, c_1, c_2, c_3) = parameter

    # calculate the derivatives of the angles alpha, beta, gamma
    alpha_dot = omega1 * cos(beta) + omega3 * sin(beta)
    beta_dot = (omega1 * sin(beta) - omega3 * cos(beta)) * tan(alpha) + omega2
    gamma_dot = (-omega1 * sin(beta) + omega3 * cos(beta)) / cos(alpha)

    # calculate k_i 
    mu1 = cos(alpha) * sin(beta)
    mu2 = -sin(alpha)
    mu3 = -cos(alpha) * cos(beta)
    p̃ = sqrt((a * mu1)^2 + (b * mu2)^2 + (c * mu3)^2)

    x1 = a^2 * mu1 / p̃
    x2 = b^2 * mu2 / p̃
    x3 = c^2 * mu3 / p̃

    # formels from rattleback/bewegungsgleichungen_rattleback_x_dot.ipynb
    x1_dot = a^2 * (alpha_dot * b^2 * mu2 * sin(beta) + b^2 * beta_dot * mu3 * cos(alpha)^2 + b^2 * beta_dot * cos(alpha) * cos(beta) + beta_dot * c^2 * cos(alpha)^3 * cos(beta)) / p̃^3
    x2_dot = b^2 * (4 * a^2 * alpha_dot * cos(2 * beta) - 4 * a^2 * alpha_dot + a^2 * beta_dot * cos(2 * alpha - 2 * beta) - a^2 * beta_dot * cos(2 * alpha + 2 * beta) - 4 * alpha_dot * c^2 * cos(2 * beta) - 4 * alpha_dot * c^2 - beta_dot * c^2 * cos(2 * alpha - 2 * beta) + beta_dot * c^2 * cos(2 * alpha + 2 * beta)) * cos(alpha) / (8 * p̃^3)
    x3_dot = c^2 * (a^2 * beta_dot * mu1 * cos(alpha)^2 + alpha_dot * b^2 * sin(alpha) * cos(beta) - b^2 * beta_dot * mu1 * cos(alpha)^2 + b^2 * beta_dot * mu1) / p̃^3

    v1 = omega3 * x2 - omega2 * (x3 + h)
    v2 = omega1 * (x3 + h) - omega3 * x1
    v3 = omega2 * x1 - omega1 * x2

    k1 = omega2 * (v3 - x3_dot) - omega3 * (v2 - x2_dot)
    k2 = omega3 * (v1 - x1_dot) - omega1 * (v3 - x3_dot)
    k3 = omega1 * (v2 - x2_dot) - omega2 * (v1 - x1_dot)



    omega1_dot = w1_dot(alpha, beta, omega1, omega2, omega3, x1, x2, x3,
        k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
    omega2_dot = w2_dot(alpha, beta, omega1, omega2, omega3, x1, x2, x3,
        k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)
    omega3_dot = w3_dot_mathematica(omega1_dot, omega2_dot, alpha, beta, omega1, omega2, omega3,
        x1, x2, x3, k1, k2, k3, g, m, h, A, B, C, F, c_1, c_2, c_3)


    dy[1] = alpha_dot
    dy[2] = beta_dot
    dy[3] = gamma_dot
    dy[4] = omega1_dot
    dy[5] = omega2_dot
    dy[6] = omega3_dot

    return nothing

end

# the calculate A, B, C, F from the principal moments of inertia and angle
function hauptachsentransformation(I, δ)
    I1, I2, I3 = I
    A = 0.5 * (I1 + I2) + 0.5 * (I1 - I2) * cos(2 * δ)
    B = 0.5 * (I1 + I2) - 0.5 * (I1 - I2) * cos(2 * δ)
    F = -0.5 * (I1 - I2) * sin(2 * δ)
    C = I3
    return [A, B, C, F]
end


# calculate the energy of the system for kuypers
function calculate_energy_kuypers(sol,parameter)
    g, m, h, A, B, C, F, a, b, c = parameter
    E_ges = zeros(length(sol))
    E_kin = zeros(length(sol))
    E_rot = zeros(length(sol))
    E_pot = zeros(length(sol))
    for i in 1:length(sol)
        alpha, beta, gamma, omega1, omega2, omega3 = sol[:, i]
        # calculate velocities
        mu1 = cos(alpha) * sin(beta)
        mu2 = -sin(alpha)
        mu3 = -cos(alpha) * cos(beta)
        p̃ = sqrt((a * mu1)^2 + (b * mu2)^2 + (c * mu3)^2)

        x1 = a^2 * mu1 / p̃
        x2 = b^2 * mu2 / p̃
        x3 = c^2 * mu3 / p̃


        v1 = omega3 * x2 - omega2 * (x3 + h)
        v2 = omega1 * (x3 + h) - omega3 * x1
        v3 = omega2 * x1 - omega1 * x2

        E_kin[i] = 0.5 * m * (v1^2 + v2^2 + v3^2)

        E_rot[i] = 0.5 * (A * omega1^2 + B * omega2^2 + C * omega3^2 - 2 * F * omega1 * omega2)

        # projektion des Vektors der von Schwerpunkt zu Auflagepunkt zeigt auf n_3^k 
        # E_pot = m * g * dot(x, n_3) 
        E_pot[i] = m * g * (x1 * cos(alpha) * sin(beta) - x2 * sin(alpha) - (x3 + h) * cos(alpha) * cos(beta))


        E_ges[i] = E_kin[i] + E_rot[i] + E_pot[i]
    end
    return E_ges, E_kin, E_rot, E_pot
end



# set up the original problem from Kuypers (used in chapter 7.1)
function set_up_problem_kuypers_wackelstein_1()
    g = 9.81
    m = 1.0
    h = 1e-3
    I = [2e-4, 1.6e-3, 1.7e-3]
    δ = 0.1
    A, B, C, F = hauptachsentransformation(I, δ)
    a, b, c = 0.2, 0.03, 0.02
    c_i = 1e-4
    c_1, c_2, c_3 = c_i, c_i, c_i
    parameter = (g=g, m=m, h=h, A=A, B=B, C=C, F=F, a=a, b=b, c=c, c_1=c_1, c_2=c_2, c_3=c_3)

    y0 = [-0.02, -0.02, 0.0,
        0.0, 0.0, -8.0]

    tspan = (0.0, 8.0)

    return y0, tspan, parameter
end

# make the plot in chapter 7.1
function plot_kuypers_wackelstein_1_Abb4(sol, parameter)
    Eges, _, _, _ = calculate_energy_kuypers(sol, parameter)

    p1 = plot(sol.t, sol[3, :], label="Gamma(t)", legend=:right)
    p2 = plot(sol.t, sol[1, :], label="Alpha(t)", legend=:topright)
    p3 = plot(sol.t, sol[2, :], label="Beta(t)", legend=:topright)
    p4 = plot(sol.t, Eges, label="Energie(t)", legend=:topright,xlabel="t",)

    fontsize = 14
    fontsize2 = 18

    p = plot(
        p1,
        p2,
        p3,
        p4,
        layout=(4, 1),
        lw = 2,
        xlims=(first(sol.t), last(sol.t)),
        grid = true,
        gridalpha = 1,
        gridstyle = :dot,
        #linecolor = :blue3,
        linecolor = RGB(18/255, 18/255, 251/255),
        #suptitle  = "Reproduzierung des Plots von Kuypers in Wackelstein_1.pdf (Abb4)",
        size = (1200, 800),
        dpi = 300,
        titlefont = font(fontsize2,"Computer Modern"),
        guidefont = font("Computer Modern"),
        legendfont = font(fontsize,"Computer Modern"),
        xtickfontsize=fontsize,
        ytickfontsize=fontsize,
        ylabelfontsize=fontsize2,
        xlabelfontsize=fontsize2,

    )
    
    return p
end



"""
    function used to calculate p_γ, which is a constant according to the Lagrange equations
    but is not constant in the simulation using the torque Kuypers approach

    the function was used during the work on the thesis when it was not clear if/what was wrong
    with the lagrange formalism
"""
function p_gamma(alpha, beta, omega1, omega2, omega3, m, h, A, B, C, F)
    return ((-(F * omega1) + B * omega2) * sin(alpha) +
            cos(alpha) * (C * omega3 * cos(beta) + (-(A * omega1) + F * omega2) * sin(beta)) +
            (m * ((a - b) * (a + b) * sin(2 * alpha) * sin(beta) *
                    (b^2 * omega1 * sin(alpha) + a^2 * omega2 * cos(alpha) * sin(beta)) -
                    cos(alpha) * sin(beta) * ((a - c) * (a + c) * cos(alpha) * cos(beta) +
                                            h * sqrt(b^2 * sin(alpha)^2 +
                                                        cos(alpha)^2 * (c^2 * cos(beta)^2 + a^2 * sin(beta)^2))) *
                    (-2 * cos(alpha) * (c^2 * omega1 * cos(beta) + a^2 * omega3 * sin(beta)) +
                    2 * h * omega1 * sqrt(b^2 * sin(alpha)^2 +
                                            cos(alpha)^2 * (c^2 * cos(beta)^2 + a^2 * sin(beta)^2))) +
                    2 * sin(alpha) * ((b - c) * (b + c) * cos(alpha) * cos(beta) +
                                    h * sqrt(b^2 * sin(alpha)^2 +
                                                cos(alpha)^2 * (c^2 * cos(beta)^2 + a^2 * sin(beta)^2))) *
                    (-(c^2 * omega2 * cos(alpha) * cos(beta)) + b^2 * omega3 * sin(alpha) +
                    h * omega2 * sqrt(b^2 * sin(alpha)^2 +
                                        cos(alpha)^2 * (c^2 * cos(beta)^2 + a^2 * sin(beta)^2))))) /
            (2 * (b^2 * sin(alpha)^2 +
                    cos(alpha)^2 * (c^2 * cos(beta)^2 + a^2 * sin(beta)^2))))
end


# set up the original problem from schoemer in the kuypers formalism
function set_up_original_problem_for_torque(ω = [0.01, -0.02, -2.0])
    # Geometric parameters of the rattleback
    a = 0.5 * (1 + √(5))
    b = 1
    c = 1


    g = 9.81
    # calculating the total mass and the height/ distance center of mass and 0
    mE = 0.1
    mp = 0.05
    m = mE + 2 * mp
    h = mE / m * 3 * c/ 8

    A, B, C, F =  0.05796875, 0.0903294297749979, 0.12236067977499791, 0.025

    c_i = 0.0
    c_1, c_2, c_3 = c_i, c_i, c_i
    parameter = (g=g, m=m, h=h, A=A, B=B, C=C, F=F, a=a, b=b, c=c, c_1=c_1, c_2=c_2, c_3=c_3)

    y0 = [0.0, 0.0, 0.0, ω...]


    return y0,  parameter
end



# plot the omega values of the kuypers solution
function plot_kuypers_omega(sol, pos = :topleft)
    lw = 2
    fontsize = 14
    fontsize2 = 18
    p = plot(sol.t, sol[4,:], label=L"\omega_1",
    xguide = L"t",
    yguide = L"\omega_i",
    title = L"\omega(t)" *" - Wackelstein Kuypers",
    legend = pos,
    size = (700, 500),
    grid = true,
    xlims=(first(sol.t), last(sol.t)),
    dpi = 300,
    lw = lw,
    titlefont = font(fontsize2,"Computer Modern"),
    guidefont = font("Computer Modern"),
    legendfont = font(fontsize,"Computer Modern"),
    xtickfontsize=fontsize,
    ytickfontsize=fontsize,
    ylabelfontsize=fontsize2,
    xlabelfontsize=fontsize2,
    )

    plot!(sol.t, sol[5,:], label=L"\omega_2", lw = lw)
    plot!(sol.t, sol[6,:], label=L"\omega_3", lw = lw)

    
    
    return p
end

# plot the angles of the kuypers solution
function plot_kuypers_angles(sol)

    p1 = plot(sol.t, sol[3,:], label="Gamma(t)")
    p2 = plot(sol.t, sol[1,:], label="Alpha(t)")
    p3 = plot(sol.t, sol[2,:], label="Beta(t)")
    p = plot(p1, p2, p3, layout=(3,1), legend = :topright, 
            xlims=(first(sol.t), last(sol.t)), 
            gridlinewidth=1,
            suptitle  = "Winkel(t)"
    )
    return p
end




# plot the difference of the omega values of the kuypers solution and the schömer red. solution
function plot_differance_omega(sol_K, sol_S)
    if length(sol_K) != length(sol_S)
        throw(ArgumentError("sol_K and sol_S must have the same length"))
    end


    omega1 = zeros(length(sol_S))
    omega2 = zeros(length(sol_S))
    omega3 = zeros(length(sol_S))

    for i in 1:length(sol_S)
        c, q, w = unpack_rattleback_reduced(sol_S[i])
        R = to_rotation_matrix(q)
        ω = R' * w

        omega1[i] = ω[1]
        omega2[i] = ω[2]
        omega3[i] = ω[3]

    end

    p = plot(sol_S.t, omega1 .- sol_K[4,:], label=L"\Delta \omega_1",
        xguide = L"t",
        yguide = L"\omega_i",
        title = L"\omega(t)" *" - Unterschied Schöemer / Kuypers",
        legend = :topleft,
        size = (800, 600),
        grid = true,
        xlims=(first(sol_S.t), last(sol_S.t)),
        )

        plot!(sol_S.t, omega2 .- sol_K[5,:], label=L"\Delta \omega_2")
        plot!(sol_S.t, omega3 .- sol_K[6,:], label=L"\Delta \omega_3")

    return p
end
