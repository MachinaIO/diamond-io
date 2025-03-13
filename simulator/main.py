#!/usr/bin/env sage -python
# import sys
# from sage.all import *
from estimator.estimator import *
from estimator.estimator.lwe_parameters import *
from estimator.estimator.nd import *
import math
import datetime
import os
from decimal import Decimal, getcontext

getcontext().prec = 100


def log_params_to_file(
    secpar: int,
    n: int,
    d: int,
    alpha: int,
    input_size: int,
    m_polys: list[list[int]],
    q: int,
    stddev_e: int,
    p: int,
    estimated_secpar: float,
):
    """
    Log parameters to params.log file
    """
    # Calculate log_q and log_p
    log_q = math.log2(q)
    log_p = math.log2(p)

    # Get current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Format m_polys as a string
    m_polys_str = str(m_polys).replace(" ", "")

    # Create log entry with key information
    log_entry = (
        f"{current_date}, "
        f"secpar={secpar}, "
        f"n={n}, "
        f"d={d}, "
        f"alpha={alpha}, "
        f"input_size={input_size}, "
        f"m_polys={m_polys_str}, "
        f"q={q}, "
        f"log_q={log_q}, "
        f"stddev_e={stddev_e}, "
        f"p={p}, "
        f"log_p={log_p}, "
        f"estimated_secpar={estimated_secpar}\n"
    )

    # Append to params.log file
    with open("params.log", "a") as f:
        f.write(log_entry)

    print(f"Parameters logged to params.log")


def find_params(
    target_secpar: int,
    n: int,
    d: int,
    input_size: int,
    m_polys: list[list[int]],
):

    min_alpha_k = -500
    max_alpha_k = -1
    min_q_k = target_secpar + 2
    max_q_k = 800
    middle_q_k = math.floor((min_q_k + max_q_k) // 2)

    found_alphas = []
    # Binary search for alpha
    while min_alpha_k + 1 < max_alpha_k:
        alpha_k = (min_alpha_k + max_alpha_k) / 2
        print(f"min_alpha_k: {min_alpha_k}")
        print(f"max_alpha_k: {max_alpha_k}")
        print(f"alpha_k: {alpha_k}")
        q = 2 ** (middle_q_k)
        stddev_e = 2 ** (middle_q_k + alpha_k)
        estimated_secpar = estimate_secpar(d * n, q, Binary, stddev_e)
        print("target_secpar:", target_secpar)
        print("estimated_secpar:", estimated_secpar)
        if target_secpar > estimated_secpar:
            print(
                f"target_secpar {target_secpar} > estimated_secpar {estimated_secpar}"
            )
            min_alpha_k = alpha_k
        else:
            found_alphas.append(alpha_k)
            print(f"found alpha_k: {alpha_k}")
            max_alpha_k = alpha_k
    if found_alphas == []:
        raise ValueError("alpha is not found after binary search")
    alpha = 2 ** min(found_alphas)
    print(f"found alpha: {alpha}")

    found_params = []
    iters = 0
    while min_q_k + 1 < max_q_k and iters < 100:
        iters += 1
        q_k = math.floor((min_q_k + max_q_k) // 2)
        print(f"min_q_k: {min_q_k}")
        print(f"max_q_k: {max_q_k}")
        print(f"q_k: {q_k}")
        q = 2**q_k
        print(f"q: {q}")
        stddev_e = alpha * q
        estimated_secpar = estimate_secpar(d * n, q, Binary, stddev_e)
        print("target_secpar:", target_secpar)
        print("estimated_secpar:", estimated_secpar)
        if target_secpar > estimated_secpar:
            print(
                f"target_secpar {target_secpar} > estimated_secpar {estimated_secpar}"
            )
            max_q_k = q_k
            continue
        try:
            p = find_p(target_secpar, n, q, d, stddev_e, input_size, m_polys)
            print(f"found p: {p}")
            max_q_k = q_k
            found_params.append((q, stddev_e, p, estimated_secpar))
        except ValueError as e:
            print(f"ValueError: {e}")
            min_q_k = q_k
    if found_params == []:
        raise ValueError("p is not found after binary search")
    # minimum q in found_params
    return min(found_params, key=lambda x: x[1]) + (alpha,)


def find_p(
    secpar: int,
    n: int,
    q: int,
    d: int,
    stddev_e: int,
    input_size: int,
    m_polys: list[list[int]],
):
    log_q = math.ceil(math.log2(q))
    if log_q - 2 < secpar:
        raise ValueError(
            f"log2(q) - 2 should be larger than secpar (given log2(q): {log_q}, secpar: {secpar})"
        )
    stddev_b = compute_stddev_b(n, log_q, d)

    final_err = bound_final_error(
        secpar, n, log_q, d, stddev_e, stddev_b, input_size, m_polys
    )

    # Convert final_err to Decimal for high precision
    final_err_decimal = Decimal(str(final_err))

    # Calculate log2 using Decimal
    # Handle infinity or very large numbers
    if math.isinf(final_err):
        # If final_err is infinity, set log_final_err to a very large value
        log_final_err = 10000  # A very large value that will likely fail the next check
        print(
            f"Warning: final_err is infinity, setting log_final_err to {log_final_err}"
        )
    elif final_err_decimal > 0:
        try:
            # Try to calculate log2 using Decimal
            log_final_err = math.ceil(
                float(final_err_decimal.ln() / Decimal("0.693147180559945"))
            )  # ln(2)
        except (OverflowError, ValueError):
            # If there's an overflow, set log_final_err to a very large value
            log_final_err = (
                10000  # A very large value that will likely fail the next check
            )
            print(
                f"Warning: final_err is too large for ln(), setting log_final_err to {log_final_err}"
            )
    else:
        raise ValueError(f"Cannot calculate log2 of non-positive value: {final_err}")

    if log_q - 2 <= log_final_err:
        raise ValueError(
            f"log2(final error) should be less than log2(q) - 2 (given log2(final error): {log_final_err}, log2(q): {log_q})"
        )

    # Use Decimal for high precision arithmetic
    # Convert to Decimal for high precision calculations
    q_decimal = Decimal(q)
    p_decimal = (q_decimal / Decimal("4")).to_integral_exact(
        rounding="ROUND_CEILING"
    ) - final_err_decimal
    p = int(p_decimal)

    # Calculate log2(p) using Decimal
    if p_decimal > 0:
        log_p = math.ceil(float(p_decimal.ln() / Decimal("0.693147180559945")))  # ln(2)
    else:
        raise ValueError(f"Cannot calculate log2 of non-positive value: {p}")

    if log_p - log_final_err < secpar:
        raise ValueError(
            f"p - error should be larger than 2^secpar (given p: {p}, secpar: {secpar})"
        )

    return p


def compute_stddev_b(
    n: int,
    log_q: int,
    d: int,
):
    c_0 = 1.3
    c_1 = 4.7
    sigma = 4.578
    return (
        c_0
        * sigma
        * (3 * sigma)
        * (sqrt_ceil((d + 1) * n * log_q) + sqrt_ceil(2 * n) + c_1)
    )


def estimate_secpar(n: int, q: int, s_dist: NoiseDistribution, stddev: int):
    params = LWEParameters(n, q, s_dist, DiscreteGaussian(stddev))
    estim = LWE.estimate.rough(params)
    # print(estim)
    min_secpar = math.ceil(math.log2(min(val["rop"] for val in estim.values())))
    # print(min_secpar)
    return min_secpar


def bound_final_error(
    secpar: int,
    n: int,
    log_q: int,
    d: int,
    stddev_e: int,
    stddev_b: int,
    input_size: int,
    m_polys: list[list[int]],
):
    # Convert all inputs to Decimal for high precision
    # secpar_d = Decimal(secpar)
    n_d = Decimal(n)
    log_q_d = Decimal(log_q)
    d_d = Decimal(d)
    stddev_e_d = Decimal(stddev_e)
    stddev_b_d = Decimal(stddev_b)

    # Calculate intermediate values with Decimal
    m_d = (d_d + Decimal(1)) * log_q_d
    m_b_d = (d_d + Decimal(1)) * (log_q_d + Decimal(2))
    sqrt_secpar_d = Decimal(sqrt_ceil(secpar))

    # Use Decimal for all calculations to maintain precision
    scale_coeff_d = (n_d * m_b_d * stddev_b_d) ** Decimal(2) * sqrt_secpar_d
    bound_p_d = stddev_e_d * sqrt_secpar_d
    bound_c_d = stddev_e_d * sqrt_secpar_d

    for _ in range(input_size):
        bound_v_d = (scale_coeff_d ** Decimal(2)) * bound_p_d
        bound_c_d = n_d * m_d * bound_c_d + bound_v_d
        bound_p_d = bound_v_d

    # Evaluate each polynomial in m_polys at the value of m using Decimal
    evaluated_polys_d = []
    for poly in m_polys:
        # Evaluate polynomial: sum(coeff * m^i for i, coeff in enumerate(poly))
        result_d = Decimal(0)
        for i, coeff in enumerate(poly):
            result_d += Decimal(coeff) * (m_d ** Decimal(i))
        evaluated_polys_d.append(result_d)

    # Find max value using Decimal
    if evaluated_polys_d:
        max_evaluated_poly_d = max(evaluated_polys_d)
    else:
        max_evaluated_poly_d = Decimal(1)  # Default if no polynomials

    bound_c_final_d = bound_c_d * max_evaluated_poly_d + stddev_b_d * sqrt_secpar_d
    bound_v_final_d = bound_p_d * scale_coeff_d

    # Return the final result as a Decimal
    return bound_c_final_d + bound_v_final_d


def sqrt_ceil(x):
    return math.ceil(math.sqrt(x))


# def bound_from_stddev(stddev: int, secpar: int):
#     return math.ceil(stddev * math.ceil(math.sqrt(secpar)))


# def derive_auto_params(
#     n: int, n_t: int, q: int, sigma_e: int, t: int = 2, target_secpar: int = 80
# ):
#     log_q = math.ceil(math.log2(q))
#     print("log_q:", log_q)
#     m = 2 * log_q
#     m_t = (n_t + 1) * log_q
#     m_b = 2 + math.ceil(log_q / math.log2(t))
#     # log_p = math.ceil(math.log2(p))
#     # m_b = n_b * log_p
#     # sigma_b = math.ceil(2 * math.sqrt(n * log_p))
#     secpar_s = math.ceil(output_secpar(n, q, Binary, sigma_e))
#     print("secpar_n:", secpar_s)
#     secpar_t = math.ceil(output_secpar(n_t, q, Binary, sigma_e))
#     print("secpar_t", secpar_t)
#     secpar = min(target_secpar, secpar_s, secpar_t)
#     print("secpar:", secpar)
#     sigma_b = math.ceil(math.sqrt(math.log(2 * n * 2 ** (80)) / 3.14))
#     print("sigma_b:", sigma_b)
#     print("sigma_b**2", sigma_b**2)
#     print("n * math.ceil(log_q / math.log2(t)))", n * math.ceil(log_q / math.log2(t)))
#     print("math.sqrt(2 * n)", math.sqrt(2 * n))
#     sigma_b = math.ceil(
#         1.3
#         * (t + 1)
#         * sigma_b**2
#         * (math.sqrt(n * math.ceil(log_q / math.log2(t))) + math.sqrt(2 * n) + 4.7)
#     )
#     print("sigma_b:", sigma_b)
#     # secpar_t = math.ceil(output_secpar(n_t, q, Binary, sigma_e))
#     # print("secpar_t", secpar_t)

#     # secpar_n_s = math.ceil(output_secpar(n_s, m_s, q, sigma_e))
#     # print("secpar_n_s:", secpar_n_s)
#     # secpar_b = math.ceil(output_secpar(n_b, q, UniformMod(q), sigma_b))
#     # print("secpar_b:", secpar_b)
#     bound_e = bound_from_stddev(sigma_e, secpar)
#     print("bound_e:", bound_e)
#     if bound_e <= 1:
#         raise ValueError("bound_e should be larger than 1")
#     bound_b = bound_from_stddev(sigma_b, secpar)
#     print("bound_b:", bound_b)
#     if bound_b <= 1:
#         raise ValueError("bound_b should be larger than 1")
#     b_c_frac = math.ceil(
#         (
#             Decimal(sigma_e)
#             * math.ceil((math.sqrt(n * m_b) * sigma_b) ** 2)
#             * math.ceil(secpar**2)
#             * math.ceil(math.sqrt(secpar))
#         )
#         / Decimal(math.ceil((math.sqrt(n * m_b) * sigma_b) ** 2) * secpar - n * m)
#     )
#     # b_c_frac = sigma_e * math.ceil(math.sqrt(secpar))
#     return {
#         "q": q,
#         "log_q": log_q,
#         "n": n,
#         "m": m,
#         "n_t": n_t,
#         "m_t": m_t,
#         "m_b": m_b,
#         "secpar": secpar,
#         "sigma_e": sigma_e,
#         "sigma_b": sigma_b,
#         "bound_e": bound_e,
#         "bound_b": bound_b,
#         "b_c_frac": b_c_frac,
#         # "p": p,
#     }


# def estimate_noise_norm(params, input: int, output: int, depth: int, scale: int = 0):
#     q = params["q"]
#     log_q = params["log_q"]
#     n = params["n"]
#     m = params["m"]
#     n_t = params["n_t"]
#     m_t = params["m_t"]
#     m_b = params["m_b"]
#     secpar = params["secpar"]
#     sigma_e = params["sigma_e"]
#     sigma_b = params["sigma_b"]
#     bound_e = params["bound_e"]
#     bound_b = params["bound_b"]
#     b_c_frac = params["b_c_frac"]
#     print("log b_c_frac", math.log2(b_c_frac))
#     print("log m", math.log2(m))
#     print(
#         "log (sigma_b * sigma_b * secpar)",
#         math.log2((n * m_b * sigma_b * sigma_b * secpar)),
#     )
#     print(
#         "log (sigma_b * sigma_b * secpar) ** input",
#         math.ceil(math.log2((n * m_b * sigma_b * sigma_b * secpar) ** input)),
#     )
#     b_c = (bound_e - b_c_frac) * (m**input) + b_c_frac * math.ceil(
#         (n * m_b * sigma_b * sigma_b * secpar) ** input
#     )
#     # print("b_c", b_c)
#     print("log b_c", math.log2(b_c) if b_c > 0 else 0)
#     input_ext = 1 + input + 256
#     # 1 + input + m * (256 + 2 * secpar) * (n + 1) * log_q
#     # b_f_exp1 = math.ceil(
#     #     depth * math.ceil(math.log2(m_s)) * math.ceil(math.log2(log_q))
#     #     + (math.ceil(math.log2(log_q)) ** 2)
#     #     + 2
#     # )
#     b_f_exp1 = depth
#     print("b_f_exp1", b_f_exp1)
#     # print("log (input_ext + n + 2)", math.log2((input_ext + n + 2)))
#     # print("log ((n + 1) * log_q)", math.log2((n + 1) * log_q))
#     # print(
#     #     "log math.ceil((m_s + 2) ** b_f_exp1)",
#     #     math.log2(math.ceil((m_s + 2) ** b_f_exp1)),
#     # )
#     b_f_term1 = (
#         b_c
#         * math.sqrt(input_ext + 2)
#         * ((n_t + 1) * log_q)
#         * math.ceil((n * m + 2) ** b_f_exp1)
#     )
#     print("log b_f_term1", math.log2(b_f_term1))
#     b_f_term2 = bound_e * log_q * ((m_t + 2) ** (depth + 1))
#     print("log b_f_term2", math.log2(b_f_term2))
#     b_f = b_f_term1 + b_f_term2
#     print("log b_f", math.log2(b_f))
#     print(
#         "log (n * m_b * sigma_b) ** (2 * input + 1)",
#         math.log2((n * m_b * sigma_b) ** (2 * input + 1)),
#     )
#     print("log secpar ** (input + 1)", math.log2(secpar ** (input + 1)))
#     b_z = b_f + sigma_e * ((math.sqrt(n * m_b) * sigma_b) ** (2 * input + 1)) * (
#         secpar ** (input + 1)
#     )
#     print("log b_z", math.log2(b_z))
#     # print("b_z", b_z)
#     # print("log2_b_z", math.log2(b_z))
#     # print("b_z + 2^**(2*secpar)", b_z * (2 ** (2 * secpar)))
#     return b_z
#     # b_f =
#     # b_c = (params.bound_b - )


# def estimate_obf_size(params, input: int, output: int, depth: int):
#     q = params["q"]
#     log_q = params["log_q"]
#     n = params["n"]
#     m = params["m"]
#     n_t = params["n_t"]
#     m_t = params["m_t"]
#     m_b = params["m_b"]
#     secpar = params["secpar"]
#     sigma_e = params["sigma_e"]
#     sigma_b = params["sigma_b"]
#     bound_e = params["bound_e"]
#     bound_b = params["bound_b"]

#     # bits
#     size = 0
#     # h (R and A matrixes are generated by a random oracle)
#     size += 256
#     # FHE encryption
#     # size += m * (256 + 2 * secpar) * (n + 1) * log_q
#     # print("FHE size", (m * (256 + 2 * secpar) * (n + 1) * log_q) / 8 / 10**9)

#     # input_ext = 1 + input + m * (256) * (n + 1) * log_q
#     input_ext = 1 + input + 256
#     # 2252 + 6000
#     # c_att
#     print("poly size", log_q * n / 8 / 10**6)
#     c_att_size = log_q * n * input_ext * m
#     size += c_att_size
#     print("c_att", c_att_size / 8 / 10**9)
#     # c_t
#     c_t_size = log_q * n * (1) * m
#     size += c_t_size
#     print("c_t", c_t_size / 8 / 10**9)

#     # p
#     p_size = log_q * n * m_b
#     size += p_size
#     print("p", p_size / 8 / 10**9)
#     # M/N size
#     m_n_size = math.log(bound_b) * n * m_b * m_b
#     print("m_n_size", m_n_size / 8 / 10**9)
#     print("m_n_size total", m_n_size * 4 / 8 / 10**9)
#     size += 2 * 2 * m_n_size
#     # K size
#     # k_size = math.log(bound_b) * m_b * (input_ext + n + 1) * m_s
#     # print("log_q", log_q)
#     # print("(2 + log_q)", m_b)
#     # print("(input_ext + n + 1)", (input_ext + n + 1))
#     # print("n_s", n_s)
#     print("log bound_b", math.log(bound_b))
#     print("n", n)
#     print("m_b", m_b)
#     print("input_ext", input_ext)
#     # print("(2 * input_ext + log_q)", (2 * input_ext + log_q))
#     k_size = math.log(bound_b) * n * m_b * (input_ext + 1) * m
#     # (input_ext + 1) * m
#     # log_q * (2 + log_q) * 2 * (input_ext + n + 1) * log_q * n_s
#     print("k_size", k_size / 8 / 10**9)
#     print("k_size total", k_size * 2 * input / 8 / 10**9)
#     size += input * 2 * k_size
#     # K_f
#     k_f_size = math.log(bound_b) * n * m_b * output
#     size += k_f_size
#     print("K_f", k_f_size / 8 / 10**9)
#     return size / 8 / (10**9)


if __name__ == "__main__":
    secpar = 100
    n = 2**10
    d = 14
    # alpha = 2 ** (-320)
    input_size = 1
    m_polys = [[0, 2720, 2176, 17408]]
    q, stddev_e, p, estimated_secpar, alpha = find_params(
        secpar, n, d, input_size, m_polys
    )
    print(f"q: {q}, log_2 q: {math.log2(q)}")
    print(f"stddev_e: {stddev_e}")
    print(f"p: {p}, log_2 p: {math.log2(p)}")
    print(f"estimated_secpar: {estimated_secpar}")

    # Log parameters to params.log file
    log_params_to_file(
        secpar, n, d, alpha, input_size, m_polys, q, stddev_e, p, estimated_secpar
    )
    # output_secpar(586, 2**32, 2 ** (-24.8) * 2**32)
    # q = 2**1024
    # alpha_log = 2 ** (-30)
    # input = 4
    # output = 1
    # depth = 1
    # n = 2**14
    # n_t = 2**14
    # q = 2**360
    # alpha = 2 ** (-350)
    # target_secpar = 80
    # params = derive_auto_params(n, n_t, q, math.ceil(alpha * q))
    # print(params)
    # error_bound = estimate_noise_norm(params, input, output, depth)
    # # print(error_bound)
    # print(math.log2(error_bound))
    # print(params["log_q"] - 2 - params["secpar"])
    # print("is less than q:", error_bound < params["q"])
    # print(
    #     "is less than (q/4)/(2^secpar):",
    #     math.log2(error_bound) < params["log_q"] - 2 - params["secpar"],
    # )
    # # print("is sigma_e < p", params["sigma_e"] < params["p"])
    # print("secpar", params["secpar"])
    # obf_size_gb = estimate_obf_size(params, input, output, depth)
    # print("obf_size_mb", obf_size_gb, "[GB]")


# input=1, output=1, depth=2, q = 2**1024, alpha = 2**(-200), n = 2**13, secpar = 105, log_e = 986
# input=1, output=1, depth=2, q = 2**256, alpha = 2**(-200), n = 2**13, secpar = 105, log_e = 206.36068319307367
# input=2, output=1, depth=2, q = 2**256, alpha = 2**(-200), n = 2**13, secpar = 105, log_e = 252.29242579686007
# input=4, output=1, depth=2, q = 2**512, alpha = 2**(-400), n = 2**14, secpar = 105, log_e = 490.20420897474514
# input=8, output=1, depth=2, q = 2**1024, alpha = 2**(-800), n = 2**15, secpar = 106, log_e = 987.5056715963805
# input=16, output=1, depth=2, q = 2**2048, alpha = 2**(-1600), n = 2**16, secpar = 106, log_e = 2027.6967419442458
# input=32, output=1, depth=2, q = 2**4096, alpha = 2**(-3320), n = 2**17, secpar = 101, log_e = 4080.1957810365197
# input=64, output=1, depth=2, q = 2**8196, alpha = 2**(-6950), n = 2**18, secpar = 94, log_e = 8186.687219336921
# input=128, output=1, depth=2, q = 2**16392, alpha = 2**(-14600), n = 2**19, secpar = 87, log_e = 16384.15472322582
# input=2, output=1, depth=4, q = 2**256, alpha = 2**(-225), n = 2**13, secpar = 88, log_e = 255.76379431529352
#
# + 2 * secpar
# input=2, output=1, depth=4, q = 2**512, alpha = 2**(-400), n = 2**14, secpar = 105, log_e = 362.4993302188647
# input=4, output=1, depth=4, q = 2**1024, alpha = 2**(-800), n = 2**15, secpar = 106, log_e = 629.7919250697457
# input=8, output=1, depth=4, q = 2**2048, alpha = 2**(-1600), n = 2**16, secpar = 106, log_e = 1263.4141811117784
# input=16, output=1, depth=4, q = 2**4096, alpha = 2**(-3200), n = 2**17, secpar = 106, log_e = 2576.3387780273965
# input=32, output=1, depth=4, q = 2**8192, alpha = 2**(-6400), n = 2**18, secpar = 106, log_e = 5296.521833879402
# input=64, output=1, depth=4, q = 2**16384, alpha = 2**(-12800), n = 2**19, secpar = 106, log_e = 10928.456653016947
# input=128, output=1, depth=4, q = 2**32768, alpha = 2**(-25600), n = 2**20, secpar = 106, log_e = 22578.237454920152
# input=256, output=1, depth=4, q = 2**65536, alpha = 2**(-51200), n = 2**21, secpar = 106, log_e = 46652.16947979474

# input=16, output=1, depth=32, q = 2**4096, alpha = 2**(-3200), n = 2**17, secpar = 106, log_e = 2744.914181351132
# input=16, output=1, depth=128, q = 2**4096, alpha = 2**(-3200), n = 2**17, secpar = 106, log_e = 5912.915237886846

# input=2, output=1, depth=2, 2**13, 2**13, 2 + 512, 2**512, 2 ** (512 - 263), 2 ** (512 - 500), secpar = 70, size = 94068.46295686903 GB

# Jan 30, 2025
# input=1, output=1, depth=1, n=2**13, q=2**256, p=2**16, sigma_e=2 ** (256 - 245), sigma_b=3.3, secpar=78, size = 802.0288007854358 GB
# Uniform secret
# input=1, output=1, depth=1, n=2**12, q=2**140, p=2**6, sigma_e=2 ** (140 - 135), sigma_b=0.01, secpar=72, size = 41.144327554 GB
# input=2, output=1, depth=1, n=2**12, q=2**140, p=2**6, sigma_e=2 ** (140 - 135), sigma_b=0.01, secpar=72, size = 41.154365204 GB
# input=4, output=1, depth=1, n=2**12, q=2**140, p=2**6, sigma_e=2 ** (140 - 135), sigma_b=0.01, secpar=72, size = 41.174440504 GB
# input=8, output=1, depth=1, n=2**12, q=2**140, p=2**6, sigma_e=2 ** (140 - 135), sigma_b=0.01, secpar=72, size = 41.214591104 GB
# input=16, output=1, depth=1, n=2**12, q=2**140, p=2**6, sigma_e=2 ** (140 - 135), sigma_b=0.01, secpar=72, size = 41.294892304 GB
# input=32, output=1, depth=1, n=2**12, q=2**140, p=2**6, sigma_e=2 ** (140 - 135), sigma_b=0.01, secpar=72, size = 41.455494704 GB
# input=64, output=1, depth=1, n=2**12, q=2**140, p=2**6, sigma_e=2 ** (140 - 135), sigma_b=0.01, secpar=72, size = 41.776699504 GB
# input=128, output=1, depth=1, n=2**12, q=2**140, p=2**6, sigma_e=2 ** (140 - 135), sigma_b=0.01, secpar=72, size = 42.419109104 GB
# The above results seem to invalid because the bound_b is just one.
# Jan 31, 2025
# m = inf in output_secpar
# input=1, output=1, depth=1, n=2**12, q=2**160, p=2**22, sigma_e=2 ** (160 - 139), sigma_b=0.2, secpar=81, size = 64.91435137376314 GB
# input=2, output=1, depth=1, n=2**12, q=2**161, p=2**22, sigma_e=2 ** (161 - 140), sigma_b=0.15, secpar=80, size = 76.92146251976315 GB
# input=4, output=1, depth=1, n=2**12, q=2**163, p=2**23, sigma_e=2 ** (163 - 141), sigma_b=0.15, secpar=80, size = 103.28398017749898 GB

# Feb 5, 2025
