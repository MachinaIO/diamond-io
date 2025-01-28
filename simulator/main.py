#!/usr/bin/env sage -python
# import sys
# from sage.all import *
from estimator.estimator import *
from estimator.estimator.lwe_parameters import *
from estimator.estimator.nd import *
import math
from decimal import Decimal, getcontext

getcontext().prec = 50


def output_secpar(n: int, q: int, stddev: int):
    params = LWEParameters(n, q, Binary, DiscreteGaussian(stddev))
    estim = LWE.estimate.rough(params)
    print(estim)
    min_secpar = math.log2(min(val["rop"] for val in estim.values()))
    print(min_secpar)
    return min_secpar


def bound_from_stddev(stddev: int, secpar: int):
    return stddev * math.floor(math.sqrt(secpar))


def derive_auto_params(n: int, n_s: int, q: int, sigma_e: int):
    log_q = math.floor(math.log2(q))
    m = (n + 1) * log_q
    m_s = (n_s + 1) * log_q
    n_b = 2 * (n_s + 1)
    m_b = n_b * log_q
    sigma_b = math.floor(math.sqrt(n_b * log_q * math.log2(m_b)))
    secpar_n = math.floor(output_secpar(n, q, sigma_e))
    print("secpar_n:", secpar_n)
    secpar_n_s = output_secpar(n_s, q, sigma_e)
    print("secpar_n_s:", secpar_n_s)
    secpar = min(secpar_n, secpar_n_s)
    print("secpar:", secpar)
    bound_e = bound_from_stddev(sigma_e, secpar)
    print("bound_e:", bound_e)
    bound_b = bound_from_stddev(sigma_b, secpar)
    print("bound_b:", bound_b)
    # b_c_frac = math.floor(
    #     (
    #         sigma_e
    #         * math.floor((m_s * sigma_b) ** 2)
    #         * math.floor(secpar**2)
    #         * math.floor(math.sqrt(secpar))
    #     )
    #     / (math.floor((m_s * sigma_b) ** 2) * secpar - m_s)
    # )
    b_c_frac = sigma_e * math.floor(math.sqrt(secpar))
    return {
        "q": q,
        "log_q": log_q,
        "n": n,
        "m": m,
        "n_s": n_s,
        "m_s": m_s,
        "n_b": n_b,
        "m_b": m_b,
        "secpar": secpar,
        "sigma_e": sigma_e,
        "sigma_b": sigma_b,
        "bound_e": bound_e,
        "bound_b": bound_b,
        "b_c_frac": b_c_frac,
    }


def estimate_noise_norm(params, input: int, output: int, depth: int):
    q = params["q"]
    log_q = params["log_q"]
    n = params["n"]
    m = params["m"]
    n_s = params["n_s"]
    m_s = params["m_s"]
    n_b = params["n_b"]
    m_b = params["m_b"]
    secpar = params["secpar"]
    sigma_e = params["sigma_e"]
    sigma_b = params["sigma_b"]
    bound_e = params["bound_e"]
    bound_b = params["bound_b"]
    b_c_frac = params["b_c_frac"]
    b_c = (bound_e - b_c_frac) * (m_s**input) + b_c_frac * (
        (sigma_b * sigma_b * secpar) ** input
    )
    print("b_c", b_c)
    input_ext = 1 + input + m * (256 + 2 * secpar) * (n + 1) * log_q
    # b_f_exp1 = math.floor(
    #     depth * math.floor(math.log2(m_s)) * math.floor(math.log2(log_q))
    #     + (math.floor(math.log2(log_q)) ** 2)
    #     + 2
    # )
    b_f_exp1 = depth
    print("b_f_exp1", b_f_exp1)
    b_f_term1 = (
        b_c * (input_ext + n + 2) * ((n + 1) * log_q) * (((m_s + 2) ** b_f_exp1))
    )
    print("b_f_term1", b_f_term1)
    b_f_term2 = bound_e * log_q * ((m + 2) ** (depth + 1))
    print("b_f_term2", b_f_term2)
    b_f = b_f_term1 + b_f_term2
    print("b_f", b_f)
    b_z = b_f + sigma_e * ((m_s * sigma_b) ** (2 * input + 1)) * (secpar ** (input + 1))

    # print("b_z", b_z)
    # print("log2_b_z", math.log2(b_z))
    return b_z
    # b_f =
    # b_c = (params.bound_b - )


if __name__ == "__main__":
    # output_secpar(586, 2**32, 2 ** (-24.8) * 2**32)
    q_log = 500
    alpha_log = -2
    params = derive_auto_params(
        2**7, 2**7, 2**q_log, math.floor(2 ** (q_log + alpha_log))
    )
    print(params)
    error_bound = estimate_noise_norm(params, 2, 1, 2)
    # print(error_bound)
    print(math.log2(error_bound))
    print("is less than q:", error_bound < params["q"])
