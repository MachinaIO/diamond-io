#!/usr/bin/env sage -python
# import sys
# from sage.all import *
from estimator.estimator import *
from estimator.estimator.lwe_parameters import *
from estimator.estimator.nd import *
import math
from decimal import Decimal, getcontext

getcontext().prec = 100


def output_secpar(n: int, q: int, s_dist: NoiseDistribution, stddev: int):
    params = LWEParameters(n, q, s_dist, DiscreteGaussian(stddev))
    estim = LWE.estimate.rough(params)
    # print(estim)
    min_secpar = math.log2(min(val["rop"] for val in estim.values()))
    # print(min_secpar)
    return min_secpar


def bound_from_stddev(stddev: int, secpar: int):
    return math.ceil(stddev * math.ceil(math.sqrt(secpar)))


def derive_auto_params(n_t: int, n_s: int, q: int, p: int, sigma_e: int, sigma_b: int):
    log_q = math.ceil(math.log2(q))
    print("log_q:", log_q)
    m_t = (n_t + 1) * log_q
    m_s = (n_s + 1) * log_q
    print("m_s:", m_s)
    log_p = math.ceil(math.log2(p))
    n_b = 2 * (n_s + 1)
    m_b = 2 + log_p
    # m_b = n_b * log_p
    # sigma_b = math.ceil(2 * math.sqrt(n * log_p))
    secpar_s = math.ceil(output_secpar(n_s, q, UniformMod(q), sigma_e))
    print("secpar_n:", secpar_s)
    # secpar_t = math.ceil(output_secpar(n_t, q, Binary, sigma_e))
    # print("secpar_t", secpar_t)

    # secpar_n_s = math.ceil(output_secpar(n_s, m_s, q, sigma_e))
    # print("secpar_n_s:", secpar_n_s)
    secpar_b = output_secpar(n_b, p, UniformMod(q), sigma_b)
    print("secpar_b:", secpar_b)
    secpar = min(secpar_s, secpar_b)
    print("secpar:", secpar)
    bound_e = bound_from_stddev(sigma_e, secpar)
    print("bound_e:", bound_e)
    if bound_e <= 1:
        raise ValueError("bound_e should be larger than 1")
    bound_b = bound_from_stddev(sigma_b, secpar)
    print("bound_b:", bound_b)
    if bound_b <= 1:
        raise ValueError("bound_b should be larger than 1")
    # b_c_frac = math.ceil(
    #     (
    #         Decimal(sigma_e)
    #         * math.ceil((m_b * sigma_b) ** 2)
    #         * math.ceil(secpar**2)
    #         * math.ceil(math.sqrt(secpar))
    #     )
    #     / Decimal(math.ceil((m_b * sigma_b) ** 2) * secpar - m_s)
    # )
    b_c_frac = sigma_e * math.ceil(math.sqrt(secpar))
    return {
        "q": q,
        "log_q": log_q,
        "n": n_t,
        "m": m_t,
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
        "p": p,
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
    print("log b_c_frac", math.log2(b_c_frac))
    print("log m_s", math.log2(m_s))
    print("log (sigma_b * sigma_b * secpar)", math.log2((sigma_b * sigma_b * secpar)))
    # print(
    #     "log (sigma_b * sigma_b * secpar) ** input",
    #     math.ceil(math.log2((sigma_b * sigma_b * secpar) ** input)),
    # )
    b_c = (bound_e - b_c_frac) * (m_s**input) + b_c_frac * math.ceil(
        (sigma_b * sigma_b * secpar) ** input
    )
    print("b_c", b_c)
    print("log b_c", math.log2(b_c) if b_c > 0 else 0)
    input_ext = 1 + input
    # 1 + input + m * (256 + 2 * secpar) * (n + 1) * log_q
    # b_f_exp1 = math.ceil(
    #     depth * math.ceil(math.log2(m_s)) * math.ceil(math.log2(log_q))
    #     + (math.ceil(math.log2(log_q)) ** 2)
    #     + 2
    # )
    b_f_exp1 = depth
    print("b_f_exp1", b_f_exp1)
    print("log (input_ext + n + 2)", math.log2((input_ext + n + 2)))
    print("log ((n + 1) * log_q)", math.log2((n + 1) * log_q))
    print(
        "log math.ceil((m_s + 2) ** b_f_exp1)",
        math.log2(math.ceil((m_s + 2) ** b_f_exp1)),
    )
    b_f_term1 = (
        b_c * (input_ext + n + 2) * ((n + 1) * log_q) * math.ceil((m_s + 2) ** b_f_exp1)
    )
    print("log b_f_term1", math.log2(b_f_term1))
    b_f_term2 = bound_e * log_q * ((m + 2) ** (depth + 1))
    print("log b_f_term2", math.log2(b_f_term2))
    b_f = b_f_term1 + b_f_term2
    print("log b_f", math.log2(b_f))
    b_z = b_f + sigma_e * ((m_b * sigma_b) ** (2 * input + 1)) * (secpar ** (input + 1))
    print("log b_z", math.log2(b_z))
    print("log b_z", math.log2(b_z))
    # print("b_z", b_z)
    # print("log2_b_z", math.log2(b_z))
    # print("b_z + 2^**(2*secpar)", b_z * (2 ** (2 * secpar)))
    return b_z
    # b_f =
    # b_c = (params.bound_b - )


def estimate_obf_size(params, input: int, output: int, depth: int):
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

    # bits
    size = 256
    # h
    size += 256
    # FHE encryption
    # size += m * (256 + 2 * secpar) * (n + 1) * log_q
    # print("FHE size", (m * (256 + 2 * secpar) * (n + 1) * log_q) / 8 / 10**9)

    # input_ext = 1 + input + m * (256) * (n + 1) * log_q
    input_ext = 1 + input
    # 2252 + 6000
    # c_att
    size += log_q * input_ext * m_s
    print("c_att", (log_q * input_ext * m_s) / 8 / 10**9)
    # c_t
    size += log_q * (n + 1) * m_s
    print("c_t", (log_q * (n + 1) * m_s) / 8 / 10**9)

    # p
    size += log_q * m_b
    print("p", (log_q * m_b) / 8 / 10**9)
    # M/N size
    m_n_size = math.log(bound_b) * m_b * m_b
    print("m_n_size", m_n_size / 8 / 10**9)
    size += 2 * input * 2 * m_n_size
    # K size
    k_size = math.log(bound_b) * m_b * (input_ext + n + 1) * m_s
    print("k_size", k_size / 8 / 10**9)
    size += input * 2 * k_size
    # K_f
    size += math.log(bound_b) * m_b * output
    print("K_f", (math.log(bound_b) * m_b * output) / 8 / 10**9)
    return size / 8 / (10**9)


if __name__ == "__main__":
    # output_secpar(586, 2**32, 2 ** (-24.8) * 2**32)
    # q = 2**1024
    # alpha_log = 2 ** (-30)
    input = 1
    output = 1
    depth = 1
    params = derive_auto_params(2**12, 2**12, 2**160, 2**160, 2 ** (160 - 139), 0.15)
    print(params)
    error_bound = estimate_noise_norm(params, input, output, depth)
    # print(error_bound)
    print(math.log2(error_bound))
    print(math.log2(params["log_q"] - 2 - params["secpar"]))
    print("is less than q:", error_bound < params["q"])
    print(
        "is less than (q/4)/(2^secpar):",
        math.log2(error_bound) < params["log_q"] - 2 - params["secpar"],
    )
    print("is sigma_e < p", params["sigma_e"] < params["p"])
    print("secpar", params["secpar"])
    obf_size_gb = estimate_obf_size(params, input, output, depth)
    print("obf_size_mb", obf_size_gb, "[GB]")


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
