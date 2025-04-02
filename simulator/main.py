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
from norms import CircuitNorms

getcontext().prec = 100


def log_params_to_file(
    secpar: int,
    n: int,
    d: int,
    input_size: int,
    norms_path: str,
    q: int,
    stddev_e_encoding: int,
    stddev_e_hardcode: int,
    stddev_e_p: int,
    p: int,
    estimated_secpar: float,
    size: int,
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
    # m_polys_str = str(m_polys).replace(" ", "")

    # Create log entry with key information
    log_entry = (
        f"{current_date}, "
        f"secpar={secpar}, "
        f"n={n}, "
        f"d={d}, "
        f"input_size={input_size}, "
        f"norms_path={norms_path}, "
        f"q={q}, "
        f"log_q={log_q}, "
        f"stddev_e_encoding={stddev_e_encoding}, "
        f"stddev_e_hardcode={stddev_e_hardcode}, "
        f"stddev_e_p={stddev_e_p}, "
        f"p={p}, "
        f"log_p={log_p}, "
        f"estimated_secpar={estimated_secpar}, "
        f"size={size} [GB]\n"
    )

    # Append to params.log file
    with open("params.log", "a") as f:
        f.write(log_entry)

    print(f"Parameters logged to params.log")


def find_params(
    target_secpar: int,
    n: int,
    d: int,
    max_log_q: int,
    input_size: int,
    norms_path: str,
):

    min_alpha_ks = [-2000 for _ in range(3)]
    max_alpha_ks = [-1 for _ in range(3)]
    min_q_k = target_secpar + 2
    max_q_k = max_log_q
    middle_q_k = math.floor((min_q_k + max_q_k) // 2)
    circuit_norms = CircuitNorms.load_from_file(norms_path, max_log_q)
    found_alphas = [[] for _ in range(3)]
    # Binary search for alphas
    for i in range(3):
        q = 2 ** (middle_q_k)
        total_n = 0
        dist = Binary
        if i == 0:
            # encoding sigma
            total_n = n * (d + 1)
        elif i == 1:
            # hardcoded key sigma
            total_n = n
            dist = UniformMod(q)
        else:
            # p sigma
            total_n = 2 * (n * (d + 1))
        while min_alpha_ks[i] + 1 < max_alpha_ks[i]:
            min_alpha_k = min_alpha_ks[i]
            max_alpha_k = max_alpha_ks[i]
            alpha_k = (min_alpha_k + max_alpha_k) / 2
            print(f"min_alpha_k: {min_alpha_k}")
            print(f"max_alpha_k: {max_alpha_k}")
            print(f"alpha_k: {alpha_k}")
            stddev_e = 2 ** (middle_q_k + alpha_k)
            estimated_secpar = estimate_secpar(total_n, q, dist, stddev_e)
            print("target_secpar:", target_secpar)
            print("estimated_secpar:", estimated_secpar)
            if target_secpar > estimated_secpar:
                print(
                    f"target_secpar {target_secpar} > estimated_secpar {estimated_secpar}"
                )
                min_alpha_ks[i] = alpha_k
            else:
                found_alphas[i].append(alpha_k)
                print(f"found alpha_k: {alpha_k}")
                max_alpha_ks[i] = alpha_k
        if len(found_alphas[i]) == 0:
            raise ValueError(f"the {i}-th alpha is not found after binary search")
    alpha_encoding = 2 ** min(found_alphas[0])
    alpha_hardcode = 2 ** min(found_alphas[1])
    alpha_p = 2 ** min(found_alphas[2])
    # print(f"found alpha: {alpha}")
    print(f"found alpha_encoding: {alpha_encoding}")
    print(f"found alpha_hardcode: {alpha_hardcode}")
    print(f"found alpha_p: {alpha_p}")

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
        stddev_e_encoding = alpha_encoding * q
        stddev_e_hardcode = alpha_hardcode * q
        stddev_e_p = alpha_p * q
        estimated_secpar_encoding = estimate_secpar(
            (d + 1) * n, q, Binary, stddev_e_encoding
        )
        estimated_secpar_hardcode = estimate_secpar(
            n, q, UniformMod(q), stddev_e_hardcode
        )
        estimated_secpar_p = estimate_secpar(2 * ((d + 1) * n), q, Binary, stddev_e_p)
        min_estimated_secpar = min(
            estimated_secpar_encoding, estimated_secpar_hardcode, estimated_secpar_p
        )
        print("target_secpar:", target_secpar)
        print("estimated_secpar:", min_estimated_secpar)
        if target_secpar > min_estimated_secpar:
            print(
                f"target_secpar {target_secpar} > estimated_secpar {min_estimated_secpar}"
            )
            max_q_k = q_k
            continue
        try:
            p = find_p(
                target_secpar,
                n,
                q,
                d,
                stddev_e_encoding,
                stddev_e_hardcode,
                stddev_e_p,
                input_size,
                circuit_norms,
            )
            print(f"found p: {p}")
            max_q_k = q_k
            found_params.append(
                (
                    q,
                    stddev_e_encoding,
                    stddev_e_hardcode,
                    stddev_e_p,
                    p,
                    min_estimated_secpar,
                )
            )
        except ValueError as e:
            print(f"ValueError: {e}")
            min_q_k = q_k
    if found_params == []:
        raise ValueError("p is not found after binary search")
    # minimum q in found_params
    (q, stddev_e_encoding, stddev_e_hardcode, stddev_e_p, p, estimated_secpar) = min(
        found_params, key=lambda x: x[0]
    )
    log_q = math.ceil(math.log2(q))
    size = compute_obf_size(
        n,
        log_q,
        d,
        (d + 1) * log_q,
        2 * (d + 1) * (log_q + 2),
        input_size,
        1,
    )
    return (
        q,
        stddev_e_encoding,
        stddev_e_hardcode,
        stddev_e_p,
        p,
        estimated_secpar,
        size,
    )


def find_p(
    secpar: int,
    n: int,
    q: int,
    d: int,
    stddev_e_encoding: int,
    stddev_e_hardcode: int,
    stddev_e_p: int,
    input_size: int,
    circuit_norms: CircuitNorms,
):
    log_q = math.ceil(math.log2(q))
    norm_b = compute_norm_b(n, log_q, d)
    final_err = bound_final_error(
        secpar,
        n,
        log_q,
        d,
        stddev_e_encoding,
        stddev_e_hardcode,
        stddev_e_p,
        norm_b,
        input_size,
        circuit_norms,
    )

    # Convert final_err to Decimal for high precision
    final_err_decimal = Decimal(str(final_err))

    # Calculate log2 using Decimal
    # Handle infinity or very large numbers
    if math.isinf(final_err):
        raise ValueError(f"Error: final_err is infinity.")
    elif final_err_decimal > 0:
        try:
            # Try to calculate log2 using Decimal
            log_final_err = math.ceil(
                float(final_err_decimal.ln() / Decimal("0.693147180559945"))
            )  # ln(2)
        except (OverflowError, ValueError):
            raise ValueError(f"Error: final_err is too large for ln().")
    else:
        raise ValueError(f"Cannot calculate log2 of non-positive value: {final_err}")

    if log_q - 2 < log_final_err + secpar:
        raise ValueError(
            f"log_q - 2 >= log_final_err + secpar should hold. log2(final_error): {log_final_err}, log2(q): {log_q})"
        )

    # Use Decimal for high precision arithmetic
    # Convert to Decimal for high precision calculations
    q_decimal = Decimal(q)
    p_decimal = (q_decimal / Decimal("4")).to_integral_exact(
        rounding="ROUND_CEILING"
    ) - final_err_decimal
    p = int(p_decimal)
    if p < 0:
        raise ValueError(f"p should be non-negative: {p}")

    # Calculate log2(p) using Decimal
    # if p_decimal > 0:
    #     log_p = math.ceil(float(p_decimal.ln() / Decimal("0.693147180559945")))  # ln(2)
    # else:
    #     raise ValueError(f"Cannot calculate log2 of non-positive value: {p}")

    # if log_p - log_final_err < secpar:
    #     raise ValueError(
    #         f"p - error should be larger than 2^secpar (given p: {p}, secpar: {secpar})"
    #     )

    return p


def compute_norm_b(
    n: int,
    log_q: int,
    d: int,
):
    c_0 = 1.8
    c_1 = 4.7
    sigma = 4.578
    return (
        c_0
        * sigma
        * (3 * sigma)
        * (sqrt_ceil(2 * (d + 1) * n * log_q) + sqrt_ceil(2 * n) + c_1)
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
    stddev_e_encoding: int,
    stddev_e_hardcode: int,
    stddev_e_p: int,
    norm_b: int,
    input_size: int,
    circuit_norms: CircuitNorms,
):
    # Convert all inputs to Decimal for high precision
    # secpar_d = Decimal(secpar)
    n_d = Decimal(n)
    log_q_d = Decimal(log_q)
    d_d = Decimal(d)
    stddev_e_encoding_d = Decimal(stddev_e_encoding)
    stddev_e_hardcode_d = Decimal(stddev_e_hardcode)
    stddev_e_p_d = Decimal(stddev_e_p)
    b_norm_d = Decimal(norm_b)
    print(f"b_norm_d: {b_norm_d}")
    m_d = (d_d + Decimal(1)) * log_q_d
    h_norms_d = [Decimal(x) for x in circuit_norms.compute_norms(m_d)]
    print(f"h_norms_d: {h_norms_d}")
    h_norm_sum_d = sum(h_norms_d)
    # Calculate intermediate values with Decimal
    m_b_d = Decimal(2) * (d_d + Decimal(1)) * (log_q_d + Decimal(2))
    sqrt_secpar_d = Decimal(sqrt_ceil(secpar))

    # Use Decimal for all calculations to maintain precision
    scale_coeff_d = b_norm_d ** Decimal(2)
    bound_p_d = stddev_e_p_d * sqrt_secpar_d
    bound_c_d = stddev_e_encoding_d * sqrt_secpar_d

    for _ in range(input_size):
        bound_v_d = (scale_coeff_d ** Decimal(2)) * bound_p_d
        bound_c_d = n_d * m_d * bound_c_d + bound_v_d
        bound_p_d = bound_v_d

    # Evaluate each polynomial in m_polys at the value of m using Decimal
    # evaluated_polys_d = []
    # for poly in m_polys:
    #     # Evaluate polynomial: sum(coeff * m^i for i, coeff in enumerate(poly))
    #     result_d = Decimal(0)
    #     for i, coeff in enumerate(poly):
    #         result_d += Decimal(coeff) * (m_d ** Decimal(i))
    #     evaluated_polys_d.append(result_d)

    # # Find max value using Decimal
    # if evaluated_polys_d:
    #     max_evaluated_poly_d = max(evaluated_polys_d)
    # else:
    #     max_evaluated_poly_d = Decimal(1)  # Default if no polynomials

    bound_c_final_d = bound_c_d * h_norm_sum_d
    bound_v_final_d = bound_p_d * b_norm_d

    # Return the final result as a Decimal
    return bound_c_final_d + bound_v_final_d + stddev_e_hardcode_d * sqrt_secpar_d


def compute_obf_size(
    n: int,
    log_q: int,
    d: int,
    m: int,
    m_b: int,
    input_size: int,
    output_size: int,
):
    size = 32
    encoding_init_size = log_q * n * input_size * m
    print("encoding_init_size GB", encoding_init_size / 8 / 10**9)
    size += encoding_init_size
    p_init_size = log_q * n * m_b
    print("p_init_size GB", p_init_size / 8 / 10**9)
    size += p_init_size
    b_norm = Decimal(compute_norm_b(n, log_q, d))
    bound_b_log = math.ceil(math.log2(b_norm))
    m_n_preimages_size = 2 * input_size * bound_b_log * n * m_b * m_b
    print("m_n_preimages_size GB", m_n_preimages_size / 8 / 10**9)
    print(
        "indiv m_n_preimages_size GB", m_n_preimages_size / 8 / 10**9 / 2 / input_size
    )
    size += m_n_preimages_size
    k_preimages_size = input_size * bound_b_log * n * m_b * (input_size * m)
    print("k_preimage_size GB", k_preimages_size / 8 / 10**9)
    print("indiv k_preimage_size GB", k_preimages_size / 8 / 10**9 / input_size)
    size += k_preimages_size
    packed_output_size = math.ceil(output_size / n)
    final_preimage_size = bound_b_log * n * m_b * packed_output_size
    print("final_preimage_size GB", final_preimage_size / 8 / 10**9)
    size += final_preimage_size
    return size / 8 / 10**9


def sqrt_ceil(x):
    return math.ceil(math.sqrt(x))


if __name__ == "__main__":
    secpar = 80
    n = 2**13
    d = 2
    max_log_q = 612
    # alpha = 2 ** (-320)
    input_size = 1
    # m_polys = [[0, 100, 200, 2000]]
    norms_path = "final_bits_norm_n_13_q_612.json"
    (
        q,
        stddev_e_encoding,
        stddev_e_hardcode,
        stddev_e_p,
        p,
        estimated_secpar,
        size,
    ) = find_params(secpar, n, d, max_log_q, input_size, norms_path)

    print(f"q: {q}, log_2 q: {math.log2(q)}")
    print(f"stddev_e_encoding: {stddev_e_encoding}")
    print(f"stddev_e_hardcode: {stddev_e_hardcode}")
    print(f"stddev_e_p: {stddev_e_p}")
    print(f"p: {p}, log_2 p: {math.log2(p)}")
    print(f"estimated_secpar: {estimated_secpar}")
    print(f"size: {size} [GB]")
    # Log parameters to params.log file
    log_params_to_file(
        secpar,
        n,
        d,
        input_size,
        norms_path,
        q,
        stddev_e_encoding,
        stddev_e_hardcode,
        stddev_e_p,
        p,
        estimated_secpar,
        size,
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
