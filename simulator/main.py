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

getcontext().prec = 300


def log_params_to_file(
    secpar: int,
    n: int,
    d: int,
    base_bits: int,
    input_size: int,
    input_width: int,
    norms_path: str,
    crt_bits: int,
    crt_depth: int,
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
    # Calculate log_q, q, and log_p
    log_q = crt_bits * crt_depth
    q = 2 ** (log_q + 1) - 1
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
        f"crt_bits={crt_bits}, "
        f"crt_depth={crt_depth}, "
        f"base_bits={base_bits}, "
        f"input_size={input_size}, "
        f"input_width={input_width}, "
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
    base_bits: int,
    crt_bits: int,
    max_crt_depth: int,
    input_size: int,
    input_width: int,
    norms_path: str,
):
    # crt_bits * depth >= target_secpar+2 => depth >= (target_secpar+2) / crt_bits
    min_crt_depth = math.ceil((target_secpar + 2) / crt_bits)
    max_log_base_q = math.ceil(crt_bits / base_bits) * max_crt_depth
    print(f"max_log_base_q: {max_log_base_q}")
    circuit_norms = CircuitNorms.load_from_file(norms_path, max_log_base_q)
    found_params = []
    iters = 0
    while min_crt_depth + 1 < max_crt_depth and iters < 100:
        iters += 1
        crt_depth = math.floor((min_crt_depth + max_crt_depth) // 2)
        q = 2 ** (crt_bits * crt_depth + 1) - 1
        print(f"min_crt_depth: {min_crt_depth}")
        print(f"max_crt_depth: {max_crt_depth}")
        print(f"crt_depth: {crt_depth}")
        print(f"q: {q}")

        min_alpha_ks = []
        for i in range(3):
            found_alpha_ks = []
            dist = Binary
            min_alpha_k = -crt_bits * crt_depth + 2
            max_alpha_k = -1
            if i == 0:
                # encoding sigma
                total_n = n * (d + 1)
            elif i == 1:
                # hardcoded key sigma
                total_n = n
                # dist = UniformMod(q)
            else:
                # p sigma
                total_n = 2 * (n * (d + 1))
            while min_alpha_k + 1 < max_alpha_k:
                alpha_k = (min_alpha_k + max_alpha_k) / 2
                print(f"min_alpha_k: {min_alpha_k}")
                print(f"max_alpha_k: {max_alpha_k}")
                print(f"alpha_k: {alpha_k}")
                stddev_e = Decimal(2 ** Decimal(crt_bits * crt_depth + alpha_k))
                estimated_secpar = estimate_secpar(total_n, q, dist, stddev_e)
                print("target_secpar:", target_secpar)
                print("estimated_secpar:", estimated_secpar)
                if target_secpar > estimated_secpar:
                    print(
                        f"target_secpar {target_secpar} > estimated_secpar {estimated_secpar}"
                    )
                    min_alpha_k = alpha_k
                else:
                    found_alpha_ks.append(alpha_k)
                    print(f"found alpha_k: {alpha_k}")
                    max_alpha_k = alpha_k
                # raise ValueError(f"the {i}-th alpha is not found after binary search")
            if len(found_alpha_ks) == 0:
                continue
            min_alpha_ks.append(min(found_alpha_ks))
        if len(min_alpha_ks) < 3:
            print("not enough alpha_ks")
            max_crt_depth = crt_depth
            continue
        alpha_encoding_k = Decimal(min_alpha_ks[0])
        alpha_hardcode_k = Decimal(min_alpha_ks[1])
        alpha_p_k = Decimal(min_alpha_ks[2])
        print(f"found alpha_encoding_k: {alpha_encoding_k}")
        print(f"found alpha_hardcode_k: {alpha_hardcode_k}")
        print(f"found alpha_p_k: {alpha_p_k}")
        alpha_encoding = Decimal(2**alpha_encoding_k)
        alpha_hardcode = Decimal(2**alpha_hardcode_k)
        alpha_p = Decimal(2**alpha_p_k)
        # print(f"found alpha: {alpha}")
        print(f"found alpha_encoding: {alpha_encoding}")
        print(f"found alpha_hardcode: {alpha_hardcode}")
        print(f"found alpha_p: {alpha_p}")

        # if q_k + alpha_encoding_k < 1:
        #     print(f"q_k + alpha_encoding < 1")
        #     min_q_k = q_k
        #     continue
        # elif q_k + alpha_hardcode_k < 1:
        #     print(f"q_k + alpha_hardcode < 1")
        #     min_q_k = q_k
        #     continue
        # elif q_k + alpha_p_k < 1:
        #     print(f"q_k + alpha_p < 1")
        #     min_q_k = q_k
        #     continue
        stddev_e_encoding = alpha_encoding * Decimal(q)
        stddev_e_hardcode = alpha_hardcode * Decimal(q)
        stddev_e_p = alpha_p * Decimal(q)
        estimated_secpar_encoding = estimate_secpar(
            (d + 1) * n, q, Binary, stddev_e_encoding
        )
        estimated_secpar_hardcode = estimate_secpar(n, q, Binary, stddev_e_hardcode)
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
            max_crt_depth = crt_depth
            continue
        try:
            p = find_p(
                target_secpar,
                n,
                crt_bits,
                crt_depth,
                d,
                base_bits,
                stddev_e_encoding,
                stddev_e_hardcode,
                stddev_e_p,
                input_size,
                input_width,
                circuit_norms,
            )
            print(f"found p: {p}")
            max_crt_depth = crt_depth
            found_params.append(
                (
                    crt_depth,
                    stddev_e_encoding,
                    stddev_e_hardcode,
                    stddev_e_p,
                    p,
                    min_estimated_secpar,
                )
            )
        except ValueError as e:
            print(f"ValueError: {e}")
            min_crt_depth = crt_depth
    if found_params == []:
        raise ValueError("p is not found after binary search")
    # minimum q in found_params
    (
        crt_depth,
        stddev_e_encoding,
        stddev_e_hardcode,
        stddev_e_p,
        p,
        estimated_secpar,
    ) = min(found_params, key=lambda x: x[0])
    size = compute_obf_size(
        n,
        crt_bits,
        crt_depth,
        d,
        2**base_bits,
        input_size,
        input_width,
        1,
    )
    return (
        crt_depth,
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
    crt_bits: int,
    crt_depth: int,
    d: int,
    base_bits: int,
    stddev_e_encoding: int,
    stddev_e_hardcode: int,
    stddev_e_p: int,
    input_size: int,
    input_width: int,
    circuit_norms: CircuitNorms,
):
    log_q = crt_bits * crt_depth
    log_base_q = math.ceil(crt_bits / base_bits) * crt_depth
    base = 2**base_bits
    norm_b = compute_norm_b(n, log_base_q, d, base)
    (final_err, bound_s) = bound_final_error(
        secpar,
        n,
        log_base_q,
        d,
        base,
        stddev_e_encoding,
        stddev_e_hardcode,
        stddev_e_p,
        norm_b,
        input_size,
        input_width,
        circuit_norms,
    )

    # Convert final_err to Decimal for high precision
    # final_err_decimal = Decimal(str(final_err))
    # Calculate log2 using Decimal
    # Handle infinity or very large numbers
    if math.isinf(final_err):
        raise ValueError(f"Error: final_err is infinity.")
    elif final_err > 0:
        log_final_err = math.ceil(math.log2(float(final_err)))
    else:
        raise ValueError(f"Cannot calculate log2 of non-positive value: {final_err}")

    if log_q - 2 < log_final_err + secpar:
        raise ValueError(
            f"log_q - 2 >= log_final_err + secpar should hold. log2(final_error): {log_final_err}, log2(q): {log_q})"
        )
    print(f"final_err: {final_err}")
    print(f"log_final_err: {log_final_err}")
    # Use Decimal for high precision arithmetic
    # Convert to Decimal for high precision calculations
    prf_bound = 2 ** (log_q - 2) - (2 ** (log_final_err + 1) - 1)
    p = math.floor(prf_bound / bound_s / n / (d + 1))
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
    log_base_q: int,
    d: int,
    base: int,
):
    c_0 = 1.8
    c_1 = 4.7
    sigma = 4.578
    return (
        6.0
        * c_0
        * sigma
        * ((base + 1) * sigma)
        * (sqrt_ceil(2 * (d + 1) * n * log_base_q) + sqrt_ceil(2 * n) + c_1)
    )


def estimate_secpar(n: int, q: int, s_dist: NoiseDistribution, stddev: int):
    params = LWEParameters(n, q, s_dist, DiscreteGaussian(stddev))
    estim = LWE.estimate.rough(params)
    # print(estim)
    vals = estim.values()
    if len(vals) == 0:
        return 0
    min_rop_log = math.log2(min(val["rop"] for val in vals))
    print(f"min_rop_log: {min_rop_log}")
    if min_rop_log == float("inf"):
        return 100000
    min_secpar = math.ceil(min_rop_log)
    print(f"min_secpar: {min_secpar}")
    # print(min_secpar)
    return min_secpar


def bound_final_error(
    secpar: int,
    n: int,
    log_base_q: int,
    d: int,
    base: int,
    stddev_e_encoding: int,
    stddev_e_hardcode: int,
    stddev_e_p: int,
    norm_b: int,
    input_size: int,
    input_width: int,
    circuit_norms: CircuitNorms,
):
    # Convert all inputs to Decimal for high precision
    # secpar_d = Decimal(secpar)
    n = Decimal(n)
    n_sqrt = sqrt_ceil(n)
    log_base_q = Decimal(log_base_q)
    d = Decimal(d)
    stddev_e_encoding = Decimal(stddev_e_encoding)
    stddev_e_hardcode = Decimal(stddev_e_hardcode)
    stddev_e_p = Decimal(stddev_e_p)
    b_norm = Decimal(norm_b)
    print(f"b_norm_d: {b_norm}")
    m = (d + Decimal(1)) * log_base_q
    m_sqrt = sqrt_ceil(m)
    # [TODO] Support outputs larger than `log_t_q`
    h_norms = [Decimal(x) for x in circuit_norms.compute_norms(m_sqrt, n_sqrt, base)]
    print(f"h_norms: {h_norms}")
    h_norm_sum = sum(h_norms)
    # Calculate intermediate values with Decimal
    m_b = Decimal(2) * (d + Decimal(1)) * (log_base_q + Decimal(2))
    sqrt_secpar = Decimal(sqrt_ceil(secpar))
    base = Decimal(base)

    # Use Decimal for all calculations to maintain precision
    bound_p = stddev_e_p * sqrt_secpar
    print(f"stddev_e_p: {stddev_e_p}")
    print(f"sqrt_secpar: {sqrt_secpar}")
    print(f"stddev_e_encoding : {stddev_e_encoding}")
    bound_c = stddev_e_encoding * sqrt_secpar
    print(f"init bound_c: {bound_c}")
    if bound_c < 0:
        raise ValueError(f"bound_c should be non-negative: {bound_c}")
    bound_s = Decimal(1.0)

    input_depth = math.ceil(input_size / input_width)
    for _ in range(input_depth):
        bound_v = n * m_b * (b_norm ** Decimal(2)) * bound_p
        bound_c = n_sqrt * (base - 1) * bound_c * m_sqrt + bound_v
        # print(f"base_d: {base_d}")
        # print(f"m_d: {m_d}")
        # print(f"bound_c_d: {bound_c_d}")
        # print(f"base-dependent error: {n_d * (base_d-1) * m_d * bound_c_d }")
        bound_p = bound_v
        bound_s = bound_s * n * d

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
    packed_input_size = Decimal(math.ceil(input_size / n))
    bound_c_final = n_sqrt * bound_c * h_norm_sum * (packed_input_size * m)
    bound_v_final = n_sqrt * m_sqrt * b_norm * bound_p
    bound_rounding = bound_s
    print(f"bound_rounding: {bound_rounding}")
    print(f"log2(bound_rounding): {math.log2(bound_rounding)}")

    # Return the final result as a Decimal
    return (
        bound_c_final
        + bound_v_final
        + stddev_e_hardcode * sqrt_secpar
        + bound_rounding,
        bound_s,
    )


def compute_obf_size(
    n: int,
    crt_bits: int,
    crt_depth: int,
    d: int,
    base: int,
    input_size: int,
    input_width: int,
    output_size: int,
):
    size = 256
    packed_input_size = math.ceil(input_size / n)
    log_q = crt_bits * crt_depth
    log_base_q = math.ceil(crt_bits / base) * crt_depth
    m = (d + 1) * log_base_q
    print("m", m)
    m_b = 2 * (d + 1) * (log_base_q + 2)
    print("m_b", m_b)
    encoding_init_size = log_q * n * packed_input_size * m
    print("encoding_init_size GB", encoding_init_size / 8 / 10**9)
    size += encoding_init_size
    p_init_size = log_q * n * m_b
    print("p_init_size GB", p_init_size / 8 / 10**9)
    size += p_init_size
    b_norm = Decimal(compute_norm_b(n, log_base_q, d, base))
    bound_b_log = math.ceil(math.log2(b_norm))
    input_depth = math.ceil(input_size / input_width)
    m_n_preimages_size = (
        2 * (2**input_width) * input_depth * bound_b_log * n * m_b * m_b
    )
    print("m_n_preimages_size GB", m_n_preimages_size / 8 / 10**9)
    size += m_n_preimages_size
    k_preimages_size = (
        (2**input_width) * input_depth * bound_b_log * n * m_b * (packed_input_size * m)
    )
    print("k_preimage_size GB", k_preimages_size / 8 / 10**9)
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
    d = 1
    base_bits = 20
    crt_bits = 51
    max_crt_depth = 23
    input_size = 56
    input_width = 7
    if input_size % input_width != 0:
        raise ValueError("input_size should be divisible by input_width")
    norms_path = "final_bits_norm_n_13_crt_51_depth_23_base_20.json"
    (
        crt_depth,
        stddev_e_encoding,
        stddev_e_hardcode,
        stddev_e_p,
        p,
        estimated_secpar,
        size,
    ) = find_params(
        secpar,
        n,
        d,
        base_bits,
        crt_bits,
        max_crt_depth,
        input_size,
        input_width,
        norms_path,
    )

    print(f"crt_depth: {crt_depth}")
    print(f"q: {2**(crt_bits * crt_depth)}, log_2 q: {crt_bits * crt_depth}")
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
        base_bits,
        input_size,
        input_width,
        norms_path,
        crt_bits,
        crt_depth,
        stddev_e_encoding,
        stddev_e_hardcode,
        stddev_e_p,
        p,
        estimated_secpar,
        size,
    )


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

# base = 2**20
# stddev_e_p_d: 3.190600825411292529310003374121151864528656005859375
# sqrt_secpar_d: 9
# base_d: 1048576
# m_d: 54
# bound_c_d: 1864799774993447876594399.948890683020389953429795654643834961427216467821921241920790635049343109131
# base-dependent error: 864998612168234998077837623509275354.7074875601141009832394687499579344347466758335940539836883544923
# final_err: 8746810258343107065957007438076535207.556248084800279398452795173157270930067594587421772665651598100
# log_final_err: 123

# base = 2**10
# stddev_e_p_d: 3.190600825411292529310003374121151864528656005859375
# sqrt_secpar_d: 9
# base_d: 1024
# m_d: 105
# bound_c_d: 3226032414224086806.362468960186188167017787808787734525376673639954715720745692664195303223095834255
# base-dependent error: 2838726834371627289030038350.712007900858086829923815457123003588900747584666817147081019356846809387
# final_err: 27825131015929632855862032467540.17738662940100282680672935342417228936988880152276760954647193919046
# log_final_err: 105

# base = 2**5
# stddev_e_p_d: 3.190600825411292529310003374121151864528656005859375
# sqrt_secpar_d: 9
# base_d: 32
# m_d: 210
# bound_c_d: 6329877527671136.947140365292906179514225175286667173341021758352098275888301509009212231227081701945
# base-dependent error: 337571862160499519700039.9098414631209992674600239374209028230911767328852610484113405675543617689982
# final_err: 109188471321800357965049573411.3224853691748555365478958752828378103720729733878934854573778100501152
# log_final_err: 97
