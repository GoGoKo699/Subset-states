from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import log2, sqrt
from typing import Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import gammaln, logsumexp, xlogy, digamma


@dataclass(frozen=True)
class SummaryStats:
    """Summary statistics for repeated entropy samples."""

    count: int
    mean: float
    std: float
    sem: float
    minimum: float
    maximum: float


def _validate_even_n(n: int) -> None:
    if n <= 0 or n % 2 != 0:
        raise ValueError(f"n must be a positive even integer; received n={n}.")


def natural_right_bits(n: int) -> tuple[int, ...]:
    """Return the natural right subsystem: the low-order n/2 bits."""

    _validate_even_n(n)
    return tuple(range(n // 2))


def _left_bits(n: int, right_bits: Sequence[int]) -> tuple[int, ...]:
    right = set(right_bits)
    if len(right) != len(right_bits):
        raise ValueError("right_bits contains repeated bit positions.")
    if any(bit < 0 or bit >= n for bit in right):
        raise ValueError(f"right_bits must be between 0 and {n - 1}.")
    return tuple(bit for bit in range(n) if bit not in right)


def random_subset(n: int, m: int, rng: np.random.Generator) -> NDArray[np.int64]:
    """Sample a uniformly random subset of {0, ..., 2**n - 1} of size m."""

    if n <= 0:
        raise ValueError("n must be positive.")
    N = 1 << n
    if not (1 <= m <= N):
        raise ValueError(f"m must satisfy 1 <= m <= 2**n={N}; received m={m}.")
    return np.sort(rng.choice(N, size=m, replace=False).astype(np.int64))


def project_bits(values: ArrayLike, bit_positions: Sequence[int]) -> NDArray[np.int64]:
    """Project integer labels onto selected bit positions.

    Bit positions are counted from the least-significant bit. The first entry of
    ``bit_positions`` becomes the least-significant bit of the projected label.
    """

    values_u = np.asarray(values, dtype=np.uint64)
    out = np.zeros(values_u.shape, dtype=np.uint64)
    for new_bit, old_bit in enumerate(bit_positions):
        out |= ((values_u >> np.uint64(old_bit)) & np.uint64(1)) << np.uint64(new_bit)
    return out.astype(np.int64, copy=False)


def matrix_from_support(
    n: int,
    support: ArrayLike,
    right_bits: Sequence[int] | None = None,
    *,
    normalize: bool = True,
    dtype: np.dtype | type = np.float64,
) -> NDArray:
    """Return the bipartite coefficient matrix Ω for an equal-amplitude support.

    The reduced density matrix of the right subsystem is Ω†Ω. For a subset state
    with M occupied basis states, the non-zero entries of Ω are 1/sqrt(M).
    """

    _validate_even_n(n)
    support_arr = np.asarray(support, dtype=np.int64)
    if support_arr.ndim != 1:
        raise ValueError("support must be a one-dimensional array of basis labels.")
    if support_arr.size == 0:
        raise ValueError("support must contain at least one basis label.")
    N = 1 << n
    if np.any(support_arr < 0) or np.any(support_arr >= N):
        raise ValueError(f"support labels must lie in [0, {N - 1}].")
    if np.unique(support_arr).size != support_arr.size:
        raise ValueError("support contains duplicate basis labels.")

    right = tuple(right_bits) if right_bits is not None else natural_right_bits(n)
    if len(right) != n // 2:
        raise ValueError(f"balanced bipartition requires {n // 2} right bits.")
    left = _left_bits(n, right)

    dim_left = 1 << len(left)
    dim_right = 1 << len(right)
    omega = np.zeros((dim_left, dim_right), dtype=dtype)
    rows = project_bits(support_arr, left)
    cols = project_bits(support_arr, right)
    amp = 1.0 / sqrt(support_arr.size) if normalize else 1.0
    omega[rows, cols] = amp
    return omega


def matrix_from_state_vector(
    n: int,
    state: ArrayLike,
    right_bits: Sequence[int] | None = None,
) -> NDArray:
    """Return the bipartite coefficient matrix Ω for an arbitrary state vector."""

    _validate_even_n(n)
    vector = np.asarray(state)
    N = 1 << n
    if vector.shape != (N,):
        raise ValueError(f"state must have shape ({N},); received {vector.shape}.")
    right = tuple(right_bits) if right_bits is not None else natural_right_bits(n)
    if len(right) != n // 2:
        raise ValueError(f"balanced bipartition requires {n // 2} right bits.")
    left = _left_bits(n, right)

    # NumPy's reshape axes correspond to |x_{n-1}> ⊗ ... ⊗ |x_0>.
    left_axes = [n - 1 - bit for bit in left]
    right_axes = [n - 1 - bit for bit in right]
    tensor = vector.reshape((2,) * n)
    return np.transpose(tensor, left_axes + right_axes).reshape(1 << len(left), 1 << len(right))


def spectrum_from_matrix(omega: ArrayLike, *, atol: float = 1e-14) -> NDArray[np.float64]:
    """Return the non-zero eigenvalues of Ω†Ω using a stable SVD."""

    mat = np.asarray(omega)
    singular_values = np.linalg.svd(mat, compute_uv=False)
    eigvals = np.real_if_close(singular_values * singular_values).astype(np.float64)
    eigvals[eigvals < atol] = 0.0
    eigvals = eigvals[eigvals > atol]
    total = eigvals.sum()
    if total <= 0:
        raise ValueError("density-matrix spectrum has zero trace.")
    # Renormalize tiny numerical trace drift, but leave exact spectra unchanged.
    eigvals = eigvals / total
    return np.sort(eigvals)


def renyi_entropy_from_spectrum(
    eigvals: ArrayLike,
    order: float = 1.0,
    *,
    atol: float = 1e-14,
) -> float:
    """Return the Rényi entropy in bits. order=1 gives von Neumann entropy."""

    p = np.asarray(eigvals, dtype=np.float64)
    p = p[p > atol]
    p = p / p.sum()
    if np.isclose(order, 1.0):
        return float(-np.sum(xlogy(p, p)) / np.log(2.0))
    if np.isinf(order):
        return float(-np.log2(np.max(p)))
    if order <= 0:
        raise ValueError("Rényi order must be positive.")
    return float(np.log2(np.sum(p**order)) / (1.0 - order))


def entropy_from_support(
    n: int,
    support: ArrayLike,
    right_bits: Sequence[int] | None = None,
    *,
    order: float = 1.0,
) -> float:
    omega = matrix_from_support(n, support, right_bits)
    return renyi_entropy_from_spectrum(spectrum_from_matrix(omega), order=order)


def entropy_from_state_vector(
    n: int,
    state: ArrayLike,
    right_bits: Sequence[int] | None = None,
    *,
    order: float = 1.0,
) -> float:
    omega = matrix_from_state_vector(n, state, right_bits)
    return renyi_entropy_from_spectrum(spectrum_from_matrix(omega), order=order)


def qft_state_from_support(n: int, support: ArrayLike) -> NDArray[np.complex128]:
    """Return the normalized QFT of an equal-amplitude subset state."""

    N = 1 << n
    support_arr = np.asarray(support, dtype=np.int64)
    vector = np.zeros(N, dtype=np.complex128)
    vector[support_arr] = 1.0 / sqrt(support_arr.size)
    return np.fft.fft(vector, norm="ortho")


def random_balanced_partition(n: int, rng: np.random.Generator) -> tuple[int, ...]:
    """Sample a unique balanced bipartition, represented by right-subsystem bits.

    The complement of a bipartition gives the same entropy. To avoid double
    counting, bit 0 is always assigned to the right subsystem.
    """

    _validate_even_n(n)
    other = rng.choice(np.arange(1, n), size=n // 2 - 1, replace=False)
    return tuple(sorted([0, *other.tolist()]))


def balanced_bipartitions(n: int) -> Iterable[tuple[int, ...]]:
    """Yield every unique balanced bipartition once.

    The count is binomial(n, n/2)/2. Bit 0 is fixed in the right subsystem to
    identify complementary bipartitions.
    """

    _validate_even_n(n)
    for combo in combinations(range(1, n), n // 2 - 1):
        yield (0, *combo)


def summary_stats(values: ArrayLike) -> SummaryStats:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0:
        raise ValueError("cannot summarize an empty array.")
    std = float(arr.std(ddof=1)) if arr.size > 1 else 0.0
    return SummaryStats(
        count=int(arr.size),
        mean=float(arr.mean()),
        std=std,
        sem=float(std / np.sqrt(arr.size)) if arr.size > 1 else 0.0,
        minimum=float(arr.min()),
        maximum=float(arr.max()),
    )


def tnm_approximation(n: int, m: int) -> float:
    """Analytical approximation T_{N,M} from Eq. (25) of the draft."""

    N = 1 << n
    if m < 0 or m > N:
        raise ValueError(f"m must satisfy 0 <= m <= {N}.")
    c = log2(np.e) / 2.0
    if m == 0:
        return n / 2.0 - c
    if m == N:
        return 0.0
    return float(((N - m) * (1.5 * n - log2(N - m) - c) + m * (n - log2(m))) / N)


def dnm_approximation(n: int, m: int) -> float:
    """Sparse balls-in-boxes approximation D_{N,M} from Appendix C.

    This evaluates the binomial expression directly in log space, avoiding the
    Monte Carlo noise present in the original script.
    """

    _validate_even_n(n)
    if m <= 0:
        return 0.0
    boxes = 1 << (n // 2)
    if m == 1:
        return 0.0
    w = np.arange(1, m + 1, dtype=np.float64)
    log_q = (
        gammaln(m + 1)
        - gammaln(w + 1)
        - gammaln(m - w + 1)
        + w * np.log(1.0 / boxes)
        + (m - w) * np.log1p(-1.0 / boxes)
    )
    weights = np.exp(log_q)
    p = w / m
    return float(-boxes * np.sum(weights * p * np.log2(p)))


def exact_page_entropy_bits(n: int) -> float:
    """Exact Page average entropy in bits for a balanced n-qubit bipartition."""

    _validate_even_n(n)
    d = 1 << (n // 2)
    # Page formula: H_{d^2} - H_d - (d-1)/(2d), converted to bits.
    return float((digamma(d * d + 1) - digamma(d + 1) - (d - 1) / (2 * d)) / np.log(2.0))


def omega_prime_factors_sieve(N: int) -> NDArray[np.uint16]:
    """Return Ω(k): the number of prime factors of k with multiplicity for k<N."""

    if N <= 1:
        return np.zeros(max(N, 0), dtype=np.uint16)
    is_prime = np.ones(N, dtype=bool)
    is_prime[:2] = False
    for p in range(2, int(np.sqrt(N - 1)) + 1):
        if is_prime[p]:
            is_prime[p * p : N : p] = False

    omega = np.zeros(N, dtype=np.uint16)
    primes = np.flatnonzero(is_prime)
    for p in primes:
        power = p
        while power < N:
            omega[power:N:power] += 1
            # Avoid integer overflow and unnecessary work.
            if power > (N - 1) // p:
                break
            power *= p
    return omega
