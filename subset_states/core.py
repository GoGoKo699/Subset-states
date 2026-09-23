from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import log2, sqrt
from numbers import Integral, Real
from typing import Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import digamma, gammaln, logsumexp, xlogy


@dataclass(frozen=True)
class SummaryStats:
    """Summary statistics for repeated entropy samples."""

    count: int
    mean: float
    std: float
    sem: float
    minimum: float
    maximum: float


def _validate_positive_n(n: int) -> None:
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral) or n <= 0:
        raise ValueError(f"n must be a positive integer; received n={n}.")


def _validate_even_n(n: int) -> None:
    _validate_positive_n(n)
    if n % 2 != 0:
        raise ValueError(f"n must be a positive even integer; received n={n}.")


def _validate_support_size(n: int, m: int, *, allow_zero: bool = False) -> None:
    N = 1 << n
    lower = 0 if allow_zero else 1
    if isinstance(m, (bool, np.bool_)) or not isinstance(m, Integral) or not (lower <= m <= N):
        raise ValueError(f"m must satisfy {lower} <= m <= 2**n={N}; received m={m}.")


def _validated_support(n: int, support: ArrayLike) -> NDArray[np.int64]:
    """Validate before conversion so fractional or unsigned labels cannot wrap."""

    _validate_positive_n(n)
    labels = np.asarray(support)
    if labels.ndim != 1 or labels.size == 0:
        raise ValueError("support must be a nonempty one-dimensional array of integer basis labels.")
    if labels.dtype.kind not in "iu":
        raise ValueError("support labels must be integers.")
    N = 1 << n
    if np.any(labels < 0) or np.any(labels >= N):
        raise ValueError(f"support labels must lie in [0, {N - 1}].")
    if np.unique(labels).size != labels.size:
        raise ValueError("support contains duplicate basis labels.")
    if np.any(labels > np.iinfo(np.int64).max):
        raise ValueError("support labels exceed the supported int64 range.")
    return labels.astype(np.int64, copy=False)


def _validate_atol(atol: float) -> None:
    if not isinstance(atol, Real) or not np.isfinite(atol) or atol < 0:
        raise ValueError("atol must be a finite nonnegative real number.")


def natural_right_bits(n: int) -> tuple[int, ...]:
    """Return the natural right subsystem: the low-order n/2 bits."""

    _validate_even_n(n)
    return tuple(range(n // 2))


def _left_bits(n: int, right_bits: Sequence[int]) -> tuple[int, ...]:
    if any(isinstance(bit, (bool, np.bool_)) or not isinstance(bit, Integral) for bit in right_bits):
        raise ValueError("right_bits must contain integer bit positions.")
    right = set(right_bits)
    if len(right) != len(right_bits):
        raise ValueError("right_bits contains repeated bit positions.")
    if any(bit < 0 or bit >= n for bit in right):
        raise ValueError(f"right_bits must be between 0 and {n - 1}.")
    return tuple(bit for bit in range(n) if bit not in right)


def random_subset(n: int, m: int, rng: np.random.Generator) -> NDArray[np.int64]:
    """Sample a uniformly random subset of {0, ..., 2**n - 1} of size m."""

    _validate_positive_n(n)
    _validate_support_size(n, m)
    return np.sort(rng.choice(1 << n, size=m, replace=False).astype(np.int64))


def project_bits(values: ArrayLike, bit_positions: Sequence[int]) -> NDArray[np.int64]:
    """Project integer labels onto selected bit positions."""

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
    """Return the bipartite coefficient matrix for an equal-amplitude support.

    The first position in each subsystem's bit list is its least significant
    bit. With ``normalize=True``, each occupied entry equals 1/sqrt(M). For
    coefficients C[a, b], the right reduced state is ``C.T @ C.conj()``
    (here C is real). Set
    ``normalize=False`` to obtain the binary support-incidence matrix.
    """

    _validate_even_n(n)
    support_arr = _validated_support(n, support)
    if normalize and np.dtype(dtype).kind not in "fc":
        raise ValueError("normalized coefficients require a floating-point or complex dtype.")

    right = tuple(right_bits) if right_bits is not None else natural_right_bits(n)
    if len(right) != n // 2:
        raise ValueError(f"balanced bipartition requires {n // 2} right bits.")
    left = _left_bits(n, right)

    matrix = np.zeros((1 << len(left), 1 << len(right)), dtype=dtype)
    rows = project_bits(support_arr, left)
    cols = project_bits(support_arr, right)
    amp = 1.0 / sqrt(support_arr.size) if normalize else 1.0
    matrix[rows, cols] = amp
    return matrix


def matrix_from_state_vector(
    n: int,
    state: ArrayLike,
    right_bits: Sequence[int] | None = None,
) -> NDArray:
    """Return C[a, b] in the same basis order as :func:`matrix_from_support`.

    The first position in each subsystem's bit list is its least significant
    bit. The right reduced state is ``C.T @ C.conj()``; ``C.conj().T @ C``
    is its transpose and has the same spectrum, but generally different entries.
    The vector need not be normalized.
    """

    _validate_even_n(n)
    vector = np.asarray(state)
    N = 1 << n
    if vector.shape != (N,):
        raise ValueError(f"state must have shape ({N},); received {vector.shape}.")
    if vector.dtype.kind not in "biufc" or not np.all(np.isfinite(vector)):
        raise ValueError("state must contain finite numeric amplitudes.")
    right = tuple(right_bits) if right_bits is not None else natural_right_bits(n)
    if len(right) != n // 2:
        raise ValueError(f"balanced bipartition requires {n // 2} right bits.")
    left = _left_bits(n, right)

    # Reshape traverses axes from most to least significant within each index.
    left_axes = [n - 1 - bit for bit in reversed(left)]
    right_axes = [n - 1 - bit for bit in reversed(right)]
    tensor = vector.reshape((2,) * n)
    return np.transpose(tensor, left_axes + right_axes).reshape(1 << len(left), 1 << len(right))


def spectrum_from_matrix(matrix: ArrayLike, *, atol: float = 1e-14) -> NDArray[np.float64]:
    """Return normalized non-zero reduced-state eigenvalues using a stable SVD."""

    _validate_atol(atol)
    mat = np.asarray(matrix)
    if mat.ndim != 2 or mat.size == 0:
        raise ValueError("matrix must be a nonempty two-dimensional array.")
    if mat.dtype.kind not in "biufc" or not np.all(np.isfinite(mat)):
        raise ValueError("matrix must contain finite numeric coefficients.")
    singular_values = np.linalg.svd(mat, compute_uv=False)
    eigvals = np.real_if_close(singular_values * singular_values).astype(np.float64)
    eigvals[eigvals < atol] = 0.0
    eigvals = eigvals[eigvals > atol]
    total = eigvals.sum()
    if total <= 0:
        raise ValueError("density-matrix spectrum has zero trace.")
    return np.sort(eigvals / total)


def renyi_entropy_from_spectrum(
    eigvals: ArrayLike,
    order: float = 1.0,
    *,
    atol: float = 1e-14,
) -> float:
    """Return the Rényi entropy in bits; order=1 gives von Neumann entropy.

    A finite, nonnegative spectrum is normalized after entries at or below
    ``atol`` are removed. Negative roundoff no smaller than ``-atol`` is allowed.
    The order must be positive; positive infinity gives the min-entropy.
    """

    _validate_atol(atol)
    if not isinstance(order, Real) or np.isnan(order) or order <= 0:
        raise ValueError("Rényi order must be a positive real number or positive infinity.")
    p = np.asarray(eigvals)
    if p.ndim != 1 or p.size == 0 or p.dtype.kind not in "biuf":
        raise ValueError("spectrum must be a nonempty one-dimensional array of real eigenvalues.")
    p = p.astype(np.float64)
    if not np.all(np.isfinite(p)) or np.any(p < -atol):
        raise ValueError("spectrum must be finite and nonnegative within atol.")
    p = p[p > atol]
    if p.size == 0:
        raise ValueError("density-matrix spectrum has zero trace at the requested atol.")
    # Scaling first avoids overflow when normalizing a finite input spectrum.
    p = p / p.max()
    p = p / p.sum()
    if order == 1.0:
        return float(-np.sum(xlogy(p, p)) / np.log(2.0))
    if np.isinf(order):
        return float(-np.log2(np.max(p)))
    log_p = np.log(p)
    if abs(order - 1.0) < 0.1:
        # log(sum(p**order)) near one, without subtracting nearly equal numbers.
        delta = np.sum(p * np.expm1((order - 1.0) * log_p))
        return float(np.log1p(delta) / ((1.0 - order) * np.log(2.0)))
    log_max = np.max(log_p)
    with np.errstate(over="ignore"):
        tail = logsumexp(order * (log_p - log_max))
    return float((order / (1.0 - order) * log_max + tail / (1.0 - order)) / np.log(2.0))


def entropy_from_support(
    n: int,
    support: ArrayLike,
    right_bits: Sequence[int] | None = None,
    *,
    order: float = 1.0,
) -> float:
    matrix = matrix_from_support(n, support, right_bits)
    return renyi_entropy_from_spectrum(spectrum_from_matrix(matrix), order=order)


def entropy_from_state_vector(
    n: int,
    state: ArrayLike,
    right_bits: Sequence[int] | None = None,
    *,
    order: float = 1.0,
) -> float:
    matrix = matrix_from_state_vector(n, state, right_bits)
    return renyi_entropy_from_spectrum(spectrum_from_matrix(matrix), order=order)


def qft_state_from_support(n: int, support: ArrayLike) -> NDArray[np.complex128]:
    """Return the positive-exponent QFT of an equal-amplitude subset state.

    NumPy's ``ifft(..., norm='ortho')`` implements the convention
    exp(+2πijk/N)/sqrt(N), matching the research notes. For the real
    input states used here, switching the sign only complex-conjugates the output
    and therefore leaves all reported entropies unchanged.
    """

    support_arr = _validated_support(n, support)
    N = 1 << n
    vector = np.zeros(N, dtype=np.complex128)
    vector[support_arr] = 1.0 / sqrt(support_arr.size)
    return np.fft.ifft(vector, norm="ortho")


def random_balanced_partition(n: int, rng: np.random.Generator) -> tuple[int, ...]:
    """Sample a unique balanced bipartition, represented by right-subsystem bits."""

    _validate_even_n(n)
    other = rng.choice(np.arange(1, n), size=n // 2 - 1, replace=False)
    return tuple(sorted([0, *other.tolist()]))


def balanced_bipartitions(n: int) -> Iterable[tuple[int, ...]]:
    """Yield every balanced bipartition once, identifying complements."""

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


def fixed_cardinality_mean_spectrum(n: int, m: int) -> tuple[float, float]:
    """Eigenvalues of E[rho] for the exact fixed-cardinality ensemble.

    For a balanced cut with d=2^(n/2) and N=d^2, the diagonal entries of E[rho]
    equal 1/d and the off-diagonal entries equal (M-1)/(d(N-1)).  The returned
    pair is ``(uniform_mode, orthogonal_mode)``; the second eigenvalue has
    multiplicity d-1.
    """

    _validate_even_n(n)
    _validate_support_size(n, m)
    d = 1 << (n // 2)
    N = 1 << n
    off_diagonal = (m - 1) / (d * (N - 1))
    uniform_mode = 1.0 / d + (d - 1) * off_diagonal
    orthogonal_mode = 1.0 / d - off_diagonal
    return float(uniform_mode), float(orthogonal_mode)


def mean_matrix_uniform_eigenvalue(n: int, m: int) -> float:
    """Uniform-mode eigenvalue of the exact mean reduced density matrix."""

    return fixed_cardinality_mean_spectrum(n, m)[0]


def dense_bulk_approximation(n: int, m: int) -> float:
    """Dense-bulk entropy ansatz using the exact mean-matrix uniform mode.

    The first moment is exact for the fixed-cardinality ensemble; treating its
    uniform-mode eigenvalue as the separated eigenvalue of a typical sample is
    still a mean-field approximation.
    """

    _validate_even_n(n)
    _validate_support_size(n, m, allow_zero=True)
    if m == 0:
        # Retain the old formal Page-limit value for backward compatibility.
        return n / 2.0 - log2(np.e) / 2.0
    if m == (1 << n):
        return 0.0

    d = 1 << (n // 2)
    lam = mean_matrix_uniform_eigenvalue(n, m)
    remainder = 1.0 - lam
    entropy_nats = (
        remainder * (np.log(d) - np.log(remainder) - 0.5)
        - lam * np.log(lam)
    )
    return float(entropy_nats / np.log(2.0))


def hypergeometric_occupancy_approximation(n: int, m: int) -> float:
    """Exact fixed-cardinality diagonal-occupancy entropy in bits.

    A column contains d of the N=d^2 support locations.  Its occupancy W is
    Hypergeometric(N, d, M), not Binomial(M, 1/d).  The returned quantity is the
    expected Shannon entropy of the diagonal of rho.  It is an upper bound on
    the von Neumann entropy and a sparse-regime approximation when off-diagonal
    coherences are negligible.
    """

    _validate_even_n(n)
    _validate_support_size(n, m, allow_zero=True)
    if m == 0:
        return 0.0
    d = 1 << (n // 2)
    N = 1 << n
    w_min = max(0, m - (N - d))
    w_max = min(d, m)
    w = np.arange(w_min, w_max + 1, dtype=np.float64)

    def log_choose(total: int, chosen: np.ndarray | float) -> np.ndarray:
        return gammaln(total + 1.0) - gammaln(chosen + 1.0) - gammaln(total - chosen + 1.0)

    log_prob = (
        log_choose(d, w)
        + log_choose(N - d, m - w)
        - (gammaln(N + 1.0) - gammaln(m + 1.0) - gammaln(N - m + 1.0))
    )
    log_prob -= logsumexp(log_prob)
    prob = np.exp(log_prob)
    positive = w > 0
    p_diag = w[positive] / m
    return float(-d * np.sum(prob[positive] * xlogy(p_diag, p_diag)) / np.log(2.0))



def fixed_cardinality_average_purity(n: int, m: int) -> float:
    """Exact ensemble-average purity for the balanced fixed-cardinality ensemble.

    For N=2**n=d**2 and a support chosen uniformly among all M-element subsets,
    this returns E[Tr(rho_R**2)] as derived in the research notes.
    """

    _validate_even_n(n)
    _validate_support_size(n, m)
    d = 1 << (n // 2)
    N = 1 << n
    return float(
        1.0 / m
        + 2.0 * (d - 1) * (m - 1) / (m * (N - 1))
        + (d - 1) ** 2 * (m - 1) * (m - 2) * (m - 3)
        / (m * (N - 1) * (N - 2) * (N - 3))
    )


def average_entropy_lower_bound_from_purity(n: int, m: int) -> float:
    """Rigorous lower bound -log2(E[Tr(rho_R**2)]) on average entropy."""

    return float(-np.log2(fixed_cardinality_average_purity(n, m)))


def residue_class_entropy_ceiling(n: int, support: ArrayLike, matched_low_bits: int) -> float:
    """Entropy ceiling from populations modulo 2**matched_low_bits.

    For the natural balanced cut, if p_r is the support population in residue
    class r modulo 2**t, then S <= n/2 - t + H(p).  The case t=1 is the
    parity ceiling.
    """

    _validate_even_n(n)
    if isinstance(matched_low_bits, (bool, np.bool_)) or not isinstance(matched_low_bits, Integral):
        raise ValueError("matched_low_bits must be an integer.")
    t = int(matched_low_bits)
    if not (0 <= t <= n // 2):
        raise ValueError(f"matched_low_bits must lie in [0, {n // 2}].")
    support_arr = _validated_support(n, support)
    modulus = 1 << t
    counts = np.bincount(support_arr % modulus, minlength=modulus).astype(float)
    probabilities = counts[counts > 0] / support_arr.size
    shannon = float(-np.sum(probabilities * np.log2(probabilities))) if probabilities.size else 0.0
    return float(n / 2 - t + shannon)

def tnm_approximation(n: int, m: int) -> float:
    """Backward-compatible alias for :func:`dense_bulk_approximation`."""

    return dense_bulk_approximation(n, m)


def dnm_approximation(n: int, m: int) -> float:
    """Backward-compatible alias for the exact hypergeometric occupancy curve."""

    return hypergeometric_occupancy_approximation(n, m)


def exact_page_entropy_bits(n: int) -> float:
    """Exact Page average entropy in bits for a balanced n-qubit bipartition."""

    _validate_even_n(n)
    d = 1 << (n // 2)
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
            if power > (N - 1) // p:
                break
            power *= p
    return omega
