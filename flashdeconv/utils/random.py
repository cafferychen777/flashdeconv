"""Random state utilities for FlashDeconv.

This module provides utilities for handling random state in a manner consistent
with scikit-learn and scanpy conventions. Using local RandomState instances
instead of global np.random.seed() ensures thread safety and avoids polluting
global state.
"""

import numpy as np
from typing import Union, Optional

# Type alias for random state parameter
RandomStateLike = Union[None, int, np.random.RandomState]


def check_random_state(seed: RandomStateLike) -> np.random.RandomState:
    """
    Turn seed into a np.random.RandomState instance.

    This function follows the scikit-learn convention for handling random state,
    ensuring consistent and reproducible random number generation without
    polluting global state.

    Parameters
    ----------
    seed : None, int, or np.random.RandomState
        If seed is None, return the global RandomState singleton used by np.random.
        If seed is an int, return a new RandomState instance seeded with seed.
        If seed is already a RandomState instance, return it unchanged.

    Returns
    -------
    np.random.RandomState
        A RandomState object for generating random numbers.

    Raises
    ------
    ValueError
        If seed is not a valid type.

    Examples
    --------
    >>> rng = check_random_state(42)
    >>> rng.randint(0, 10)
    1
    >>> rng2 = check_random_state(42)
    >>> rng2.randint(0, 10)  # Same result
    1

    Notes
    -----
    This approach is preferred over using np.random.seed() directly because:

    1. Thread safety: Local RandomState instances are independent
    2. Reproducibility: Each call with the same seed produces identical sequences
    3. No side effects: Does not affect other code using np.random

    See Also
    --------
    sklearn.utils.check_random_state : The original implementation this is based on.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (int, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError(
        f"'{seed}' cannot be used to seed a numpy.random.RandomState instance. "
        f"Expected None, int, or np.random.RandomState, got {type(seed)}."
    )
