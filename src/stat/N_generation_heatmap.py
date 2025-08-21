"""negative hypergeometric-like"""


def expected_batches(S: int, B: int) -> float:
    """
    
    """
    max_batches = S // B
    expectation = 0.0
    for t in range(1, max_batches + 1):
        prob = 1.0
        for i in range(1, t):
            prob *= (1 - B / (S - (i - 1) * B))
        prob *= (B / (S - (t - 1) * B))
        expectation += t * prob
    return expectation

