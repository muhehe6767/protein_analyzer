import numpy as np

def predict_secondary_structure(sequence: str) -> dict:
    """Mock/Baseline secondary structure prediction (H=Helix, E=Sheet, C=Coil)."""
    rng = np.random.default_rng(len(sequence))
    pred = "".join(rng.choice(["H", "E", "C"], size=len(sequence)))

    seq_len = len(sequence)
    if seq_len == 0:
        return {
            "predicted_sequence": "",
            "helix_pct": 0.0,
            "sheet_pct": 0.0,
            "coil_pct": 0.0,
        }

    return {
        "predicted_sequence": pred,
        "helix_pct": round((pred.count("H") / seq_len) * 100, 2),
        "sheet_pct": round((pred.count("E") / seq_len) * 100, 2),
        "coil_pct": round((pred.count("C") / seq_len) * 100, 2)
    }
