# src/sequence.py

def parse_fasta(fasta_str: str):
    """Parses FASTA string and returns (header, sequence)."""
    lines = fasta_str.strip().split("\n")
    if lines and lines[0].startswith(">"):
        header = lines[0][1:].strip()
        sequence = "".join([line.strip() for line in lines[1:]]).upper()
    else:
        header = "User Input Sequence"
        sequence = "".join([line.strip() for line in lines]).upper()
    return header, sequence


def validate_sequence(seq: str):
    """Validates if sequence contains valid amino acid characters."""
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    if not seq:
        return False, "Sequence is empty."
    invalid_chars = set(seq) - valid_aa
    if invalid_chars:
        return False, f"Invalid characters found: {', '.join(invalid_chars)}"
    return True, "Valid sequence"


def calculate_sequence_properties(seq: str):
    """Calculates basic properties for the given protein sequence."""
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    analyzed_seq = ProteinAnalysis(seq)
    return {
        "molecular_weight": analyzed_seq.molecular_weight(),
        "isoelectric_point": analyzed_seq.isoelectric_point(),
        "amino_acid_counts": analyzed_seq.count_amino_acids(),
    }