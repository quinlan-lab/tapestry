"""Shared decoding of per-sample VCF PS values (including cyvcf2 sentinels)."""


def read_ps(ps_all, sample_index):
    """Return a positive PS as a string, or None for missing/invalid values."""
    if ps_all is None:
        return None
    val = ps_all[sample_index]
    try:
        value = int(val[0]) if hasattr(val, '__len__') else int(val)
    except (TypeError, ValueError, IndexError):
        return None
    return str(value) if value > 0 else None
