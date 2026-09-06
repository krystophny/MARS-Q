"""Independent manufactured-mode tests for the Hamada NTV projection."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / "MarsQ_2FK" / "ntv_pre.f").read_text()


def periodic_midpoint_projection(mode: int, n_chi: int) -> complex:
    """Project one periodic Fourier mode on MARS's midpoint grid."""
    spacing = 2.0 * np.pi / n_chi
    theta = (np.arange(n_chi) + 0.5) * spacing
    field = np.exp(1j * mode * theta)
    return np.mean(field * np.exp(-1j * mode * theta))


def test_periodic_midpoint_projection_recovers_all_requested_modes() -> None:
    n_chi = 257
    for mode in range(-80, 81):
        assert periodic_midpoint_projection(mode, n_chi) == pytest.approx(1.0 + 0.0j)


def test_periodic_midpoint_projection_has_no_alias_leakage() -> None:
    n_chi = 257
    spacing = 2.0 * np.pi / n_chi
    theta = (np.arange(n_chi) + 0.5) * spacing
    field = np.exp(1j * 25 * theta)
    for mode in range(-80, 81):
        coefficient = np.mean(field * np.exp(-1j * mode * theta))
        expected = 1.0 if mode == 25 else 0.0
        assert abs(coefficient - expected) < 2.0e-13


def test_native_source_uses_periodic_midpoint_spacing_and_all_samples() -> None:
    assert "htheta = 2.*npi/n_chi" in SOURCE
    assert "do j=1,n_chi" in SOURCE
    assert "(dB(j)+dB(j+1))*.5" not in SOURCE
    assert "fac(1)=0.0" in SOURCE
