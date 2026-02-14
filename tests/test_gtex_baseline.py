"""Tests for src/6_gtex_baseline.py — classification logic."""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from importlib import import_module
gtex_mod = import_module('6_gtex_baseline')


# ===================================================================
# classify_silencing
# ===================================================================

class TestClassifySilencing:
    """Test expression classification against GTEx baselines."""

    # --- No GTEx data ---
    def test_none_gtex_returns_no_data(self):
        assert gtex_mod.classify_silencing(10.0, None) == "no GTEx data"

    # --- Tissue-normal silence: both < 1 TPM ---
    def test_both_low(self):
        assert gtex_mod.classify_silencing(0.1, 0.2) == "tissue-normal silence"

    def test_both_zero(self):
        assert gtex_mod.classify_silencing(0.0, 0.0) == "tissue-normal silence"

    # --- Tumour-specific silencing: tumour low, GTEx high ---
    def test_tumour_silenced_gtex_expressed(self):
        assert gtex_mod.classify_silencing(0.3, 5.0) == "tumour-specific silencing"

    def test_tumour_zero_gtex_expressed(self):
        assert gtex_mod.classify_silencing(0.0, 10.0) == "tumour-specific silencing"

    # --- Tumour over-expression: tumour ≥ 4× GTEx ---
    def test_overexpression(self):
        # 20 / 4 = 5.0 → ratio = 5.0 ≥ 4.0
        assert gtex_mod.classify_silencing(20.0, 4.0) == "tumour over-expression"

    def test_overexpression_exact_threshold(self):
        # 4.0 / 1.0 = 4.0 → exactly at threshold → over-expression
        assert gtex_mod.classify_silencing(4.0, 1.0) == "tumour over-expression"

    def test_overexpression_large_ratio(self):
        assert gtex_mod.classify_silencing(100.0, 2.0) == "tumour over-expression"

    # --- Comparable: both expressed, tumour < 4× GTEx ---
    def test_comparable(self):
        assert gtex_mod.classify_silencing(5.0, 4.0) == "comparable"

    def test_comparable_equal(self):
        assert gtex_mod.classify_silencing(10.0, 10.0) == "comparable"

    def test_comparable_tumour_slightly_lower(self):
        assert gtex_mod.classify_silencing(3.0, 10.0) == "comparable"

    # --- Edge: threshold boundaries ---
    def test_tumour_at_threshold_gtex_below(self):
        """Tumour = 1.0 (expressed), GTEx = 0.5 (not expressed)."""
        # tumour_expressed=True, gtex_expressed=False
        # Falls to the over-expression / comparable branch
        # tumour / gtex = 1.0 / 0.5 = 2.0 < 4.0 → comparable
        assert gtex_mod.classify_silencing(1.0, 0.5) == "comparable"

    def test_tumour_below_threshold_gtex_at_threshold(self):
        """Tumour = 0.9 (not expressed), GTEx = 1.0 (expressed)."""
        assert gtex_mod.classify_silencing(0.9, 1.0) == "tumour-specific silencing"

    def test_gtex_zero_tumour_expressed(self):
        """GTEx = 0, tumour expressed — gtex_expressed is False."""
        # Both expressed check: tumour yes, gtex no
        # Falls to: tumour_expressed=True, gtex_expressed=False
        # Then: gtex_tpm > 0 is False → skip over-expression → comparable
        assert gtex_mod.classify_silencing(10.0, 0.0) == "comparable"
