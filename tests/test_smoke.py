"""Existing native input-deck smoke test."""
from support import NativeTest


class Smoke(NativeTest):
    def test_existing_input_deck(self):
        self.smoke()
