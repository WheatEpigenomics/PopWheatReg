# Resolve a sibling frozen module by name (used for lazy, in-function imports).
import importlib


def load(name):
    return importlib.import_module("." + name, __package__)
