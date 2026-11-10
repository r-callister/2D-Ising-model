"""
pytest configuration: use a non-interactive matplotlib backend so that tests never open
plot windows.
"""
import matplotlib

matplotlib.use("Agg")
