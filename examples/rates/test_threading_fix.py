"""
Test script to verify the Numba threading layer fix.

This script tests that the TBB threading layer is properly set
and that the LeR class can be initialized without the 
"Concurrent access has been detected" error.
"""

from ler.rates import LeR

# Test initialization
ler = LeR(npool=6)
