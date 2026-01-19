import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

try:
    data = np.array(0.5)  # 0-d array
    kde = gaussian_kde(data)
except Exception as e:
    print(f"0-d array error: {e}")

try:
    data = np.array([0.5])  # 1-d array, 1 element
    kde = gaussian_kde(data)
except Exception as e:
    print(f"1-d array 1 element error: {e}")

try:
    data = np.array([0.5, 0.6])  # 1-d array, 2 elements
    kde = gaussian_kde(data)
    print("1-d array 2 elements: success")
except Exception as e:
    print(f"1-d array 2 elements error: {e}")
