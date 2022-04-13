# Imported Packages
import numpy as np
import matplotlib.pyplot as plt

# User Defined Function
a = np.linspace(-np.sqrt(10), np.sqrt(10), 11)
b = np.sign(a)*a**2
x, y = np.meshgrid(b, b)
plt.plot(x.flatten(), y.flatten(), '*')
plt.show()
