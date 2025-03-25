import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(6, 6))

a = np.array([3 * np.sqrt(3) / 2, 3 / 2])
b = np.array([2 / 2, 2 * np.sqrt(3) / 2])
k = np.dot(a, b) / np.dot(a, a)
p = k * a
e = np.array([b[0] - p[0], b[1] - p[1]])
r = b - 2 * p

plt.quiver(0, 0, a[0], a[1], angles="xy", scale_units="xy", scale=1, color="blue", label="a")
plt.quiver(0, 0, b[0], b[1], angles="xy", scale_units="xy", scale=1, color="red", label="b")
plt.quiver(0, 0, p[0], p[1], angles="xy", scale_units="xy", scale=1, color="green", label="projection")
plt.quiver(0, 0, e[0], e[1], angles="xy", scale_units="xy", scale=1, color="green", label="error")
plt.quiver(0, 0, r[0], r[1], angles="xy", scale_units="xy", scale=1, color="red", label="reflection")

plt.text(a[0], a[1], "a", color="blue", fontsize=12, ha="right")
plt.text(b[0], b[1], "b", color="red", fontsize=12, ha="right")
plt.text(p[0], p[1], "projection", color="green", fontsize=12, ha="left")
plt.text(e[0], e[1], "error", color="green", fontsize=12, ha="right")
plt.text(r[0], r[1], "reflection", color="red", fontsize=12, ha="right")

plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.axhline(0, color="black", linewidth=0.5)
plt.axvline(0, color="black", linewidth=0.5)
plt.show()
