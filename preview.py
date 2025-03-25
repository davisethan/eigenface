import numpy as np
from sklearn.datasets import fetch_olivetti_faces
import matplotlib.pyplot as plt


def preview(images):
    """
    Preview Olivetti faces
    """
    _, axes = plt.subplots(2, 5, figsize=(10, 5))
    for i, ax in enumerate(axes.flat):
        ax.imshow(images[i], cmap="gray")
        ax.axis("off")
    plt.tight_layout()
    plt.show()


def read_input():
    """
    Read random 10 Olivetti faces from archive
    """
    faces = fetch_olivetti_faces(shuffle=True, random_state=42)
    images = faces.images
    indexes = np.random.choice(images.shape[0], size=10, replace=False)
    return images[indexes]


def write_input(images):
    """
    Write preprocessed Olivetti faces to disk
    """
    flattened = images.astype(np.float64).reshape(len(images), -1)
    binary = flattened.tobytes()
    filename = "input.bin"
    with open(filename, "wb") as file:
        file.write(binary)


def read_output(filename):
    """
    Read processed Olivetti faces from disk
    """
    with open(filename, "rb") as file:
        binary_data = file.read()
    images = np.frombuffer(binary_data, dtype=np.float64)
    images = images.reshape(10, *(64, 64))
    return images


# input = read_input()
# preview(input)
# write_input(input)

# outputEvdPm = read_output("outputEvdPm.bin")
# preview(outputEvdPm)
# outputSvdPm = read_output("outputSvdPm.bin")
# preview(outputSvdPm)
# outputEvdQrm = read_output("outputEvdQrm.bin")
# preview(outputEvdQrm)
# outputSvdQrm = read_output("outputSvdQrm.bin")
# preview(outputSvdQrm)
