import numpy as np
from sklearn.datasets import fetch_olivetti_faces
import matplotlib.pyplot as plt


def preview(images):
    """
    Preview Olivetti faces
    """
    _, axes = plt.subplots(2, 5, figsize=(10, 4))
    for i, ax in enumerate(axes.flat):
        ax.imshow(images[i], cmap="gray")
        ax.axis("off")
    plt.tight_layout()
    plt.show()


def read_archive():
    """
    Read random 10 Olivetti faces from archive
    """
    faces = fetch_olivetti_faces(shuffle=True, random_state=42)
    images = faces.images
    indexes = np.random.choice(images.shape[0], size=10, replace=False)
    return images[indexes]


def pipe_archive(images):
    """
    Write preprocessed Olivetti faces to disk
    """
    flattened = images.astype(np.float64).reshape(len(images), -1)
    binary = flattened.tobytes()
    filename = "input.bin"
    with open(filename, "wb") as file:
        file.write(binary)


def read_binary(filename):
    """
    Read binary Olivetti faces from disk
    """
    with open(filename, "rb") as file:
        binary_data = file.read()
    images = np.frombuffer(binary_data, dtype=np.float64)
    images = images.reshape(10, *(64, 64))
    return images


# images = read_archive()
# preview(images)
# pipe_archive(images)

# outputEvdPm = read_binary("outputEvdPm.bin")
# preview(outputEvdPm)
# outputSvdPm = read_binary("outputSvdPm.bin")
# preview(outputSvdPm)
# outputEvdQrm = read_binary("outputEvdQrm.bin")
# preview(outputEvdQrm)
# outputSvdQrm = read_binary("outputSvdQrm.bin")
# preview(outputSvdQrm)
