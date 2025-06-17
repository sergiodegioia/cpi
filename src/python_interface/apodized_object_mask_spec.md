# Apodized Object Input to the CPI Simulator

This document specifies the format and structure required for defining **apodized transmissive planar objects** as input to the CPI Simulator. These masks are stored in TIFF files and represent complex-valued transmission functions.

## Format Summary

- **File type**: TIFF
- **Bit depth**: 16-bit unsigned integer (`uint16`)
- **Channels**: 2 (real and imaginary parts)
- **Shape**: Square (e.g., 512×512 or 1024×1024)
- **Compression**: None (uncompressed TIFF)

## Representation

Each pixel in the image represents a complex number:

- Channel 0: Real part
- Channel 1: Imaginary part

The complex number \( z = A \cdot e^{i\phi} \) is mapped to two unsigned 16-bit integers:

- **Real part**: \( 	ext{Re}(z) \in [-1, 1] \) → mapped to \( [0, 65535] \)
- **Imaginary part**: \( 	ext{Im}(z) \in [-1, 1] \) → mapped to \( [0, 65535] \)

This mapping is achieved using:

```python
mapped = ((value + 1.0) * 0.5 * 65535).astype(np.uint16)
```

## Generation Example (Python)

```python
import numpy as np
import tifffile
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = False

def generate_apodized_object_mask(shape=(512, 512), type='binary', directory="."):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    filename = f"mask_{type}.tiff"
    save_path = directory / filename

    x = np.linspace(-1, 1, shape[1])
    y = np.linspace(-1, 1, shape[0])
    X, Y = np.meshgrid(x, y)
    
    amplitude = np.exp(-((X**2 + Y**2) / 0.1))           # range [0, 1]
    phase = np.pi * np.sin(2 * np.pi * X)               # range [-π, π]
    complex_mask = amplitude * np.exp(1j * phase)       # complex64/128

    # Map real and imaginary from [-1, 1] → [0, 65535]
    re_scaled = ((np.real(complex_mask) + 1.0) * 0.5 * 65535).astype(np.uint16)
    im_scaled = ((np.imag(complex_mask) + 1.0) * 0.5 * 65535).astype(np.uint16)

    mask_2ch = np.stack((re_scaled, im_scaled), axis=-1)  # shape: (H, W, 2), dtype=uint16

    tifffile.imwrite(str(save_path), mask_2ch, compression=None)
    print(f"[INFO] 16-bit 2-channel TIFF mask saved to: {Path(save_path).resolve()}")
    amplitude = np.abs(complex_mask)

    plt.imshow(amplitude, cmap='gray')
    plt.title("Amplitude of Apodized Mask")
    plt.colorbar()
    plt.show()

    print(f"[INFO] Mask saved to: {save_path.resolve()}")

if __name__ == "__main__":
    generate_apodized_object_mask()
```

