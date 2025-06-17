# Object Input to the CPI Simulator

## Title:
**Specification for Planar Object Representation as TIFF Input**

## Purpose:
This document specifies the format, type, and semantics of the object input file required by the **CPI Simulator**. The input describes a transmissive planar object through a 2D grayscale image mask.

---

## File Format

- **File type**: TIFF (Tagged Image File Format)
- **Channel configuration**: Single channel
- **Data type**: Unsigned 8-bit integers (`uint8`)
- **Compression**: **No compression**
- **Color depth**: 8 bits per pixel
- **File extension**: `.tiff`

---

## Image Dimensions

- **Shape**: Square
- **Size**: \( N \times N \), where \( N \) is defined by the simulator configuration
- Example: 512×512, 1264×1264, etc.

---

## Pixel Values and Semantics

- The object is modeled as a **purely amplitude-modulating mask** with two distinct transmission levels.

| Pixel Value | Meaning                              |
|-------------|---------------------------------------|
| `0`         | Fully **opaque** (no light passes)    |
| `255`       | Fully **transparent** (light passes)  |

- **Intermediate values** (1–254) are not allowed in this mode and will be considered invalid unless explicitly enabled by a different configuration.

---

## Interpretation

- The mask defines which regions of the incoming optical field are transmitted or blocked.
- The image is interpreted as a **2D function** over the field of view of the object plane.
- The pixel grid is aligned with the simulator’s spatial sampling, and each pixel represents a uniformly transmissive square patch of the object.

---

## Implementation Example (Reference)

Below is a Python example that creates a simple circular aperture mask compatible with this specification:

```python
import numpy as np
import tifffile
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = False

def generate_binary_object_mask(shape=(512, 512), type='binary', directory="."):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    filename = f"mask_{type}.tiff"
    save_path = directory / filename

    if type == 'binary':
        mask = np.zeros(shape, dtype=np.uint8)
        cx, cy = shape[1] // 2, shape[0] // 2
        r = min(shape) // 4
        for y in range(shape[0]):
            for x in range(shape[1]):
                if (x - cx) ** 2 + (y - cy) ** 2 < r ** 2:
                    mask[y, x] = 1
        mask_out = mask * 255
        tifffile.imwrite(str(save_path), mask_out, compression=None)

        plt.imshow(mask_out, cmap='gray')
        plt.title("Binary Mask")
        plt.colorbar()
        plt.show()

if __name__ == "__main__":
    generate_binary_object_mask()
```

## Requirements

You must have the following Python packages installed in your environment (e.g., conda environment named `cpi`):

- `numpy`
- `opencv-python` (for general image operations)
- `tifffile` (used to write 16-bit TIFF images)
- `matplotlib` (optional, for visualization/debug)

### Installation (conda recommended)

To ensure compatibility, we recommend using a dedicated conda environment.
You can install all required dependencies with:

```bash
conda install -c conda-forge numpy opencv tifffile matplotlib
```

Make sure to activate the correct environment before running the scripts:

```bash
conda activate cpi
```

