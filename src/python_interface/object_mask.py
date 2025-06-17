import numpy as np
import tifffile
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = False

def generate_mask(shape=(512, 512), type='binary', directory="."):
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

    elif type == 'apodized':
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
    else:
        raise ValueError("Unknown mask type. Use 'binary' or 'apodized'.")

    print(f"[INFO] Mask saved to: {save_path.resolve()}")

if __name__ == "__main__":
    generate_mask( shape=(1264,1264))
