import numpy as np
import tifffile
from pathlib import Path

def generate_triple_slit_mask(N=1264, w_ratio=35, h_ratio=7, slits=3,
                               w_offset_ratio=0.0, directory="."):

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    save_path = directory / f"triple_slit_mask_{N}.tiff"

    mask = np.zeros((N, N), dtype=np.uint8)
    hnw = N // w_ratio
    hnw -= hnw % 2  # ensure even
    hnh = N // h_ratio
    hnh -= hnh % 2  # ensure even
    onw = int(hnw * w_offset_ratio)
    nw = hnw // slits
    nw -= nw % 2
    nw = max(nw, 2)
    offset = nw // 2

    line = np.zeros(N, dtype=np.uint8)

    if slits % 2 == 1:
        line[N // 2 : N // 2 + nw // 2] = 1
        offset += nw
        slits -= 1

    for i in range(slits // 2):
        start = N // 2 + offset
        line[start : start + nw] = 1
        offset += 2 * nw

    line = line + line[::-1]  # symmetric

    # Apply vertical offset
    if onw != 0:
        temp = np.zeros_like(line)
        if onw > 0:
            temp[onw:] = line[:N - onw]
        else:
            temp[:N + onw] = line[-onw:]
        line = temp

    # Create 2D mask
    mask_section = np.outer(np.ones(hnh, dtype=np.uint8), line)
    mask[(N - hnh)//2 : (N + hnh)//2, :] = mask_section
    mask[:(N - hnh)//2, :] = 0
    mask[(N + hnh)//2:, :] = 0

    mask_out = mask * 255  # binarize
    tifffile.imwrite(str(save_path), mask_out, dtype=np.uint8)
    print(f"[INFO] Corrected triple slit mask saved to: {save_path.resolve()}")

if __name__ == "__main__":
    generate_triple_slit_mask()
