import subprocess
import tempfile
import os
from dataclasses import dataclass, fields, field
from typing import Any, Dict
from pathlib import Path

# === CONFIG ===

@dataclass
class StorageConfig:
    root_dir: str = ""
    overwrite: str = ""
    log_file: str = ""
    armA: str = ""
    armB: str = ""

@dataclass
class TIFFConfig:
    pixel_format: int = 0

@dataclass
class ObjectConfig:
    use_external: str = "Y"
    root_dir: str = ""
    filename: str = ""
    image_crop_factor: float = 1.0

@dataclass
class IlluminationConfig:
    beam_shape: str = ""
    beam_size: float = 0.0
    transversal_coherence_length: float = 0.0

@dataclass
class LensConfig:
    focal_length: float = 0.0
    aperture_diameter: float = 0.0

@dataclass
class CommonArmConfig:
    object_to_lens: float = 0.0
    lens_to_detectorA: float = 0.0
    lens_to_detectorB: float = 0.0

@dataclass
class CPIConfig:
    storage: StorageConfig = field(default_factory=StorageConfig)
    tiff: TIFFConfig = field(default_factory=TIFFConfig)
    experiment: str = ""
    frames: int = 0
    N: int = 0
    side_length_in_meter: float = 0.0
    lambda_: float = 0.0
    illumination: IlluminationConfig = field(default_factory=IlluminationConfig)
    object: ObjectConfig = field(default_factory=ObjectConfig)
    lens: LensConfig = field(default_factory=LensConfig)
    common_arm: CommonArmConfig = field(default_factory=CommonArmConfig)

    def to_flat_dict(self) -> Dict[str, Any]:
        def flatten(prefix: str, obj: Any) -> Dict[str, Any]:
            flat = {}
            if hasattr(obj, "__dataclass_fields__"):
                for f in fields(obj):
                    val = getattr(obj, f.name)
                    key = f"{prefix}.{f.name}" if prefix else f.name
                    flat.update(flatten(key, val))
            else:
                flat[prefix] = obj
            return flat

        flat_dict = flatten("", self)
        if "lambda_" in flat_dict:
            flat_dict["lambda"] = flat_dict.pop("lambda_")
        return flat_dict

    def to_cfg(self) -> str:
        flat = self.to_flat_dict()
        return "\n".join(f"{k} = {v}" for k, v in flat.items())


# === RUN LOGIC ===

def run_cpi(cfg: CPIConfig, exe_path: Path):
    # Scrive file cfg temporaneo
    cfg_text = cfg.to_cfg()
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.cfg', delete=False) as f:
        f.write(cfg_text)
        f.flush()
        os.fsync(f.fileno())
        cfg_path = f.name

    print(f"[INFO] Lanciando: {exe_path} {cfg_path}")

    env = os.environ.copy()
    env["PATH"] = r"C:\msys64\ucrt64\bin;" + env["PATH"]

    command = f'"{exe_path}" "{cfg_path}"'
    result = subprocess.run(
        command,
        shell=True,
        env=env,
        capture_output=True,
        text=True
    )
    print("--- STDOUT ---")
    print(result.stdout)
    print("--- STDERR ---")
    print(result.stderr)

    return result.returncode


# === MAIN ===

if __name__ == "__main__":
    PROPAGATION_EXE_PATH = Path(
        "/Users/sdegi/home/pers/fisica/ricerca/cpisimulation/dev/cpi_root/build/bin/propagation/propagation.exe"
    )

    # Configurazione base
    base_storage_dir = "/Users/sdegi/home/pers/fisica/ricerca/cpisimulation/dev/cpi_root/analysis/data/python_interface"

    cfg = CPIConfig()
    cfg.storage.overwrite = "Y"
    cfg.storage.log_file = "cpi.log"
    cfg.storage.armA = "framesA"
    cfg.storage.armB = "framesB"

    cfg.tiff.pixel_format = 16

    cfg.experiment = "CPI"
    cfg.frames = 1
    cfg.N = 1264
    cfg.side_length_in_meter = 4e-3
    cfg.lambda_ = 632.8e-9

    cfg.illumination.beam_shape = "GAUSSIAN"
    cfg.illumination.beam_size = 0.125e-3
    cfg.illumination.transversal_coherence_length = 3.165e-6

    cfg.object.use_external = "Y"  # Imposto a Y sempre
    cfg.object.root_dir = "/Users/sdegi/home/pers/fisica/ricerca/cpisimulation/dev/cpi_root/python/cpipy_root/cpipy"
    cfg.object.image_crop_factor = 0.02857142857142857142857142857143

    cfg.lens.focal_length = 10.0e-3
    cfg.lens.aperture_diameter = 2e-3

    cfg.common_arm.object_to_lens = 10.0e-3
    cfg.common_arm.lens_to_detectorA = 12.0e-3
    cfg.common_arm.lens_to_detectorB = 0.0e-3

    # Loop sui file .tiff in object.root_dir
    mask_dir = Path(cfg.object.root_dir)
    for tiff_file in sorted(mask_dir.glob("*.tiff")):
        cfg.object.filename = tiff_file.name
        sub_output_dir = Path(base_storage_dir) / tiff_file.stem
        cfg.storage.root_dir = str(sub_output_dir)

        print(f"\n=== PROCESSING {tiff_file.name} ===")
        run_cpi(cfg, PROPAGATION_EXE_PATH)
