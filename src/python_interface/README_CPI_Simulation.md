# 📘 CPI Simulation Runner – User Guide

This Python script automates running the `propagation.exe` simulator on a collection of input mask files. It is designed for use in batch simulation workflows, where each input mask defines a separate simulation case.

---

## 🧾 Purpose

The simulator requires a configuration file to run. This script takes care of generating those configuration files automatically based on a Python object model. Simulation parameters are exposed as Python variables and must be set appropriately to match the desired physical and experimental setup.

---

## 🗂️ Input: TIFF Masks

To run the simulations, you must first prepare one or more input mask files in `.tiff` format. These files should be placed in a **single directory**, for example:

```
/path/to/masks/
├── mask1.tiff
├── mask2.tiff
└── ...
```

In the script, specify the path to this directory with:

```python
cfg.object.root_dir = "/path/to/masks"
```

The script will loop through all `.tiff` files in this directory and process each one in turn.

---

## 📤 Output: Per-Run Directories

Each `.tiff` file will produce output in its own subdirectory, named after the mask file (without the `.tiff` extension), under a common parent directory.

Example output structure:

```
/path/to/output/
├── mask1/
│   └── [output files]
├── mask2/
│   └── [output files]
└── ...
```

Specify the common output base directory via:

```python
base_storage_dir = "/path/to/output"
```
Note: The directory specified by base_storage_dir must already exist before running the script. It will not be created automatically.

Each simulation’s output path will then be:
```python
cfg.storage.root_dir = base_storage_dir / <mask-name>
```

---

## ⚙️ Configuration Parameters

The simulation parameters (such as wavelength, beam shape, slit geometry, lens settings, etc.) are exposed as Python dataclass fields. These should be configured **explicitly** in the script to suit the desired simulation scenario:

```python
cfg.lambda_ = 632.8e-9
cfg.N = 1264
cfg.illumination.beam_shape = "GAUSSIAN"
cfg.object.w_ratio = 35
...
```

---

## 🚀 Execution Model

- Each simulation is launched as a **separate subprocess**, one at a time.
- The script **waits for each simulation to finish** before starting the next one.
- Standard output and error messages from the simulator are printed to the console.

---

## ▶️ Running the Script

From a terminal or Python environment (e.g., Spyder):

```python
%run cpipy.py
```

Make sure to adjust:
- The path to `propagation.exe`
- The directory with `.tiff` files
- The output directory
- Any relevant simulation parameters
