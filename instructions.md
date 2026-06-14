# Project: Spinodoid Topology Optimization Framework
**Status:** Active Porting (MATLAB -> Python)
**License:** MIT

## 1. Overarching Goal
We are porting a multi-scale topology optimization framework (focusing on spinodoid metamaterials) from a legacy MATLAB/Gibbon/FEBio pipeline into a fully open-source, high-performance Python package. The final product will feature a PyTorch surrogate model, a zero-crash SciPy sparse finite element solver, and an interactive Streamlit/PyVista 3D UI.

## 2. Rigid Repository Structure
Do not generate or place files in the root directory. Adhere strictly to this layout:
- `matlab_framework/`: Contains all legacy MATLAB scripts, datasets, and exported intermediate files.
- `python_framework/.venv/`: The active virtual environment (pre-installed with PyTorch, SciPy, etc.).
- `python_framework/src/`: All generated Python packages and modules go here.
- `python_framework/tests/`: All `pytest` verification suites go here.

## 3. Data & File Interoperability Rules
- **Neural Network Weights:** MATLAB proprietary objects cannot be read by Python. Read trained weights exclusively from the flattened dictionaries inside `matlab_framework/exported_weights.mat` using `scipy.io.loadmat`.
- **Testing:** All Python math and PyTorch predictions must be verified against MATLAB baseline data to a precision of 4 decimal places via automated `pytest` suites.

## 4. Agent Instructions (Auto-Updating)
Claude, as you operate in this repository:
1. Always activate the `python_framework/.venv/` environment before running tests or scripts.
2. Read this file at the start of every session to establish context.
3. When you complete a sub-task below, you are explicitly instructed to edit this `CLAUDE.md` file, change the `[ ]` to `[x]`, and add a 1-sentence note on the implementation.

## 5. Master Task Tracker
*Agent Instruction: You are encouraged to dynamically add new `[ ]` subtasks under these main phases as you encounter necessary architectural steps or debugging requirements. Ensure UI components are built modularly so they can be combined into a single master application later.*

### Phase 0: Project Setup & Foundation
- [x] Establish repository structure and MIT license.
- [x] Set up Python virtual environment and pre-install heavy dependencies.

### Phase 1: Neural Network Surrogate Migration
- [x] Export MATLAB neural network parameters to a flat `.mat` dictionary.
- [ ] Implement PyTorch neural network architecture in `python_framework/src/`.
- [ ] Load `exported_weights.mat` into the PyTorch `state_dict` and pass the 4-decimal `pytest` suite.
- [ ] **UI Task:** Build a Streamlit data explorer module to visualize the homogenized spinodoid datasets and compare PyTorch property predictions against target values.

### Phase 2: Topology Optimization Engine
- [ ] Translate the core topology optimization loop (e.g., OC/MMA) from MATLAB to Python.
- [ ] Integrate the PyTorch surrogate model to evaluate localized material properties during the optimization loop.
- [ ] **UI Task:** Build a Streamlit dashboard module that plots the objective function convergence and volume fraction constraints in real-time as the loop runs.

### Phase 3: High-Performance FEM Integration
- [ ] Audit the legacy dense FEM matrix assembly that caused memory crashes on large grids.
- [ ] Implement a memory-efficient sparse solver pipeline using `scipy.sparse` for the macro-scale topology optimization evaluations.
- [ ] **UI Task:** Add a module to the dashboard visualizing the live macro-scale density distribution updates during the FEM optimization steps.

### Phase 4: Spinodoid Rendering & Export
- [ ] Translate the spinodoid phase-field generation and micro-structure rendering mathematics to vectorized Python arrays.
- [ ] Build an automated export pipeline to generate clean, production-ready STL/OBJ files.
- [ ] **UI Task:** Integrate PyVista into the Streamlit app to render interactive 3D structures, allowing users to click macro-elements to preview the generated micro-scale spinodoid geometry.

### Phase 5: Native Python Meshing & Verification 
- [ ] Phase out manual Gibbon steps by integrating automated Python-native meshing tools (e.g., PyGmsh or Trimesh).
- [ ] Implement a native Python Finite Element Modelling framework (e.g., `SfePy`) for detailed structural verification and stress analysis.
- [ ] **FEBio Fallback:** If native Python solvers fail to converge on highly non-linear metamaterial verification, default back to FEBio. Build a Python wrapper to automatically write `.feb` XML files and call the FEBio executable in the background. 
  - *Note for Claude:* The path to the FEBio executable must be parameterized at the top of the verification script so it can be easily updated. (e.g., `FEBIO_PATH = "/mnt/c/Program Files/FEBioStudio/bin/FEBio4.exe"` or equivalent Linux path).
- [ ] **UI Task:** Create the final "Verification" tab in the UI that visualizes the high-fidelity von Mises stress distributions on the fully meshed, final optimized structure.

## 6. Environment Library Stack
The active virtual environment (`python_framework/.venv/`) contains the following pre-installed core libraries. Use these exclusively before requesting new dependencies:
- **`torch` (PyTorch):** Machine learning surrogate model generation and evaluation.
- **`scipy` (`scipy.sparse`):** Core dense/sparse matrix operations and linear solvers for the topology optimization loop.
- **`numpy` & `pandas`:** Data ingestion, vectorization, and array manipulation.
- **`pytest`:** Automated mathematical and structural verification suites.
- **`streamlit`:** Interactive web-based user interface and dashboard generation.
- **`pyvista` & `matplotlib`:** 3D and 2D data visualization, stress plotting, and rendering.
- **`trimesh`:** 3D geometry manipulation, isosurface extraction, and STL/OBJ export.
- **`sfepy`:** Simple Finite Elements in Python for native non-linear structural verification.