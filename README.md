# RNMG - Reaction Network Matrix Generator

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://rnmgapp.streamlit.app)

A comprehensive Streamlit application for parsing chemical reaction mechanisms and automatically generating stoichiometric matrices, atomic matrices, and thermodynamic analysis for heterogeneous catalysis research.

## 🚀 Features

**Core Functionality:**
- **Stoichiometric Matrix Generation** - Species vs. reactions with proper coefficient handling
- **Atomic Matrix Creation** - Atomic composition tracking for all species
- **Parameter Management** - Editable kinetic and thermodynamic parameters
- **Mass Balance Verification** - Automated atom conservation checking using A × ν = 0
- **Equilibrium Analysis** - Calculate equilibrium constants and Gibbs free energy changes

**Advanced Capabilities:**
- Supports multiple surface site types (`*`, `_`, `#`, `@`, `&`)
- Handles complex chemical formulas with parentheses (e.g., Ca(OH)₂, Al₂(SO₄)₃)
- Auto-orders species: **Gas Phase → Surface Species → Empty Sites**
- Comprehensive thermodynamic analysis with visual plots
- Matrix consistency verification
- CSV export functionality for all generated data

## 🛠️ Prerequisites

- Python 3.8+
- Required packages:
  - streamlit
  - pandas
  - numpy
  - plotly
  - re (built-in)

## 📦 Installation

1. Clone this repo:
   ```bash
   git clone https://github.com/kkusima/RNMG.git
   cd reaction_network_app
   ```
2. Create and activate a virtual environment (optional but recommended):
   ```bash
   python -m venv .venv
   source .venv/bin/activate     # macOS/Linux
   .venv\Scripts\activate      # Windows
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

---
## ▶️ Running the App

```bash
streamlit run RNMG_app.py
```

This will open a browser window (or show a local URL) where you can:

1. Enter reactions in the sidebar (e.g. `CO(g) + * ⇌ CO*`).
2. View/download the **Stoichiometric Matrix**, **Atomic Matrix**, and **Parameter Table**.
3. Run the **Mass‑Balance Check** to ensure atom conservation.

---

### Upload Mode
1. **Upload existing matrices** (CSV format)
2. **Upload parameter files** for analysis
3. **Perform comprehensive analysis** on uploaded data

## 📊 Analysis Types

### Mass Balance Analysis
Verifies atom conservation across all reactions using the fundamental equation:

A × ν = 0
Where A is the atomic matrix and ν is the stoichiometric vector.

### Equilibrium Analysis
Calculates equilibrium constants and thermodynamic properties:
- Individual equilibrium constants (K_eq = k_forward / k_reverse)
- Overall equilibrium constant (K_overall = ∏K_i^σ_i)
- Gibbs free energy changes (ΔG = -RT ln K)
- Visual analysis with interactive plots

## 📝 Input Guidelines

**Reaction Format:**
- Use **only ⇌ arrows** (copy from the input box)
- Enter **one reaction at a time**
- Coefficients: prefix species with numbers (e.g., `2O*`, `3H2O(g)`)

**Supported Species:**
- **Gas phase**: `CO(g)`, `H2O(g)`, `O2(g)` etc.
- **Surface species**: `CO*`, `H2O*`, `OH*` etc.
- **Surface sites**: `*`, `_`, `#`, `@`, `&`
- **Complex formulas**: `Ca(OH)2`, `Al2(SO4)3` etc.

**Phase Labels** (optional): `(g)`, `(l)`, `(s)` are automatically handled

## 🔬 Scientific Applications

RNMG is specifically designed for **heterogeneous catalysis research** and supports:
- Surface reaction mechanism analysis
- Catalyst active site modeling
- Thermodynamic feasibility studies
- Mass balance verification for complex reaction networks
- Kinetic parameter organization

## 📄 File Formats

**Export Options:**
- Atomic Matrix: CSV format with atoms as rows, species as columns
- Stoichiometric Matrix: CSV format with reactions as rows, species as columns
- Parameters: CSV format with reaction descriptions, parameter names, values, and units

## ⚠️ Limitations

- **Generator Mode**: Limited to gas-solid (heterogeneous) interfaces
- **Surface Sites**: Currently supports 5 different site types
- **Reaction Types**: Designed primarily for elementary surface reactions

## 📚 Version History

- **v1.1**: Added multi-site capability, thermodynamic analysis
- **v1.0**: Initial release with basic matrix generation and mass balance checking

## 📄 License

This project is licensed under the MIT License.

## 👨‍💻 Author

**Kenneth Kusima** - Chemical Engineering Researcher
- Specialized in heterogeneous catalysis and reaction network analysis
- Contact: [GitHub Profile](https://github.com/kkusima)

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests to improve RNMG's functionality for the catalysis research community.

---

**© 2025 Kenneth Kusima** | Built using Streamlit
