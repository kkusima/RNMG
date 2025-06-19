# RNMG - Reaction Network Matrix Generator

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://rnmgapp.streamlit.app)

A comprehensive Streamlit application for parsing chemical reaction mechanisms and automatically generating stoichiometric matrices, atomic matrices, and thermodynamic analysis for heterogeneous catalysis research. RNMG features an intuitive welcome screen interface with dedicated modes for matrix generation and file analysis.

## 🚀 Features

### **Interactive Welcome Screen**
- **Comprehensive Introduction** - Detailed overview of RNMG capabilities and features
- **Mode Selection Interface** - Easy navigation between Generator and Upload modes
- **Quick Start Buttons** - Direct access to working modes from the welcome screen
- **Visual Feature Comparison** - Side-by-side comparison of Generator vs Upload modes

### **Core Functionality**
- **Stoichiometric Matrix Generation** - Species vs. reactions with proper coefficient handling
- **Atomic Matrix Creation** - Atomic composition tracking for all species
- **Parameter Management** - Editable kinetic and thermodynamic parameters
- **Mass Balance Verification** - Automated atom conservation checking using A × ν = 0
- **Equilibrium Analysis** - Calculate equilibrium constants and Gibbs free energy changes

### **Advanced Capabilities**
- **Multi-Mode Interface** - Welcome screen, Generator mode, and Upload mode
- **Smart Species Ordering** - Auto-orders species: **Gas Phase → Surface Species → Empty Sites**
- **Multiple Surface Sites** - Supports various surface site types (`*`, `_`, `#`, `@`, `&`)
- **Complex Formula Parsing** - Handles formulas with parentheses (e.g., Ca(OH)₂, Al₂(SO₄)₃)
- **Comprehensive Analysis** - Thermodynamic analysis with interactive visual plots
- **Matrix Consistency Verification** - Ensures data integrity across matrices
- **CSV Export Functionality** - Download all generated data for external use

## 🛠️ Prerequisites

- Python 3.8+
- Required packages:
  - streamlit
  - pandas
  - numpy
  - plotly
  - re (built-in)

## 📦 Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/kkusima/RNMG.git
   cd RNMG
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


This opens RNMG in your browser with the **Welcome Screen** where you can:

1. **Learn about RNMG** - Comprehensive introduction to features and capabilities
2. **Choose your mode** - Select Generator Mode or Upload Mode based on your needs
3. **Get quick help** - Access mode-specific guidance and workflows

## 🎯 Application Modes

### **🏠 Welcome Mode**
The default landing screen providing:
- **Feature Overview** - Detailed comparison of Generator vs Upload modes
- **Getting Started Guide** - Step-by-step workflows and best practices
- **Important Guidelines** - Rules, limitations, and supported formats
- **Quick Navigation** - Direct buttons to switch to working modes

### **🔬 Generator Mode**
Build reaction networks from scratch:
1. **Enter reactions** using ⇌ arrows (e.g., `CO(g) + * ⇌ CO*`)
2. **Auto-generate matrices** - Stoichiometric and atomic matrices
3. **Create parameters** - Default kinetic and thermodynamic parameters
4. **Perform analysis** - Mass balance and equilibrium analysis
5. **Export results** - Download CSV files for external modeling

### **📁 Upload Mode**
Analyze existing data:
1. **Upload CSV matrices** - Import atomic and stoichiometric matrices
2. **Upload parameters** - Import kinetic parameter files
3. **Run comprehensive analysis** - Mass balance and thermodynamic analysis
4. **Visualize results** - Interactive plots and detailed reports

## 📊 Analysis Capabilities

### **Mass Balance Analysis**
Verifies atom conservation using the fundamental equation:

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

RNMG is specifically designed for **heterogeneous catalysis research**:

- **Surface Reaction Mechanisms** - CO oxidation, water-gas shift, methane reforming
- **Catalyst Active Site Modeling** - Multiple site types and surface coverage effects
- **Thermodynamic Studies** - Reaction feasibility and equilibrium analysis
- **Mass Balance Verification** - Ensuring atom conservation in complex networks
- **Kinetic Parameter Organization** - Structured data for microkinetic modeling

## 📄 File Formats and Export

### **Export Options**
All matrices and parameters can be exported as CSV files:

**1. Atomic Matrix** (`atomic_matrix.csv`)
- Rows: Chemical elements and surface sites
- Columns: All species in the reaction network
- Values: Stoichiometric coefficients for atomic composition

**2. Stoichiometric Matrix** (`stoichiometric_matrix.csv`)
- Rows: Reactions (r1, r2, r3, ...)
- Columns: Species with proper notation (P_species for gas, theta_species for surface)
- Values: Net stoichiometric coefficients

**3. Parameters** (`parameters.csv`)
- Temperature, gas constant, pressures
- Forward/reverse rate constants
- Additional modeling constants

### **Import Formats**
Upload Mode accepts CSV files with the same structure as the exported files.

## ⚠️ Current Limitations

- **Reaction Types**: Optimized for gas-solid (heterogeneous) interfaces
- **Surface Sites**: Supports 5 different site types (`*`, `_`, `#`, `@`, `&`)
- **Mechanism Scope**: Designed for elementary surface reactions
- **Phase Systems**: Currently limited to gas-surface interactions

**Future Development:**
- Liquid-phase reaction support
- Multi-phase system capabilities
- Complex reaction intermediate handling

## 📚 Version History

- **v1.1**: Added multi-site capability, thermodynamic analysis
- **v1.0**: Initial release with basic matrix generation and mass balance checking

## 📖 Citation

If you use RNMG in your research, please cite:

### **BibTeX**

```@software{kusima2025rnmg,
author = {Kenneth Kusima},
title = {RNMG - Reaction Network Matrix Generator: A Comprehensive Tool for Heterogeneous Catalysis Research},
year = {2025},
version = {1.1},
url = {https://github.com/kkusima/RNMG},
note = {Streamlit application for stoichiometric and atomic matrix generation from chemical reaction mechanisms}
}
```
### **APA Format**
Kusima, K. (2025). *RNMG - Reaction Network Matrix Generator* (Version 1.1) [Computer software]. https://github.com/kkusima/RNMG

## 🤝 Contributing

Contributions are welcome! Please feel free to:

- **Submit Issues** - Report bugs or request features
- **Feature Requests** - Suggest improvements for catalysis research
- **Pull Requests** - Contribute code improvements
- **Documentation** - Help improve user guides and examples

## 📞 Support and Contact

- **GitHub Issues**: [Report bugs or request features](https://github.com/kkusima/RNMG/issues)
- **GitHub Profile**: [https://github.com/kkusima](https://github.com/kkusima)
- **Email**: [mailto:kennethlcs5+githubRNMG@gmail.com](kennethlcs5+githubRNMG@gmail.com)

## 👨‍💻 Author

**Kenneth Kusima** - Chemical Engineering Researcher
- Specialized in heterogeneous catalysis and reaction network analysis
- Contact: [GitHub Profile](https://github.com/kkusima)

## 📄 License
This project is licensed under the MIT License.

---

**© 2025 Kenneth Kusima** | Built using Streamlit
