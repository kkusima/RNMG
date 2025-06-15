# RNMG
RNMG - Reaction Network Matrix Generator

**A Streamlit app to parse chemical reaction mechanisms and automatically generate:

- **Stoichiometric Matrix** (species vs. reactions)
- **Atomic Matrix** (atom counts per species)
- **Parameter Table** (editable kinetics parameters)
- **Mass‑balance checker** to verify atom conservation

---

## 🚀 Features

- Input arbitrary reactions using `⇌`, `->`, `<->`, etc.
- Supports surface sites (`*`) and gas‑phase species (e.g. `CO(g)`)
- Auto‑orders columns: **Gas → Surface → Empty Sites (*)**
- Downloadable CSVs for matrices and parameter tables
- Consistency check for species ordering
- Mass‑balance verification via \(A \times v = 0\)

---

## 🛠️ Prerequisites

- Python 3.8+
- Streamlit
- pandas
- NumPy

---

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
streamlit run reaction_network_app.py
```

This will open a browser window (or show a local URL) where you can:

1. Enter reactions in the sidebar (e.g. `CO(g) + * ⇌ CO*`).
2. View/download the **Stoichiometric Matrix**, **Atomic Matrix**, and **Parameter Table**.
3. Run the **Mass‑Balance Check** to ensure atom conservation.

---

## 📝 Usage Tips

- **Arrow notation**: you can use `⇌`, `↔`, `<->`, or `->`.
- **Surface sites**: denote with `*`.
- **Coefficients**: prefix species with a number (e.g. `2O*`).
- **Phase labels** (optional): `(g)`, `(l)`, `(s)` are stripped when computing atomic composition.
- To clear all reactions, click **Clear All** in the sidebar.

---

## 🤝 Contributing

1. Fork the repo.
2. Create a feature branch: `git checkout -b feature/YourFeature`
3. Commit your changes: `git commit -m "Add YourFeature"`
4. Push to your branch: `git push origin feature/YourFeature`
5. Open a pull request.

---

## 📄 License

This project is licensed under the MIT License.

---

## © 2025 Kenneth Kusima
**
