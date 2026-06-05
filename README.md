# 🧬 PrimersQuest Pro

**PrimersQuest Pro** is a free, web-based qPCR primer design tool built with [Streamlit](https://streamlit.io). It combines the Primer3 thermodynamic engine with Harvard PrimerBank integration, interactive visualisations, and NCBI Primer-BLAST validation links — all in a single Python script.

> Developed by **Dr. Ahmed bey Chaker**, King's College London.

---

## ✨ Features

- **Primer3-powered design** — nearest-neighbour Tm calculations with full thermodynamic salt correction (Mv²⁺, Dv²⁺, dNTP, DNA concentration)
- **Secondary structure checks** — hairpin ΔG, homodimer ΔG, and heterodimer ΔG filtering via `primer3-py`
- **Smart fallback cascade** — relaxes product-size → Tm → GC constraints progressively if no primers are found at first pass
- **Chemistry presets** — SYBR Green and TaqMan (probe-based) parameter sets out of the box
- **RT-qPCR / exon-junction spanning** — automatically extracts CDS splice junctions from GenBank files to force exon-boundary crossing
- **PrimerBank integration** — search Harvard PrimerBank by NCBI gene symbol with canary validation and quality scoring
- **Multi-species support** — Human, Mouse, Rat, and Zebrafish (NCBI taxonomy IDs for Primer-BLAST)
- **Flexible sequence input** — paste raw sequence or FASTA, upload FASTA / GenBank (`.gb`) files, or enter up to 50 kb sequences
- **Interactive visualisations** (Plotly)
  - Genome-browser-style primer map with GC-content track
  - Tm correlation scatter plot
  - Confidence score bar chart
  - Product-size distribution strip
- **Confidence scoring** — composite 0–100 score weighting Tm balance, GC content, and product size
- **NCBI Primer-BLAST links** — one-click validation for every designed primer pair
- **CSV export** — download all designed primer pairs with full metrics
- **SEO-ready** — includes `static/sitemap.xml` and `<meta>` tags for web deployment

---

## 📋 Requirements

- Python ≥ 3.9
- Dependencies listed in `requirements.txt`:

```
streamlit>=1.30.0
pandas>=2.0.0
requests>=2.31.0
beautifulsoup4>=4.12.0
primer3-py>=2.0.0
plotly>=5.18.0
numpy>=1.24.0
biopython>=1.83
```

> **Note:** BioPython is required for GenBank (`.gb`) file parsing and exon-junction extraction. Without it, only FASTA and plain-text sequences are supported.

---

## 🚀 Installation & Running Locally

```bash
# 1. Clone the repository
git clone https://github.com/Gastrodzz/primerquest.git
cd primerquest

# 2. (Recommended) Create a virtual environment
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Launch the app
streamlit run prime6.py
```

The app will open at `http://localhost:8501` in your default browser.

---

## 🖥️ App Structure

The UI is organised into four tabs:

| Tab | Description |
|-----|-------------|
| **Design** | Enter/upload a sequence, choose chemistry preset, set species, and generate primer pairs |
| **PrimerBank** | Search Harvard PrimerBank by gene symbol; results include quality scoring and Primer-BLAST links |
| **Analysis** | Interactive charts — primer map, Tm scatter, confidence bars, product-size strip |
| **About** | Tool information, links to NCBI Primer-BLAST and PrimerBank |

---

## ⚙️ Design Parameters

Default values follow community-accepted qPCR standards and can be adjusted in the sidebar.

### SYBR Green preset

| Parameter | Value |
|-----------|-------|
| Primer length | 18 – 25 bp (opt. 20) |
| Tm range | 58 – 62 °C (opt. 60) |
| GC content | 40 – 60% |
| Product size | 80 – 150 bp |
| Max poly-X run | 3 |
| Primers returned | 10 |

### TaqMan preset

Similar constraints with a tighter product-size window to accommodate an internal probe.

Salt conditions used for all Tm calculations:

| Parameter | Default |
|-----------|---------|
| Monovalent cations (Mv²⁺) | 50 mM |
| Divalent cations (Dv²⁺) | 1.5 mM |
| dNTP concentration | 0.2 mM |
| DNA concentration | 50 nM |

---

## 🔬 Confidence Score

Each primer pair receives a composite **confidence score (0–100%)** based on:

- **Tm balance** — difference between forward and reverse Tm
- **GC content** — proximity to 50% for both primers
- **Product size** — preference for amplicons near 100 bp

| Score | Badge |
|-------|-------|
| ≥ 85% | ✅ EXCELLENT |
| 70 – 84% | ⚠️ GOOD |
| < 70% | ❌ NEEDS OPTIMISATION |

---

## 📂 Repository Structure

```
primerquest/
├── prime6.py          # Main Streamlit application (all logic and UI)
├── requirements.txt   # Python dependencies
└── static/
    └── sitemap.xml    # SEO sitemap for web deployments
```

---

## 🌐 Deploying to Streamlit Community Cloud

1. Fork this repository to your GitHub account.
2. Go to [share.streamlit.io](https://share.streamlit.io) and sign in with GitHub.
3. Select your fork, set the main file to `prime6.py`, and click **Deploy**.

No additional configuration is required — all dependencies are specified in `requirements.txt`.

---

## 🔗 External Resources

- [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) — specificity validation
- [Harvard PrimerBank](https://pga.mgh.harvard.edu/primerbank/) — curated qPCR primer database
- [Primer3](https://primer3.org/) — underlying thermodynamic design engine
- [BioPython](https://biopython.org/) — GenBank file parsing

---

## 📄 License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## 👤 Author

**Dr. Ahmed bey Chaker**  
King's College London

---

*PrimersQuest Pro v3.0 — Advanced qPCR primer design for the modern molecular biology lab.*
