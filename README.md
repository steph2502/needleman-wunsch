# Needleman-Wunsch Sequence Analyzer

A full-stack (local) web application for **global sequence alignment** using the **Needleman–Wunsch** algorithm.

- **Backend**: Flask (Python)
- **Frontend**: Plain HTML/CSS/JavaScript
- **Styling**: Bootstrap 5 + custom CSS

## Features

- Global alignment for **DNA**, **RNA**, and **Protein** sequences
- Configurable scoring: match, mismatch penalty, gap penalty
- Accepts plain sequences or **FASTA** input (FASTA headers are ignored)
- Optional support for **ambiguity codes** (toggle in the UI)
- Alignment output with match markers
- **Scoring matrix** visualization with highlighted **traceback path** and final optimal score
- Metrics: matches, mismatches, gaps, identity %
- Copy alignment + download as TXT
- Reset form, loading spinner, dark/light mode toggle

## Run locally

Install dependencies:

```bash
pip install -r requirements.txt
```

Start the server:

```bash
python app.py
```

Then open the URL shown in your terminal (usually `http://127.0.0.1:5000`).

## Project structure

```
project/
├── app.py
├── requirements.txt
├── README.md
├── utils/
│   └── needleman_wunsch.py
├── templates/
│   ├── base.html
│   ├── index.html
│   ├── results.html
│   └── about.html
└── static/
    ├── css/
    │   └── style.css
    └── js/
        └── script.js
```

