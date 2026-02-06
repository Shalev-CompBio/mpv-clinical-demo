# How to Run the MPV Clinical Decision Support App

Follow these simple steps to run the application on your computer.

## Step 0: Open Terminal

1. Press `Win` + `R` on your keyboard.
2. Type `powershell` and press **Enter**.
3. A blue/black window will open.

## Step 1: Navigate to the Project Folder

Copy and paste this exact command into the window and press **Enter**:

```powershell
cd "c:\Users\Shalev\OneDrive - huji.ac.il\bio_projects\multiomics_IRD_genes_clustering\IRD_phenotype_genes_network\MPV"
```

## Step 2: Install Requirements (One time only)

If this is your first time running it, install the necessary libraries by pasting this command and pressing **Enter**:

```powershell
pip install -r requirements.txt
```

*Wait for it to finish installing.*

## Step 3: Run the App

Paste this command and press **Enter**:

```powershell
streamlit run app.py
```

## Step 4: Open in Browser

- The browser should open automatically.
- If not, open your web browser (Chrome/Edge) and type:
  **http://localhost:8501**

---

## Troubleshooting

- **"Port is already in use"**: If 8501 is taken, Streamlit will use the next available port (e.g., 8502 or 8503). Look at the text in the terminal window to see the correct URL.
- **To stop the app**: Go back to the PowerShell window and press `Ctrl` + `C`.
