# Assg1_Bio2

## 📖 Description
This is an interactive web-based tool built with **Streamlit** to visualize and perform biological sequence alignment using Dynamic Programming. It supports both **Global Alignment** (Needleman-Wunsch) and **Local Alignment** (Smith-Waterman).

The application calculates the optimal alignment score, generates the scoring matrix, visualizes the traceback path, and produces the final aligned sequences.

## ✨ Key Features
* **Dual Algorithms:**
    * **Needleman-Wunsch:** For aligning two sequences from beginning to end (Global).
    * **Smith-Waterman:** For finding the most similar region within two sequences (Local).
* **Custom Scoring:** Users can define specific values for:
    * Match Score
    * Mismatch Score
    * Gap Penalty
* **Visual Matrices:**
    * Displays the full Scoring Matrix.
    * Displays the Traceback Direction Matrix (using arrows).
* **Path Highlighting:** The optimal path taken by the algorithm is highlighted in **yellow** (and **blue** for the starting point in local alignment) directly on the dataframes.

## 🛠️ Requirements
To run this application, you need Python installed along with the following libraries:

* `streamlit`
* `numpy`
* `pandas`

## 🚀 Installation & Usage

1.  **Clone or Download** this repository/file.

2.  **Install Dependencies**:
    Open your terminal or command prompt and run:
    ```bash
    pip install streamlit numpy pandas
    ```

3.  **Run the Application**:
    Navigate to the folder containing the script (e.g., `app.py`) and run:
    ```bash
    streamlit run app.py
    ```

4.  **Interact**:
    * The application will open in your default web browser (usually at `http://localhost:8501`).
    * Use the **Sidebar** to input your DNA/Protein sequences and adjust scoring parameters.

## 🧬 How It Works

### 1. Matrix Initialization
The system creates a matrix $(m+1) \times (n+1)$ based on the lengths of Sequence 1 and Sequence 2.
* **Global:** The first row and column are initialized with gap penalties.
* **Local:** The first row and column are initialized to 0.

### 2. Scoring (Fill Step)
The matrix is filled using the following logic:
* **Diagonal:** Score + (Match or Mismatch)
* **Up:** Score + Gap Penalty
* **Left:** Score + Gap Penalty
* *(Local Alignment only)*: If a calculated score is negative, it is reset to 0.

### 3. Traceback
* **Global:** Starts from the bottom-right cell and works backward to the top-left (0,0).
* **Local:** Starts from the cell with the **highest score** in the entire matrix and works backward until it hits a cell with a score of 0.
