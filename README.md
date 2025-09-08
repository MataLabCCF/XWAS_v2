# XWAS Version 2
## Overview

This pipeline is an updated version of the XWAS pipeline described in Leal et al., 2023
(DOI: 10.1002/mds.29508), originally available at https://github.com/MataLabCCF/XWAS

In the original repository, the XWAS analysis was performed in three main steps:
- harmonizationAdapted.py
- makeAllStepsToRegression.py (or the autosomal PCA-based version: makeAllStepsToRegressionAutosomalPCA.py)
- metaAnalysisGWAMA.py

## What's New in Version 2

In this updated pipeline, we kept the original harmonizationAdapted.py unchanged, but we merged and improved the 
functionalities of makeAllStepsToRegression.py and metaAnalysisGWAMA.py.

Key updates include:
- A simplified workflow by combining multiple steps into a single process.
- A new configuration file that allows you to easily specify all analysis parameters.
- Increased flexibility based on valuable feedback from collaborators who found the previous implementation too rigid.

We sincerely thank our collaborators for their input, which was essential for making this tool more user-friendly and adaptable.

---

## Configuration File Overview

The **Configuration File** is the central component of the XWAS pipeline.  
It defines **all the parameters** required to perform the analysis, including:

- Input data sources  
- Reference datasets  
- Outlier detection parameters  
- PCA settings  
- Regression model configuration  
- Paths to external tools  

Each line in the configuration file represents **one instruction** for the pipeline.  
The file follows a **tab-delimited format** and is composed of **flags** and their respective parameters.

### General Rules
- **Comments** → Any line starting with `#` will be **ignored** by the pipeline.  
- **Order of flags** → The flags can be declared **in any order**.  
- **Case sensitivity** → Flags and keywords (e.g., `Input`, `Reference`, `covar`) are not **case-sensitive** (I think).  
- **Paths** → Always provide **absolute paths** to files and tools to avoid unexpected errors.  
- **Consistency** → Make sure chromosome types (`X` or `Autosomal`) and file formats match the required specifications.


## How to Configure Your Analysis

All the configuration for your analysis is done through the **Configuration File**.  
It is composed of six main flags: **Input**, **Reference**, **Outlier**, **PCA**, **Model**, and **Programs**.

> **Note:** Any line starting with `#` is treated as a **comment** and will be ignored by the code.

### **1. Input**

The **Input** section specifies the data required for the analysis.  
It can be divided into two types: **genetic data** and **non-genetic data**.

#### **1.1. Input for Genetic Data**

To declare input genetic data, you need to provide four columns:

| Column | Description                           | Example                              |
|--------|-------------------------------------|--------------------------------------|
| 1      | **Flag** (always `Input`)           | `Input`                             |
| 2      | **Chromosome type** (`X` or `Autosomal`) | `X`                          |
| 3      | **Data type** (`Typed` or `Imputed`) | `Typed`                            |
| 4      | **Path to file**                    | `/home/path/target/genotypedX.vcf.gz` |

**Example:**

```
Input   X   Typed   /home/path/target/genotypedX.vcf.gz 
Input   Autosomal   Typed   /home/path/target/genotypedAuto.vcf.gz
Input   X   Imputed /home/path/target/imputedX.vcf.gz
```

> **Important:** The **Autosomal Typed** data will be used **only for PCA purposes**.  
> **Imputed autosomal data** is **not expected** and should **not** be provided.

#### **1.2. Input for Non-Genetic Data**

You can also provide a table containing **demographic or clinical variables** (e.g., `AGE`, `SEX`, etc.).  
For this input, the columns should be:

| Column | Description                  |
|--------|----------------------------|
| 1      | **Flag** (always `Input`)  |
| 2      | Must be set to `covar`     |
| 3      | Path to the **covariates file** |
| 4      | Name of the **dependent variable** |

**Example:**
```
Input   covar   COVAR.txt   disease
```

### **2. Reference**

This section provides the **reference data** used for **projected PCA**.  
It consists of three columns:

| Column | Description                           | Example                         |
|--------|---------------------------------------|---------------------------------|
| 1      | **Flag** (always `Reference`)         | `Reference`                    |
| 2      | **Chromosome type** (`X` or `Autosomal`) | `X`                         |
| 3      | **Path to file**                      | `/home/path/ref/refX.vcf.gz`  |

**Example:**

```
Reference	X	/home/path/ref/refX.vcf.gz
Reference	Autosomal	/home/path/ref/refAuto.vcf.gz
```

### **3. Outlier**

This section is used to select which PCA(s) should be used to detect **genetic outliers**.  
It consists of **four columns**:

| Column | Description                            | Example   |
|--------|----------------------------------------|-----------|
| 1      | **Flag** (always `Outlier`)            | `Outlier` |
| 2      | **Chromosome type** (`X` or `Autosomal`) | `X`    |
| 3      | **Projected PCA?** (`Yes` or `No`)     | `Yes`    |

**Example:**

```
Outlier	X	No
Outlier	X	Yes
Outlier	Autosomal	No
Outlier Autosomal	Yes
```

To declare the number of PCs to be used for **X-PCA** and **Autosomal PCA**,  
add additional lines using the `Outlier` flag following this pattern:

| Column | Description                          | Example   |
|--------|--------------------------------------|-----------|
| 1      | **Flag** (always `Outlier`)          | `Outlier` |
| 2      | **PC flag** (always `PC`)            | `PC`      |
| 3      | **Chromosome type** (`X` or `Autosomal`) | `X`   |
| 4      | **Number of PCs to use**            | `2`       |

**Example:**
```
Outlier PC	X	5
Outlier PC	Autosomal	10
```

> **Note:** All inferences described in this section will be executed. In the example we will run the 4 PCs to detect 
> outliers

### **4. PCA**

This section selects which **PCA** will be used in the regression analysis.  
The pipeline was designed for **LARGE-PD data**, which has a high **Native American** component.  
Because **projected PCA** is **not suitable** in datasets with **low Native American representation**,  
this pipeline does **not** allow projected PCA to be used for regression.

> If you need projected PCA for regression, please contact the developer.

This section consists of two columns:

| Column | Description                           | Example   |
|--------|-------------------------------------|-----------|
| 1      | **Flag** (always `PCA`)            | `PCA`     |
| 2      | **Chromosome type** (`X` or `Autosomal`) | `X` |

### **5. Model**

This section defines the **regression model** to be used or provides the path to `selectModel.R`  
to perform a **stepwise regression** for model selection.

| Column            | Description                            | Example                            |
|-------------------|----------------------------------------|------------------------------------|
| 1                 | **Flag** (always `Model`)             | `Model`                            |
| 2                 | **Script or regression formula**      | `Script`                           |
| 3 (if script)     | **Path to selectModel.R**             | `/home/path/scripts/selectModel.R` |
| 3 (if regression) | **Model formula**                     | `PHENO~SEX+PC1+PC2`                |

### **6. Programs**

This section tells the pipeline **how to call the external tools** used in the analysis.  
The pipeline relies on **Plink 1**, **Plink 2**, **GCTA**, **GWAMA**, and **RScript**.  
This section consists of three columns:

| Column | Description                     | Example                |
|--------|---------------------------------|------------------------|
| 1      | **Flag** (always `Programs`)    | `Programs`            |
| 2      | **Tool name**                   | `plink`               |
| 3      | **Path or command to call tool** | `/home/programs/plink` |



## Command line

Because of the Config file, the command line is simple:

```
usage: main.py [-h] -c CONFIG [-t THREADS] -f FOLDER

Mata Lab XWAS v2.1

options:
  -h, --help            show this help message and exit

Input arguments:
  -c CONFIG, --config CONFIG
                        Path to configuration file
  -t THREADS, --threads THREADS
                        Number of threads to be used by GCTA (default = 1)

Output arguments:
  -f FOLDER, --folder FOLDER
```

**Example:**

```commandline
python main.py -c Config.txt -t 4 -f /home/path/XWAS_results
```

 ## Acknowledgements

This work is supported by NIH Grant R01 1R01NS112499-01A1, MJFF Grant ID: 18298, ASAP-GP2 and Parkinson's Foundation

 
### Contact
 
 Developer: Thiago Peixoto Leal. PhD (PEIXOTT@ccf.org or thpeixotol@hotmail.com)