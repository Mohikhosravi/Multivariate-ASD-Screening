# Multivariate Analysis of ASD Adult Screening Data

This project applies a comprehensive set of multivariate statistical methods to the **Autism Spectrum Disorder (ASD) Adult Screening Dataset**, exploring the structure of ASD-related behavioral features and evaluating their capacity to discriminate between ASD-positive and ASD-negative individuals.

The analysis was conducted as part of a graduate-level **Multivariate Analysis** course, with the goal of moving beyond univariate summaries and leveraging the joint distribution of screening items to uncover latent structure and group differences.

---

## Research Questions

1. Do the ASD screening items (A1–A10) and age jointly differ between individuals diagnosed with ASD and those without?
2. What is the underlying factor structure of the screening instrument — do the items cluster into interpretable latent dimensions?
3. How much of the total variance can be captured by a reduced set of principal components, and which items contribute most?
4. Can a linear discriminant function accurately classify individuals into ASD / non-ASD groups using these variables alone?

---

## Dataset

| Property | Detail |
|---|---|
| Source | [Autism Adult Screening Data — UCI / Kaggle](https://raw.githubusercontent.com/shaheennamboori/CSV_dataset_for_autism_diagnostics/master/Autism-Adult_Data.csv) |
| Observations | 704 adults (after removing missing age values) |
| Variables | 21 total; 11 used in analysis (A1–A10 scores + age) |
| Target variable | `Class/ASD` — binary (YES / NO) |

### Key Variables

- **A1–A10 Score**: Binary responses (0/1) to the Autism Spectrum Quotient (AQ-10) screening questionnaire items
- **age**: Participant age in years
- **Class/ASD**: Clinical ASD diagnosis (YES = ASD-positive)

> **Note:** The `result` variable (sum of A1–A10) was excluded from all analyses because it is a perfect linear combination of the item scores, causing a singular covariance matrix and making multivariate inference impossible.

---

## Methods

### 1. Multivariate Descriptive Analysis
Computed the mean vector, covariance matrix, and correlation matrix for all numeric variables. Visualized pairwise relationships using a scatterplot matrix colored by ASD class, and identified multivariate outliers using **Mahalanobis distance** with a chi-square cutoff at the 97.5th percentile.

### 2. Multivariate Normality — Mardia's Test
Assessed multivariate normality using **Mardia's skewness and kurtosis statistics**, implemented manually without reliance on a package. Complemented with a chi-square Q-Q plot of Mahalanobis distances. This step was necessary to evaluate the validity of subsequent parametric tests.

### 3. Hotelling's T² Test
Tested whether the mean vectors of the ASD-positive and ASD-negative groups differ significantly across all variables simultaneously. The T² statistic was converted to an F-statistic and tested against the F-distribution.

### 4. MANOVA
Extended the group comparison using **Multivariate Analysis of Variance (MANOVA)** with both Wilks' Lambda and Pillai's Trace criteria, providing robustness against potential violations of normality.

### 5. Principal Component Analysis (PCA)
Performed PCA on the standardized data to reduce dimensionality and identify the directions of maximum variance. Selected components based on the **eigenvalue > 1 rule** and a scree plot, and examined component loadings to interpret each PC.

### 6. Factor Analysis
Applied **exploratory factor analysis** with maximum likelihood estimation and Varimax rotation (3 factors) to uncover the latent structure underlying the AQ-10 items. Loadings above 0.3 were retained for interpretation.

### 7. Linear Discriminant Analysis (LDA)
Built a **linear discriminant function** to classify individuals as ASD or non-ASD based on the 10 screening items and age. Evaluated classification performance using a confusion matrix and overall accuracy rate.

---

## Key Findings

- **Hotelling's T² and MANOVA** both yielded highly significant results, confirming that the multivariate mean profiles of ASD and non-ASD groups are substantially different.
- **Mardia's test** indicated departure from multivariate normality, which was expected given the binary nature of the screening items — this informed the interpretation of parametric results with appropriate caution.
- **PCA** revealed that the first few components capture a meaningful proportion of variance, with social communication items loading together on the leading component.
- **Factor Analysis** identified interpretable latent dimensions aligned with known ASD symptom domains (social interaction, communication, and attention to detail).
- **LDA** achieved high classification accuracy, demonstrating that the AQ-10 items carry strong discriminative signal for ASD screening.

---

## Repository Structure

```
multivariate-asd-screening/
├── asd_analysis.R          # Full R analysis script
├── autism_adult_ASD.csv    # Dataset
└── README.md
```

---

## How to Run

1. Clone the repository:
```bash
git clone https://github.com/Mohikhosravi/Multivariate-ASD-Screening.git
cd Multivariate-ASD-Screening
```

2. Open `asd_analysis.R` in RStudio and install required packages if needed:
```r
install.packages(c("ggplot2", "corrplot", "MASS", "psych",
                   "GGally", "mvnormtest", "biotools"))
```

3. Run the script. The dataset is loaded directly from the local CSV file.

---

## Tools & Environment

- **Language:** R
- **Key packages:** `ggplot2`, `corrplot`, `MASS`, `psych`, `GGally`
- **IDE:** RStudio

---

## Author

**Mohadese Khosravi**  
Statistics Student, Alzahra University  
[github.com/Mohikhosravi](https://github.com/Mohikhosravi)
