# ML Biomarker Discovery Script Analysis

## 1. Script Purpose
This script uses **machine learning** to identify which gut bacteria can serve as **biomarkers** to distinguish between your four study groups (normal, obesity, diabetic, and T2D+obesity patients). Think of it as finding bacterial "fingerprints" that are characteristic of each condition.

## 2. Key Functions

### Data Preparation
- Loads bacterial abundance data from previous taxonomic analysis
- Filters out rare bacteria (keeps only those present in >10% of samples)
- Transforms data for machine learning (log transformation, standardization)

### Feature Selection
- **Univariate Selection**: Picks bacteria that show the strongest statistical differences between groups
- **Recursive Feature Elimination**: Uses a smart algorithm to find the best combination of bacteria
- **Importance Ranking**: Identifies bacteria that contribute most to group differences
- **Differential Abundance Integration**: Uses results from previous statistical tests to prioritize important bacteria

### Machine Learning Classification
Trains four different AI models to classify samples:
- **Random Forest**: Uses multiple decision trees (like a committee of experts)
- **Gradient Boosting**: Builds models that learn from previous mistakes
- **Support Vector Machine (SVM)**: Finds the best boundaries between groups
- **Logistic Regression**: Uses statistical relationships for classification

### Model Evaluation
- **Cross-validation**: Tests models multiple times to ensure reliability
- **Performance metrics**: Measures accuracy, precision, and recall
- **ROC curves**: Shows how well models distinguish between groups

## 3. Input/Output

### Input Required:
- **Taxonomic abundance tables** (from previous analysis steps)
- **Sample metadata** with group classifications
- **Differential abundance results** (optional, for feature prioritization)

### Output Generated:
- **Trained ML models** saved for future use
- **Biomarker rankings** showing most important bacteria
- **Classification reports** with model performance statistics
- **Visualization plots** (accuracy comparisons, confusion matrices, ROC curves)
- **Dimensionality reduction plots** (PCA and t-SNE visualizations)

## 4. Methods Used

### Statistical Methods:
- **F-statistic** for univariate feature selection
- **Cross-validation** (5-fold stratified) for robust model evaluation
- **Principal Component Analysis (PCA)** for data visualization
- **t-SNE** for non-linear data visualization

### Machine Learning Algorithms:
- **Random Forest Classifier** with balanced class weights
- **Gradient Boosting Classifier** with learning rate optimization
- **Support Vector Machine** with RBF kernel
- **Logistic Regression** with L2 regularization

### Feature Selection Techniques:
- **SelectKBest** with F-statistic scoring
- **Recursive Feature Elimination with Cross-Validation (RFECV)**
- **Variance thresholding** to remove low-variance features
- **Ensemble feature selection** combining multiple methods

## 5. Results Generated

### Performance Metrics:
- **Accuracy scores** for each model (test and cross-validation)
- **Precision, recall, and F1-scores** for comprehensive evaluation
- **Confusion matrices** showing classification errors
- **ROC curves** with Area Under Curve (AUC) values

### Biomarker Identification:
- **Ranked lists** of bacterial taxa by importance score
- **Feature importance plots** showing top 20 biomarkers per taxonomic level
- **Combined biomarker panel** across all taxonomic levels
- **Model coefficients** indicating direction of association

### Visualizations:
- **Model comparison charts** showing relative performance
- **Biomarker importance plots** (horizontal bar charts)
- **PCA plots** showing sample clustering patterns
- **t-SNE plots** for complex pattern visualization
- **Combined biomarker visualization** color-coded by taxonomic level
