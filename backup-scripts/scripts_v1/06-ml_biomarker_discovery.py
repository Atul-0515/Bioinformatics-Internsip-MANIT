#!/usr/bin/env python3

"""
Machine Learning Biomarker Discovery Script for 16S rRNA Microbiome Data
Identifies potential microbial biomarkers using various ML approaches
Integrates taxonomic abundance data with differential abundance results
Author: Bioinformatics Pipeline
Date: 2025
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.feature_selection import SelectKBest, f_classif, RFECV
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import warnings

warnings.filterwarnings("ignore")

# Set style for plots
plt.style.use("default")
sns.set_palette("husl")


class MicrobiomeBiomarkerDiscovery:
    def __init__(self, results_dir="results"):
        self.results_dir = results_dir
        self.ml_dir = os.path.join(results_dir, "ml_biomarker_discovery")
        self.create_directories()

        # Initialize data containers
        self.abundance_data = {}
        self.metadata = None
        self.differential_results = {}
        self.ml_results = {}

    def create_directories(self):
        """Create necessary output directories"""
        dirs = [
            self.ml_dir,
            os.path.join(self.ml_dir, "models"),
            os.path.join(self.ml_dir, "classification_reports"),
            os.path.join(self.ml_dir, "dimensionality_reduction"),
            os.path.join(self.ml_dir, "biomarker_panels"),
        ]

        for directory in dirs:
            os.makedirs(directory, exist_ok=True)

    def load_taxonomic_data(self):
        """Load taxonomic abundance data and metadata"""
        print("Loading taxonomic abundance data...")

        # Load metadata
        metadata_file = os.path.join(self.results_dir, "taxonomy", "sample_categories.csv")
        if not os.path.exists(metadata_file):
            raise FileNotFoundError("Sample categories file not found. Run taxonomic profiling first.")

        self.metadata = pd.read_csv(metadata_file)
        print(f"Loaded metadata for {len(self.metadata)} samples")
        print(f"Categories: {', '.join(self.metadata['Category'].unique())}")

        # Load abundance data for different taxonomic levels
        taxonomy_dir = os.path.join(self.results_dir, "taxonomy", "merged_tables")
        levels = ["species", "genus", "family", "order", "class", "phylum"]

        for level in levels:
            abundance_file = os.path.join(taxonomy_dir, f"{level}_relative_abundance.csv")
            if os.path.exists(abundance_file):
                df = pd.read_csv(abundance_file, index_col=0)
                # Filter low abundance taxa (present in >10% samples with >0.1% abundance)
                min_samples = int(len(df.columns) * 0.1)
                abundant_taxa = (df > 0.001).sum(axis=1) >= min_samples
                self.abundance_data[level] = df[abundant_taxa]
                print(f"Loaded {level} data: {len(df[abundant_taxa])} taxa, {len(df.columns)} samples")

        if not self.abundance_data:
            raise ValueError("No abundance data found after filtering")

    def load_differential_results(self):
        """Load differential abundance results to prioritize features"""
        print("Loading differential abundance results...")

        diff_dir = os.path.join(self.results_dir, "differential_abundance", "tables")
        if not os.path.exists(diff_dir):
            print("Warning: Differential abundance results not found. Proceeding without feature prioritization.")
            return

        for level in self.abundance_data.keys():
            self.differential_results[level] = {}

            # Load DESeq2 results
            deseq_files = [f for f in os.listdir(diff_dir) if f.startswith(f"{level}_deseq2_")]
            for file in deseq_files:
                comparison = file.replace(f"{level}_deseq2_", "").replace(".csv", "")
                df = pd.read_csv(os.path.join(diff_dir, file))
                # Filter significant results
                significant = df[(df["padj"] < 0.05) & (abs(df["log2FoldChange"]) > 1)]
                if not significant.empty:
                    self.differential_results[level][f"deseq2_{comparison}"] = significant

            # Load ALDEx2 results
            aldex_files = [f for f in os.listdir(diff_dir) if f.startswith(f"{level}_aldex2_")]
            for file in aldex_files:
                comparison = file.replace(f"{level}_aldex2_", "").replace(".csv", "")
                df = pd.read_csv(os.path.join(diff_dir, file))
                significant = df[(df["p_value_bh"] < 0.05) & (abs(df["effect_size"]) > 1)]
                if not significant.empty:
                    self.differential_results[level][f"aldex2_{comparison}"] = significant

        print(f"Loaded differential results for {len(self.differential_results)} taxonomic levels")

    def prepare_ml_data(self, level, target_column="Category"):
        """Prepare data for machine learning"""
        print(f"Preparing ML data for {level} level...")

        # Get abundance data
        X = self.abundance_data[level].T  # Transpose so samples are rows

        # Align with metadata
        common_samples = list(set(X.index) & set(self.metadata["SampleID"]))
        X = X.loc[common_samples]
        y_metadata = self.metadata.set_index("SampleID").loc[common_samples, target_column]

        # Encode labels
        le = LabelEncoder()
        y = le.fit_transform(y_metadata)

        # Log transform abundance data
        X_log = np.log10(X + 1e-6)

        # Remove features with zero variance
        var_filter = X_log.var() > 0
        X_filtered = X_log.loc[:, var_filter]

        print(f"Data shape: {X_filtered.shape}")
        print(f"Classes: {dict(zip(le.classes_, range(len(le.classes_))))}")

        return X_filtered, y, le.classes_

    def feature_selection(self, X, y, level, method="all", k=50):
        """Perform feature selection with improved methods"""
        print(f"Performing feature selection for {level} level...")

        selected_features = {}

        # 1. Univariate feature selection
        if method in ["all", "univariate"]:
            selector = SelectKBest(score_func=f_classif, k=min(k, X.shape[1]))
            # REMOVE THIS LINE: X_selected = selector.fit_transform(X, y)  # UNUSED VARIABLE
            selector.fit(X, y)  # Just fit, don't assign unused variable
            selected_idx = selector.get_support(indices=True)
            selected_features["univariate"] = X.columns[selected_idx].tolist()

        # 2. Recursive Feature Elimination with Cross-Validation (FIXED)
        if method in ["all", "rfecv"]:
            rf = RandomForestClassifier(n_estimators=100, random_state=42)
            # Use RFECV instead of RFE for better feature selection
            cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)  # Now using StratifiedKFold
            selector = RFECV(estimator=rf, cv=cv, scoring="accuracy", min_features_to_select=5)
            selector.fit(X, y)
            selected_features["rfecv"] = X.columns[selector.support_].tolist()
            print(f"    RFECV selected {selector.n_features_} features (optimal)")

        # 3. Feature importance from Random Forest
        if method in ["all", "importance"]:
            rf = RandomForestClassifier(n_estimators=100, random_state=42)
            rf.fit(X, y)
            importances = pd.Series(rf.feature_importances_, index=X.columns)
            top_features = importances.nlargest(min(k, len(importances))).index.tolist()
            selected_features["importance"] = top_features

        # 4. Features from differential abundance (IMPROVED)
        if level in self.differential_results and method in ["all", "differential"]:
            diff_features = set()
            for comparison, results in self.differential_results[level].items():
                if "taxon" in results.columns:
                    # Prioritize by effect size/fold change
                    if "log2FoldChange" in results.columns:
                        # DESeq2 results - sort by absolute log2FoldChange
                        sorted_results = results.reindex(results["log2FoldChange"].abs().sort_values(ascending=False).index)
                    elif "effect_size" in results.columns:
                        # ALDEx2 results - sort by absolute effect size
                        sorted_results = results.reindex(results["effect_size"].abs().sort_values(ascending=False).index)
                    else:
                        sorted_results = results

                    # Take top differential features
                    valid_features = [f for f in sorted_results["taxon"][:30] if f in X.columns]
                    diff_features.update(valid_features)

            if diff_features:
                selected_features["differential"] = list(diff_features)
                print(f"    Added {len(diff_features)} differentially abundant features")

        # 5. Variance-based selection (NEW)
        if method in ["all", "variance"]:
            # Remove low-variance features
            from sklearn.feature_selection import VarianceThreshold

            variance_threshold = VarianceThreshold(threshold=0.01)  # Remove features with very low variance
            variance_threshold.fit(X)
            high_var_features = X.columns[variance_threshold.get_support()].tolist()
            selected_features["variance"] = high_var_features[:k]

        # Combine all selected features with scoring
        feature_scores = {}
        for method_name, feature_list in selected_features.items():
            for feature in feature_list:
                if feature not in feature_scores:
                    feature_scores[feature] = 0
                feature_scores[feature] += 1

        # Select features that appear in multiple methods (ensemble approach)
        if len(selected_features) > 1:
            # Prioritize features selected by multiple methods
            sorted_features = sorted(feature_scores.items(), key=lambda x: x[1], reverse=True)
            combined_features = [f[0] for f in sorted_features[: min(100, len(sorted_features))]]
        else:
            # If only one method used, take all features
            all_selected = set()
            for feature_list in selected_features.values():
                all_selected.update(feature_list)
            combined_features = list(all_selected)[: min(100, len(all_selected))]

        selected_features["combined"] = combined_features

        print(f"Selected {len(combined_features)} features using combined methods")
        print(f"Feature selection methods used: {list(selected_features.keys())}")

        return selected_features

    def train_classifiers(self, X, y, class_names, level):
        """Train multiple classifiers with improved cross-validation"""
        print(f"Training classifiers for {level} level...")

        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Split data
        X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.3, random_state=42, stratify=y)

        # Define classifiers with better parameters
        classifiers = {
            "Random Forest": RandomForestClassifier(
                n_estimators=100,
                random_state=42,
                class_weight="balanced",  # Handle class imbalance
            ),
            "Gradient Boosting": GradientBoostingClassifier(n_estimators=100, random_state=42, learning_rate=0.1),
            "SVM": SVC(
                probability=True,
                random_state=42,
                class_weight="balanced",  # Handle class imbalance
                kernel="rbf",
            ),
            "Logistic Regression": LogisticRegression(
                random_state=42,
                max_iter=1000,
                class_weight="balanced",  # Handle class imbalance
            ),
        }

        results = {}

        # Use StratifiedKFold for better cross-validation (NOW USING THE IMPORTED CLASS)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

        for name, clf in classifiers.items():
            print(f"  Training {name}...")

            # Train classifier
            clf.fit(X_train, y_train)

            # Predictions
            y_pred = clf.predict(X_test)
            y_pred_proba = clf.predict_proba(X_test)

            # Cross-validation scores with multiple metrics
            cv_accuracy = cross_val_score(clf, X_scaled, y, cv=cv, scoring="accuracy")
            cv_precision = cross_val_score(clf, X_scaled, y, cv=cv, scoring="precision_macro")
            cv_recall = cross_val_score(clf, X_scaled, y, cv=cv, scoring="recall_macro")
            cv_f1 = cross_val_score(clf, X_scaled, y, cv=cv, scoring="f1_macro")

            # Store results
            results[name] = {
                "classifier": clf,
                "scaler": scaler,
                "y_test": y_test,
                "y_pred": y_pred,
                "y_pred_proba": y_pred_proba,
                "cv_accuracy": cv_accuracy,
                "cv_precision": cv_precision,
                "cv_recall": cv_recall,
                "cv_f1": cv_f1,
                "test_accuracy": (y_pred == y_test).mean(),
                "cv_accuracy_mean": cv_accuracy.mean(),
                "cv_accuracy_std": cv_accuracy.std(),
                "cv_precision_mean": cv_precision.mean(),
                "cv_recall_mean": cv_recall.mean(),
                "cv_f1_mean": cv_f1.mean(),
            }

            print(f"    Test Accuracy: {results[name]['test_accuracy']:.3f}")
            print(f"    CV Accuracy: {results[name]['cv_accuracy_mean']:.3f} ± {results[name]['cv_accuracy_std']:.3f}")
            print(f"    CV F1-Score: {results[name]['cv_f1_mean']:.3f}")

        return results

    def evaluate_models(self, results, class_names, level):
        """Evaluate and visualize model performance"""
        print(f"Evaluating models for {level} level...")

        # Create comparison plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f"Model Comparison - {level.title()} Level", fontsize=16)

        # 1. Accuracy comparison - FIXED
        models = list(results.keys())
        test_acc = [results[m]["test_accuracy"] for m in models]
        cv_acc_mean = [results[m]["cv_accuracy_mean"] for m in models]  # FIXED: use cv_accuracy_mean
        cv_std = [results[m]["cv_accuracy_std"] for m in models]  # FIXED: use cv_accuracy_std

        x_pos = np.arange(len(models))
        axes[0, 0].bar(x_pos - 0.2, test_acc, 0.4, label="Test Accuracy", alpha=0.8)
        axes[0, 0].errorbar(
            x_pos + 0.2,
            cv_acc_mean,
            yerr=cv_std,
            fmt="o",
            label="CV Accuracy",
            capsize=5,
        )
        axes[0, 0].set_xlabel("Models")
        axes[0, 0].set_ylabel("Accuracy")
        axes[0, 0].set_title("Model Accuracy Comparison")
        axes[0, 0].set_xticks(x_pos)
        axes[0, 0].set_xticklabels(models, rotation=45)
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)

        # 2. Confusion matrix for best model
        best_model = max(results.keys(), key=lambda x: results[x]["cv_accuracy_mean"])
        cm = confusion_matrix(results[best_model]["y_test"], results[best_model]["y_pred"])

        axes[0, 1].imshow(cm, interpolation="nearest", cmap=plt.cm.Blues)
        axes[0, 1].set_title(f"Confusion Matrix - {best_model}")
        axes[0, 1].set_xlabel("Predicted Label")
        axes[0, 1].set_ylabel("True Label")

        # Add text annotations to confusion matrix
        thresh = cm.max() / 2.0
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                axes[0, 1].text(
                    j,
                    i,
                    format(cm[i, j], "d"),
                    ha="center",
                    va="center",
                    color="white" if cm[i, j] > thresh else "black",
                )

        axes[0, 1].set_xticks(range(len(class_names)))
        axes[0, 1].set_yticks(range(len(class_names)))
        axes[0, 1].set_xticklabels(class_names, rotation=45)
        axes[0, 1].set_yticklabels(class_names)

        # 3. ROC curves (for binary or multiclass)

        # Replace the ROC curve section in evaluate_models method:

        # 3. ROC curves (FIXED for both binary and multiclass)
        # 3. ROC curves (FIXED for both binary and multiclass)
        if len(class_names) == 2:
            # Binary classification ROC
            for name, result in results.items():
                fpr, tpr, _ = roc_curve(result["y_test"], result["y_pred_proba"][:, 1])
                auc = roc_auc_score(result["y_test"], result["y_pred_proba"][:, 1])
                axes[1, 0].plot(fpr, tpr, label=f"{name} (AUC = {auc:.3f})")

            axes[1, 0].plot([0, 1], [0, 1], "k--", alpha=0.5)
            axes[1, 0].set_xlabel("False Positive Rate")
            axes[1, 0].set_ylabel("True Positive Rate")
            axes[1, 0].set_title("ROC Curves")
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)
        else:
            # Multiclass ROC - One vs Rest approach
            from sklearn.preprocessing import label_binarize
            # REMOVE THESE CONFLICTING IMPORTS - they're already imported at top:
            # from sklearn.metrics import roc_curve, auc

            # Select best model for multiclass ROC
            best_model = max(results.keys(), key=lambda x: results[x]["cv_accuracy_mean"])
            result = results[best_model]

            # Binarize the output
            y_test_bin = label_binarize(result["y_test"], classes=range(len(class_names)))
            y_score = result["y_pred_proba"]

            # Compute ROC curve and ROC area for each class
            fpr = dict()
            tpr = dict()
            roc_auc = dict()

            for i in range(len(class_names)):
                fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], y_score[:, i])
                # Use sklearn.metrics.auc instead of the local auc import
                from sklearn.metrics import auc as sklearn_auc

                roc_auc[i] = sklearn_auc(fpr[i], tpr[i])
                axes[1, 0].plot(fpr[i], tpr[i], label=f"{class_names[i]} (AUC = {roc_auc[i]:.3f})")

            axes[1, 0].plot([0, 1], [0, 1], "k--", alpha=0.5)
            axes[1, 0].set_xlabel("False Positive Rate")
            axes[1, 0].set_ylabel("True Positive Rate")
            axes[1, 0].set_title(f"Multiclass ROC Curves - {best_model}")
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)

        # 4. Cross-validation scores distribution
        cv_data = [results[m]["cv_accuracy"] for m in models]  # This is the array of CV scores
        axes[1, 1].boxplot(cv_data, labels=models)
        axes[1, 1].set_title("Cross-Validation Score Distribution")
        axes[1, 1].set_xlabel("Models")
        axes[1, 1].set_ylabel("Accuracy")
        axes[1, 1].tick_params(axis="x", rotation=45)
        axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(
            os.path.join(self.ml_dir, f"{level}_model_comparison.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

        # Save detailed classification report
        report_file = os.path.join(self.ml_dir, "classification_reports", f"{level}_classification_report.txt")
        with open(report_file, "w") as f:
            f.write(f"Classification Report - {level.title()} Level\n")
            f.write("=" * 50 + "\n\n")

            for name, result in results.items():
                f.write(f"\n{name}:\n")
                f.write("-" * 30 + "\n")
                f.write(f"Test Accuracy: {result['test_accuracy']:.4f}\n")
                f.write(
                    f"CV Accuracy: {result['cv_accuracy_mean']:.4f} ± {result['cv_accuracy_std']:.4f}\n"  # FIXED
                )
                f.write(f"CV Precision: {result['cv_precision_mean']:.4f}\n")  # ADD THESE
                f.write(f"CV Recall: {result['cv_recall_mean']:.4f}\n")  # ADD THESE
                f.write(f"CV F1-Score: {result['cv_f1_mean']:.4f}\n")  # ADD THESE
                f.write("\nDetailed Classification Report:\n")
                f.write(classification_report(result["y_test"], result["y_pred"], target_names=class_names))
                f.write("\n" + "=" * 50 + "\n")

    def extract_biomarkers(self, results, feature_names, class_names, level):
        """Extract and rank potential biomarkers"""
        print(f"Extracting biomarkers for {level} level...")

        biomarkers = []

        # Ensure feature_names is a list
        if hasattr(feature_names, "tolist"):
            feature_names = feature_names.tolist()
        elif not isinstance(feature_names, list):
            feature_names = list(feature_names)

        # Get feature importances from Random Forest and Gradient Boosting
        for model_name in ["Random Forest", "Gradient Boosting"]:
            if model_name in results:
                clf = results[model_name]["classifier"]
                if hasattr(clf, "feature_importances_"):
                    importances = clf.feature_importances_
                    # Ensure we don't exceed the feature_names length
                    n_features = min(len(importances), len(feature_names))
                    for i in range(n_features):
                        biomarkers.append(
                            {
                                "feature": feature_names[i],
                                "importance": importances[i],
                                "model": model_name,
                                "level": level,
                            }
                        )

        # Get coefficients from Logistic Regression
        if "Logistic Regression" in results:
            clf = results["Logistic Regression"]["classifier"]
            if hasattr(clf, "coef_"):
                coef = clf.coef_
                if len(class_names) == 2:
                    # Binary classification
                    n_features = min(len(coef[0]), len(feature_names))
                    for i in range(n_features):
                        biomarkers.append(
                            {
                                "feature": feature_names[i],
                                "importance": abs(coef[0][i]),
                                "coefficient": coef[0][i],
                                "model": "Logistic Regression",
                                "level": level,
                            }
                        )
                else:
                    # Multiclass - use mean absolute coefficient
                    mean_coef = np.mean(np.abs(coef), axis=0)
                    n_features = min(len(mean_coef), len(feature_names))
                    for i in range(n_features):
                        biomarkers.append(
                            {
                                "feature": feature_names[i],
                                "importance": mean_coef[i],
                                "model": "Logistic Regression",
                                "level": level,
                            }
                        )

        # Convert to DataFrame and aggregate
        biomarker_df = pd.DataFrame(biomarkers)
        if not biomarker_df.empty:
            # Aggregate importances across models
            agg_biomarkers = biomarker_df.groupby("feature").agg({"importance": ["mean", "std", "count"], "level": "first"}).round(4)

            agg_biomarkers.columns = [
                "mean_importance",
                "std_importance",
                "n_models",
                "level",
            ]
            agg_biomarkers = agg_biomarkers.sort_values("mean_importance", ascending=False)

            # Save biomarker rankings
            output_file = os.path.join(self.ml_dir, "biomarker_panels", f"{level}_biomarker_ranking.csv")
            agg_biomarkers.to_csv(output_file)

            # Create biomarker importance plot
            top_biomarkers = agg_biomarkers.head(20)

            plt.figure(figsize=(12, 8))
            y_pos = np.arange(len(top_biomarkers))
            plt.barh(
                y_pos,
                top_biomarkers["mean_importance"],
                xerr=top_biomarkers["std_importance"],
                alpha=0.8,
            )
            plt.yticks(y_pos, top_biomarkers.index)
            plt.xlabel("Mean Feature Importance")
            plt.title(f"Top 20 Biomarkers - {level.title()} Level")
            plt.gca().invert_yaxis()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(
                os.path.join(self.ml_dir, "biomarker_panels", f"{level}_top_biomarkers.png"),
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()

            return agg_biomarkers

        return pd.DataFrame()

    def dimensionality_reduction(self, X, y, class_names, level):
        """Perform dimensionality reduction and visualization"""
        print(f"Performing dimensionality reduction for {level} level...")

        # Standardize data
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        fig, axes = plt.subplots(1, 2, figsize=(15, 6))

        # PCA
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(X_scaled)

        for i, class_name in enumerate(class_names):
            mask = y == i
            axes[0].scatter(X_pca[mask, 0], X_pca[mask, 1], label=class_name, alpha=0.7, s=50)

        axes[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)")
        axes[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)")
        axes[0].set_title("PCA Visualization")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        # t-SNE
        if len(X) > 50:  # Only run t-SNE if we have enough samples
            tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(X) // 4))
            X_tsne = tsne.fit_transform(X_scaled)

            for i, class_name in enumerate(class_names):
                mask = y == i
                axes[1].scatter(X_tsne[mask, 0], X_tsne[mask, 1], label=class_name, alpha=0.7, s=50)

            axes[1].set_xlabel("t-SNE 1")
            axes[1].set_ylabel("t-SNE 2")
            axes[1].set_title("t-SNE Visualization")
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)
        else:
            axes[1].text(
                0.5,
                0.5,
                "Not enough samples\nfor t-SNE\n(minimum 50 required)",
                ha="center",
                va="center",
                transform=axes[1].transAxes,
            )
            axes[1].set_title("t-SNE Visualization")

        plt.tight_layout()
        plt.savefig(
            os.path.join(self.ml_dir, "dimensionality_reduction", f"{level}_dim_reduction.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

        return X_pca, pca.explained_variance_ratio_

    def save_models(self, results, level):
        """Save trained models for future use"""
        import joblib

        models_dir = os.path.join(self.ml_dir, "models")

        for name, result in results.items():
            model_file = os.path.join(models_dir, f"{level}_{name.lower().replace(' ', '_')}_model.pkl")
            scaler_file = os.path.join(models_dir, f"{level}_{name.lower().replace(' ', '_')}_scaler.pkl")

            # Save model and scaler
            joblib.dump(result["classifier"], model_file)
            joblib.dump(result["scaler"], scaler_file)

        print(f"Models saved for {level} level")

    def run_analysis(self):
        """Run complete biomarker discovery analysis"""
        print("Starting ML Biomarker Discovery Analysis...")
        print("=" * 50)

        # Load data
        self.load_taxonomic_data()
        self.load_differential_results()

        # Analysis summary
        summary_results = {}
        all_biomarkers = []

        # Process each taxonomic level
        for level in self.abundance_data.keys():
            print(f"\n--- Processing {level.title()} Level ---")

            try:
                # Prepare data
                X, y, class_names = self.prepare_ml_data(level)

                # Feature selection
                selected_features = self.feature_selection(X, y, level)

                # Use combined selected features
                if "combined" in selected_features and selected_features["combined"]:
                    X_selected = X[selected_features["combined"]]
                    feature_names = selected_features["combined"]
                else:
                    X_selected = X
                    feature_names = X.columns.tolist()

                print(f"Using {len(feature_names)} features for ML")

                # Train classifiers
                results = self.train_classifiers(X_selected, y, class_names, level)

                # Evaluate models
                self.evaluate_models(results, class_names, level)
                self.save_models(results, level)

                # Extract biomarkers
                biomarkers = self.extract_biomarkers(results, feature_names, class_names, level)
                if not biomarkers.empty:
                    all_biomarkers.append(biomarkers.head(10))  # Top 10 per level

                # Dimensionality reduction
                X_pca, variance_ratio = self.dimensionality_reduction(X_selected, y, class_names, level)

                # Store summary
                best_model = max(
                    results.keys(),
                    key=lambda x: results[x]["cv_accuracy_mean"],  # FIXED
                )
                summary_results[level] = {
                    "n_samples": len(X),
                    "n_features": len(feature_names),
                    "n_classes": len(class_names),
                    "best_model": best_model,
                    "best_accuracy": results[best_model]["cv_accuracy_mean"],  # FIXED
                    "pca_variance": variance_ratio[:2].sum(),
                }

            except Exception as e:
                print(f"Error processing {level} level: {str(e)}")
                continue

        # Create overall summary
        self.create_summary_report(summary_results, all_biomarkers)

        print("\n" + "=" * 50)
        print("ML Biomarker Discovery Analysis Complete!")
        print(f"Results saved in: {self.ml_dir}")
        print("=" * 50)

    def create_summary_report(self, summary_results, all_biomarkers):
        """Create a comprehensive summary report"""
        print("Creating summary report...")

        # Summary statistics
        summary_file = os.path.join(self.ml_dir, "analysis_summary.txt")
        with open(summary_file, "w") as f:
            f.write("ML Biomarker Discovery Analysis Summary\n")
            f.write("=" * 50 + "\n\n")

            f.write("Dataset Overview:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Sample Categories: {', '.join(self.metadata['Category'].unique())}\n")
            f.write(f"Total Samples: {len(self.metadata)}\n")
            category_counts = self.metadata["Category"].value_counts()
            for category, count in category_counts.items():
                f.write(f"  {category}: {count} samples\n")
            f.write("\n")

            f.write("Taxonomic Level Results:\n")
            f.write("-" * 25 + "\n")
            for level, results in summary_results.items():
                f.write(f"\n{level.title()} Level:\n")
                f.write(f"  Samples: {results['n_samples']}\n")
                f.write(f"  Features: {results['n_features']}\n")
                f.write(f"  Classes: {results['n_classes']}\n")
                f.write(f"  Best Model: {results['best_model']}\n")
                f.write(f"  Best Accuracy: {results['best_accuracy']:.4f}\n")
                f.write(f"  PCA Variance (PC1+PC2): {results['pca_variance']:.4f}\n")

            f.write("\n" + "=" * 50 + "\n")
            f.write("Key Findings:\n")
            f.write("-" * 15 + "\n")

            if summary_results:
                best_level = max(
                    summary_results.keys(),
                    key=lambda x: summary_results[x]["best_accuracy"],
                )
                f.write(f"• Best performing taxonomic level: {best_level.title()}\n")
                f.write(f"• Highest accuracy achieved: {summary_results[best_level]['best_accuracy']:.4f}\n")
                f.write(f"• Best model type: {summary_results[best_level]['best_model']}\n")

            f.write("\n• Check individual level results for detailed biomarker rankings\n")
            f.write("• Heatmaps and visualizations available in respective folders\n")

        # Combined biomarker ranking - FIXED
        if all_biomarkers:
            # Reset index to make 'feature' a regular column instead of index
            combined_biomarkers_list = []
            for biomarker_df in all_biomarkers:
                if not biomarker_df.empty:
                    # Reset index to convert the grouped feature names back to a column
                    df_reset = biomarker_df.reset_index()
                    combined_biomarkers_list.append(df_reset)

            if combined_biomarkers_list:
                combined_biomarkers = pd.concat(combined_biomarkers_list, ignore_index=True)
                combined_biomarkers = combined_biomarkers.sort_values("mean_importance", ascending=False)

                combined_file = os.path.join(self.ml_dir, "biomarker_panels", "combined_top_biomarkers.csv")
                combined_biomarkers.to_csv(combined_file, index=False)

                # Create combined biomarker visualization
                plt.figure(figsize=(14, 10))
                top_combined = combined_biomarkers.head(30)

                # Color by taxonomic level
                level_colors = {
                    "species": "red",
                    "genus": "blue",
                    "family": "green",
                    "order": "orange",
                    "class": "purple",
                    "phylum": "brown",
                }
                colors = [level_colors.get(level, "gray") for level in top_combined["level"]]

                y_pos = np.arange(len(top_combined))
                plt.barh(y_pos, top_combined["mean_importance"], color=colors, alpha=0.7)
                plt.yticks(y_pos, top_combined["feature"], fontsize=8)
                plt.xlabel("Mean Feature Importance")
                plt.title("Top 30 Biomarkers Across All Taxonomic Levels")
                plt.gca().invert_yaxis()

                # Create legend
                legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.7, label=level.title()) for level, color in level_colors.items() if level in top_combined["level"].values]
                plt.legend(handles=legend_elements, loc="lower right")

                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.savefig(
                    os.path.join(self.ml_dir, "biomarker_panels", "combined_top_biomarkers.png"),
                    dpi=300,
                    bbox_inches="tight",
                )
                plt.close()


def main():
    """Main function to run the analysis"""
    import argparse

    parser = argparse.ArgumentParser(description="ML Biomarker Discovery for Microbiome Data")
    parser.add_argument(
        "--results_dir",
        type=str,
        default="results",
        help="Path to results directory (default: results)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=None,
        help="Custom output directory (default: results/ml_biomarker_discovery)",
    )

    args = parser.parse_args()

    try:
        # Initialize analyzer
        if args.output_dir:
            analyzer = MicrobiomeBiomarkerDiscovery(args.output_dir)
        else:
            analyzer = MicrobiomeBiomarkerDiscovery(args.results_dir)

        # Run analysis
        analyzer.run_analysis()

    except FileNotFoundError as e:
        print(f"Error: Required input files not found - {e}")
        print("Please ensure taxonomic profiling has been completed first.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
