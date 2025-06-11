import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys


def load_pvalues(file_path):
    """
    Load p-values from the gene-based results file.
    """
    df = pd.read_csv(file_path, sep='\t')

    pvalues = pd.to_numeric(df['pvalue'], errors='coerce').dropna()
    # pvalues = df['pvalue'].astype(float)
    return pvalues


def calculate_inflation_factor(pvalues):
    """
    Calculate the Genomic Inflation Factor (λ) using the median chi-squared approach.
    """
    # Convert p-values to chi-squared statistics (1 degree of freedom)
    chisq_observed = np.quantile(-2 * np.log(pvalues), 0.5)
    print(chisq_observed)
    chisq_expected = np.quantile(np.random.chisquare(1, size=len(pvalues)), 0.5)
    print(chisq_expected)
    # Genomic Inflation Factor (λ)
    inflation_factor = chisq_observed / chisq_expected
    return inflation_factor


def make_qq_plot(pvalues, inflation_factor, output_file="qq_plot.png"):
    """
    Create QQ plot for the p-values and display the inflation factor (λ).
    """
    # Sort p-values and calculate expected p-values
    observed = -np.log10(np.sort(pvalues))
    expected = -np.log10(np.linspace(1 / len(pvalues), 1, len(pvalues)))

    # Create the QQ plot
    plt.figure(figsize=(8, 8))
    plt.scatter(expected, observed, c='blue', alpha=0.6, edgecolor='k')
    plt.plot([0, max(expected)], [0, max(expected)], color='red', linestyle='--')
    plt.title('QQ Plot of p-values')
    plt.xlabel('Expected -log10(p-value)')
    plt.ylabel('Observed -log10(p-value)')
    plt.grid(True)
    plt.tight_layout()

    # Display Inflation Factor (λ)
    plt.text(0.1, max(observed) * 0.9, f'λ = {inflation_factor:.3f}', fontsize=12, color='red')

    # Save and show the plot
    plt.savefig(output_file)
    plt.show()
    print(f"QQ plot saved to {output_file}")


def main(gene_file, output_file="qq_plot.png"):
    # Load p-values from the gene-based results
    pvalues = load_pvalues(gene_file)

    # Calculate Genomic Inflation Factor (λ)
    inflation_factor = calculate_inflation_factor(pvalues)
    print(f"Genomic Inflation Factor (λ): {inflation_factor:.3f}")

    # Make QQ plot
    make_qq_plot(pvalues, inflation_factor, output_file)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python qq_plot_with_inflation.py <gene_file> [output_file]")
        sys.exit(1)

    gene_file = sys.argv[1]  # Input gene-based results file
    output_file = sys.argv[2] if len(sys.argv) > 2 else "qq_plot.png"

    main(gene_file, output_file)
