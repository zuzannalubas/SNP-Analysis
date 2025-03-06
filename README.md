# 🧬 SNP Analysis  

## 📌 Description  
This project analyzes **Single Nucleotide Polymorphisms (SNPs)** to identify significant genetic variations associated with specific phenotypes.  
We perform various statistical and visualization techniques to uncover potential genetic markers.  

## 🎯 Objectives  
- **Data Preprocessing**: Clean and filter SNP data (remove low-frequency variants).  
- **Minor Allele Frequency (MAF) Analysis**: Identify the distribution of rare alleles.  
- **Statistical Tests**: Use Chi-Square tests to check SNP-phenotype associations.  
- **Manhattan Plot**: Visualize the significance of SNPs across the genome.  
- **PCA Analysis**: Explore clustering of SNPs based on genetic variation.  

## 📂 Dataset  
The dataset consists of SNPs with the following attributes:  
- `SNP_ID`: Unique identifier for each SNP  
- `Chromosome`: Chromosome number where SNP is located  
- `Position`: Genomic position of the SNP  
- `Reference_Allele`: The standard nucleotide at this position  
- `Alternate_Allele`: The variant nucleotide  
- `MAF`: Minor Allele Frequency (frequency of the least common allele)  
- `p-value`: Statistical significance of association with phenotype  
- `Phenotype`: Binary trait (0 = control, 1 = affected)  

## 🛠 Technologies  
- **Python** 🐍  
- **Pandas** 📊 (Data processing)  
- **NumPy** 🔢 (Numerical analysis)  
- **Matplotlib & Seaborn** 📈 (Data visualization)  
- **SciPy** 🏛️ (Statistical tests)  
- **scikit-learn** 🤖 (PCA analysis)  

## 📊 Results & Visualizations  
- **MAF Histogram**: Shows allele frequency distribution.  
- **Chi-Square Test**: Identifies statistically significant SNPs.  
- **Manhattan Plot**: Highlights significant SNPs with p-values.  
- **PCA Scatter Plot**: Reveals clustering patterns of SNPs.  

## 🚀 Usage  
1️⃣ **Install required libraries**  
```bash
pip install pandas numpy matplotlib seaborn scipy scikit-learn
