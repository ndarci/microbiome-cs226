# microbiome-cs226
Final project for Bioinformatics 226 at UCLA. Analyzing the effect of gut microbiome species composition on LDL and HDL cholesterol. 

John Randazzo, Maddie Murphy, Niko Darci-Maher

# Merge data and select baseline samples
```{bash}
Rscript cleanup.R
```

# Run dimensionality reduction
```{bash}
Rscript PCA.R
```

# Run penalized regression analysis
```{bash}
Rscript model_comparison.R
```

# Compute dysbiosis scores and regress on them
```{bash}
Rscript dysbiosis_analysis.R
```
