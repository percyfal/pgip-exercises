---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: pgip
---


```{note}
About: exercises on variant calling. Input data consists of simulated and/or experimental data.
```

(sec_variantcalling)=

# Variant calling


## Best practice pipeline

Briefly describe best practice pipeline (e.g. GATK, bcftools,
freebayes), highlighting some of the problems with natural
populations:
- many workflows are human-centric and don't scale to non-model
  organisms (e.g. chromosome numbers in GATK)
  

## Variant calling exercise

```{tabbed} bcftools
```{code-block} shell
bcftools | head -2
```

```{tabbed} gatk
```{code-block} shell
gatk HaplotypeCaller
```



```{code-cell} ipython3
:"tags": ["remove-input", "hide-output"]
import subprocess
subprocess.run("bcftools | head -2", shell=True, check=False)
```

```{code-cell} ipython3
:"tags": ["remove-input", "hide-output"]
subprocess.run("gatk HaplotypeCaller", shell=True, check=False)
```

Optional: run bcftools / HaplotypeCaller / freebayes on data set or at
least show commands how files were generated

# References

https://link.springer.com/protocol/10.1007/978-1-4939-3578-9_11
