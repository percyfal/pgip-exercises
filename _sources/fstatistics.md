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

(sec_fstatistics)=

# f-statistics #

# Exercises #

```{tabbed} Dsuite
```{code-block} shell
pwd
```

```{tabbed} sgkit
```{code-block} ipython3
import sgkit
```


```{code-block} shell
pwd
```

```{code-cell} ipython3
:"tags": ["remove-input"]
import subprocess as sp
_ = sp.check_call(["pwd"])
```


NOTE: Dsuite uses allele frequency estimates. Start example by
calculating stuff by hand.


Want also to relate the allele frequencies to branch lengths; cf
Peters
