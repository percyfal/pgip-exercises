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

(sec_intro)=

# Welcome!

This site contains a number of tutorials to develop your understanding of
[succinct tree sequences](https://tskit.dev/learn.html#what) and software programs,
such as [msprime](https://tskit.dev/msprime), that use them.


```{code-cell} python
import tskit
import msprime
a = 3+3
print(a)
```



```{math}
\pi = \sum_{ij} x_i x_j \pi_{ij}
```

```{code-cell} ipython3
:"tags": ["hide-input"]
import os
print(os.getcwd())
```

```{code-cell} ipython3
:"tags": ["remove-input"]
# This cell deliberately removed (not just hidden via a toggle) as it's not helpful
# for understanding tskit code (it's merely plotting code)
from IPython.display import SVG
```
