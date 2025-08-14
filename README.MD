Segmentation Difference is algorithm for comparison of two segmentations to determine how close they are.


# Windows

For Windows, use Visual Studio with SD3.sln and build the solution.
You can build 'Reindex' project to reindex second segmentation to match segment numbers in best effort approach.


# Python tool

In a Unix shell:

```shell
pixi install
pixi shell
source compile.cmd
```

And then in a Python3 interpreter:

```python
from SD3 import sd3
metric_value = sd3(numpy_2d_array, numpy_2d_array)
```
