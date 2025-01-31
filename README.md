# wgcna-alternatives

 The repository for a report done in 2023 on WGCNA alternatives.

 See the final report [here](./final%20report/WGCNA%20Final%20Report.pdf)!

 Find the main workflow in the `workflow/` folder.

## How to run in Google Colab

1. Copy a Jupyter Notebook file: In Colab, on the top left corner, click `File` -> `Open notebook`. From there, click `GitHub` on the left, enter this repo's URL (https://github.com/benkoppe/wgcna-alternatives/) and open a Notebook in another tab.
2. On the top right corner, click the downwards-facing arrow next to the `Connect` button. Click `Connect to a local runtime` and follow the steps.
3. In the first code block, replace the lines:
```
import os
from pathlib import Path
dir = Path(os.path.abspath('../workflow')) # sets dir variable
```
with these:
```
from pathlib import Path
!git clone https://github.com/benkoppe/wgcna-alternatives.git
dir = Path("wgcna-alternatives/workflow")

!pip install rpy2
%load_ext rpy2.ipython
```
Then run the first code block with the attached forwards-facing arrow.
4. Next, click the `+ Code` button on the top. Inside the block, paste in this code:
```
%%R
install.packages("BiocManager")
BiocManager::install()
BiocManager::install("GWENA")
BiocManager::install("CEMiTool")
BiocManager::install("limma")
```

That's it! Your notebook should now be fully working within Colab. Be sure to delete the files saved to your machine in this process when you're done. Unfortunately, files must be saved to your machine due to problems loading Bioconductor packages in Colab.
