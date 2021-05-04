# OMG: Omics Mock Generator

Generates a mock dataset of _omics_ data (importable in EDD using the new import format): transcriptomics, proteomics, and metabolomics.

### Dependencies
- pip
- Python 3.8
- cobra
- pandas
- [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio?)


### Package Installation

Please follow the next steps after cloning the repository:
```bash
$ pip install pipenv
$ pipenv install
```
You may want to use the `--user` flag in `pip` and prepend `~/.local/bin` to your path.

### Jupyter Notebook

An example Jupyter Notebook, using the OMG module, is included here. To run,
build the Docker image, and run it to access a Jupyter Lab server:

```bash
DOCKER_BUILDKIT=1 docker build -t omg .
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$(pwd):/home/jovyan/work" omg
```

After the server launches, the terminal will have instructions on how to open
the Jupyter Lab server from your browser, looking something like the below:

```
    To access the server, open this file in a browser:
        file:///home/jovyan/.local/share/jupyter/runtime/jpserver-8-open.html
    Or copy and paste one of these URLs:
        http://<container_id>:8888/lab?token=<random token value>
        http://127.0.0.1:8888/lab?token=<random token value>
```

Copy the last URL and open in a web browser. The notebook `example.ipynb` will
be in the `work` folder. Run the example as-is, or make modifications to cells
in the notebook or uploaded input files.

### Contact
- For questions contact Somtirtha Roy at [somtirtharoy@lbl.gov](somtirtharoy@lbl.gov) or Jose M. Martí at [jmm@lbl.gov](jmm@lbl.gov).

### License
Omics Mock Generator Library (OMG) Copyright (c) 2021, The
Regents of the University of California, through Lawrence
Berkeley National Laboratory (subject to receipt of any required
approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative
works, and perform publicly and display publicly, and to permit others to do so.
