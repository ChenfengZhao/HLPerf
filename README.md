# HLPerf

## Overview
The development of FPGA-based applications using HLS is fraught with performance pitfalls and large design space exploration times. Specifically, (1) RTL simulation is very time consuming for large input dataset, while tuning with a small dataset is less meaningful for irregular GNN kernels; (2) on-board measurement is also not attractive because the the hardware compilation usually takes several hours for each tuning iteration; (3) the notion of dataflow architectures further mystifies the optimization process, because it induces more challenges, such as task partitioning, FIFO depth tuning, and bottleneck identification. These issues are exacerbated when the application is complicated and its performance is dependent on the input dataset, as is often the case with graph neural network approaches to machine learning. Here, we introduce HLPerf, an simulation-based, approximately-cycle-accurate performance evaluation framework for dataflow architectures that both supports early exploration of the design space and shortens the performance evaluation cycle. HLPerf is a different method, effectively between the approaches of static estimation and cycle-accurate simulation, to investigate the impact of irregularity of data and algorithms on performance. We apply the methodology to [GNNHLS](https://github.com/ChenfengZhao/GNNHLS), an HLS-based graph neural network benchmark containing six commonly used graph neural network models and four datasets with distinct topologies and scales. The results show that HLPerf achieves dramatically better simulation speed than both RTL simulation and more recently developed cycle-accurate simulators at resonable accuracy cost. This acceleration positions HLPerf as a viable component in the design cycle.
More details are discussed in the [paper](https://dl.acm.org/doi/abs/10.1145/3655627).

## Requirements
- Recommended OS: Linux (Ubuntu >= 16.04). macOS and Windows also work.
- pypy3 >= 3.9 (recommended) or Python3 >= 3.8
- Simpy >= 4.0.1

## Installation Guide
HLPerf is based on Python syntax and [Simpy](https://simpy.readthedocs.io/en/4.0.1/simpy_intro/installation.html). We recommend pypy3 to execute the program because it is several times faster than Python3.

1. Install pypy3 (for Linux)
  ```
  sudo add-apt-repository ppa:pypy/ppa
  sudo apt update
  sudo apt install pypy3
  ```
2. Install Simpy
  ```
  pypy3 -m pip install simpy
  ```

We also provide a docker container with pypy3 and Simpy installed.

## Usage Example

1. The input graph data is a list of pointer pairs to the in-coming neighbors of target vertices. Uses can either directly use the datasets at ./data/ (e.g., ./data/OGBG-MOLTOX21/csr_indptr_trans.txt) or use the software stack of [GNNHLS](https://github.com/ChenfengZhao/GNNHLS) to convert any other graphs to the csr_indptr_trans.txt file.
   
2. Uncomment the command  Run the shell script run.sh which lists the commands for different models. In the shell script, $DATASET means the input graph dataset, $OUTPATH is the path to place output files.
  ```
  ./run.sh
  ```

Or run the following command for GCN with OGBG-MOLTOX21 in the terminal
```
mkdir -p ./result/GIN/OGBG-MOLTOX21
pypy3 gcn.py --ds_fn OGBG-MOLTOX21 | tee ./result/GCN/OGBG-MOLTOX21/result.log
```


3. After finishing the execution, an report (i.e., report.log) is generated. Enther the output directory and open result.log, where "Actual running time" represents the time it takes to run the program, "Simulation ends at" denotes the estimated cycle number, and "Calculated simulation time" represent the estimated performance in the unit of second based on the estimated cycle number.
  ```
  vi <output directory>/result.log
  ```


## License
[MIT_license]: https://spdx.org/licenses/MIT.html

The input data set is in the public domain. The source code of this project is released under the [MIT License][MIT_license]

## Citation
If you think it is helpful for your research, please cite the following [paper](https://dl.acm.org/doi/abs/10.1145/3655627):

<!-- Chenfeng Zhao, Zehao Dong, Yixin Chen, Xuan Zhang, and Roger D. Chamberlain. 2023. GNNHLS: Evaluating Graph Neural Network Inference via High-Level Synthesis. In Proc. of 41st IEEE International Conference on Computer Design (ICCD), November 6-8, 2023, Washington, DC, USA

[Arxiv's paper](https://arxiv.org/abs/2309.16022) -->

```
@article{zhao2024hlperf,
  title={HLPerf: Demystifying the Performance of HLS-based Graph Neural Networks with Dataflow Architectures},
  author={Zhao, Chenfeng and Faber, Clayton J and Chamberlain, Roger D and Zhang, Xuan},
  journal={ACM Transactions on Reconfigurable Technology and Systems},
  year={2024},
  publisher={ACM New York, NY}
}
```