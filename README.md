[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ningz97/MLE-DSBM) 
MATLAB code for directed graph clustering algorithms proposed in our paper:\
[Maximum Likelihood Estimation on Stochastic Blockmodels for Directed Graph Clustering](http://arxiv.org/abs/2403.19516)
* You will need optimization toolboxes to implement the SDP:
  - **[recommended]** If you implement the Burer-Monteiro method, please install [Manopt](https://www.manopt.org/tutorial.html).
   - If you choose to solve SDP directly, please install [CVX](https://cvxr.com/cvx/).

- Contact: ning.zhang@stats.ox.ac.uk
[I plan to upload a Python version of the algorithms later!]


* Scripts included:
<pre>
   'Embedding.m': visualize MLE-driven embedding space
   'DSBM_visualize_adj': visualize adjacency relation before & after clustering on DSBM synthetic dataset
   'DSBM_comp.m': compare our algorithms (MLE-SC and MLE-SDP) with existing directed clustering methods
   'Iter_Convg_Rand.m': visualize how the iterative algorithm (Algorithm 4 in our paper) updates the DSBM parameters
   'cluster_algs/f_IT_MLE_sc.m': MLE-SC (Algorithm 4 + Algorithm 1) in our paper
   'cluster_algs/f_IT_MLE_bm.m': MLE-SDP (Algorithm 4 + Algorithm 3) in our paper (can replace Algorithm 3 with Algorithm 2 if you prefer using SDP solver)
</pre>

* Cite our paper:
```
  @misc{MLE_DSBM,
        title={Maximum Likelihood Estimation on Stochastic Blockmodels for Directed Graph Clustering}, 
        author={Mihai Cucuringu and Xiaowen Dong and Ning Zhang},
        year={2024},
        eprint={2403.19516},
        archivePrefix={arXiv},
        primaryClass={stat.ML}
        }
```
