
MATLAB code for directed graph clustering algorithms proposed in our paper [Maximum Likelihood Estimation on Stochastic Blockmodels
for Directed Graph Clustering](http://arxiv.org/abs/2403.19516)

- If you implement the SDP for directed clustering, please install [CVX](https://cvxr.com/cvx/).\
- If you implement the Burer-Monteiro method for directed clustering, please install [Manopt](https://www.manopt.org/tutorial.html).

- Contact: ning.zhang@stats.ox.ac.uk

* scripts included:
```
   Embedding.m: visualize MLE-driven embedding space\
   DSBM_visualize_adj: visualize adjacency relation before & after clustering on DSBM synthetic dataset\
   DSBM_comp: compare our algorithms (MLE-SC and MLE-SDP) with existing directed clustering methods\
   Iter_Convg_Rand: visualize how the iterative algorithm (Algorithm 4 in our paper) updates the DSBM parameters\
    cluster_algs/f_IT_MLE_sc: MLE-SC (Algorithm 4 + Algorithm 1) in our paper\
  cluster_algs/f_IT_MLE_bm: MLE-SDP (Algorithm 4 + Algorithm 3) in our paper (can replace Algorithm 3 with Algorithm 2 if you prefer using SDP solver)\
```

* Cite our paper:
<pre>
  @misc{MLE_DSBM,
        title={Maximum Likelihood Estimation on Stochastic Blockmodels for Directed Graph Clustering}, 
        author={Mihai Cucuringu and Xiaowen Dong and Ning Zhang},
        year={2024},
        eprint={2403.19516},
        archivePrefix={arXiv},
        primaryClass={stat.ML}
        }
<pre>
