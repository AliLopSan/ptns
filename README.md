# ptns
A repository for the generation of Perfect Transfer Networks (PTNs)

## What are Perfect Transfer Networks?
**Perfect Transfer Networks** also called PTNs [1](https://link.springer.com/article/10.1186/s13015-023-00242-2?utm_source=rct_congratemailt&utm_medium=email&utm_campaign=oa_20240206&utm_content=10.1186/s13015-023-00242-2) are a character-based model for the inference of Horizontal Gene Transfer events. How do they work? Well, it all starts with:
![](https://github.com/AliLopSan/ptns/blob/main/README_files/panel_1.gif)
### How do we find Horizontal Gene Transfers?
![](https://github.com/AliLopSan/ptns/blob/main/README_files/panel_2.gif)
### The intuition behind PTNs:
![](https://github.com/AliLopSan/ptns/blob/main/README_files/panel_3.gif)

## Dependencies
`ptns` has the following dependencies:

| package name | version |
|:------------:|:-------:|
|  asymmetree  |  2.2.1  |
|  tralda      |  1.1.0  |
|  requests    | 2.32.3  |

## Datastructures
In the `\src` folder, you will find the `datastructures.py` file. Since tree-based networks are just trees with additional edges, the main datastructure consists of a [tralda] TreeNode (https://github.com/david-schaller/tralda/tree/main) the main classes are:
- TB_Node
- TB_Network

## Tests
There are two ways to use the library, either on simulated or on real-life data. Both tests are contained in the `\tests` subfolder.

## References:
If you use this repository, please cite our work:
> López Sánchez, A., Lafond, M. Predicting horizontal gene transfers with perfect transfer networks. Algorithms Mol Biol 19, 6 (2024). https://doi.org/10.1186/s13015-023-00242-2
> López Sánchez, A., Lafond, M. (2024). Galled Perfect Transfer Networks. In: Scornavacca, C., Hernández-Rosales, M. (eds) Comparative Genomics. RECOMB-CG 2024. Lecture Notes in Computer Science(), vol 14616. Springer, Cham. https://doi.org/10.1007/978-3-031-58072-7_2
