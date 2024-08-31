[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# A Bilevel Optimization Approach for a Class of Combinatorial Problems with Disruptions and Probing
This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The purpose of this repository is to share the codes, instances, and results used in the paper "A Bilevel Optimization Approach for a Class of Combinatorial Problems with Disruptions and Probing", authored by Leonardo Lozano and Juan Borrero.

## Cite
To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.
[https://doi.org/10.1287/ijoc.2024.0629](https://doi.org/10.1287/ijoc.2024.0629)

[https://doi.org/10.1287/ijoc.2024.0629.cd](https://doi.org/10.1287/ijoc.2024.0629.cd)

Below is the BibTex for citing this snapshot of the repository.
```
@article{combinatorialProbingGITHUB,
  author =        {Lozano, Leonardo and Borrero, Juan},
  publisher =     {INFORMS Journal on Computing},
  title =         {A Bilevel Optimization Approach for a Class of Combinatorial Problems with Disruptions and Probing},
  year =          {2024},
  doi =           {10.1287/ijoc.2024.0629.cd},
  note =          {Available for download at https://github.com/INFORMSJoC/2024.0629},
}  
```

## Dataset 
The "data" folder inside the "combinatorialProbing.zip" file contains all the datasets used in the paper, including the shortest path networks in DIMACS format and the project selection (knapsack) problem instances. For the project selection problem instances, the files display the nominal profits, the weigths, and the right-hand-side coefficient. The name of the files present the size (network size for shortest path and number of items for knapsack) as well as the random seed used for its generation (a number between 0 and 9 corresponding to each one of the ten replicas per configuration).   

## Results 
The results folder inside the "combinatorialProbing.zip" file the raw outputs from the execution of the models. The file "results2024" contains all the organized data extracted from the outputs, as well as all the tables presented in the manuscript. 

## Replicating
The "combinatorialProbing.zip" file contains all the source code in Java with the implementation of all the optimization models and heuristics described in the manuscript. 

To run the code and fully replicate the experiments, you will need to make sure that you have a valid license for <code>CPLEX</code>. The project contains two packages named "KP" and "SP". Package "KP" contains all the source code related to project selection problems modeled as knapsack problems while package "SP" contains all the source code related to the shorstest path problems. 

