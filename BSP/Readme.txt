First, use 

>make clean
>make all 

to compile.

For BSP without copula, use 

>./bsp BSP_NC <datafilename> <outputfilename>

The output file shows the estimated joint distribution, each row corresponds to a subregion, its density and the number of points in the subregion.


For BSP with copula (usually used when dim>3)

>./bsp BSP_C <datafilename> <ofilename> <level1> <level2>

level1 and level 2 are optional. level1 is the maximum number of levels for joint distribution and level2 is the maximum number of levels for marginal distribution.
There are d+1 output files: "ofilename_0.txt",..., "ofilename_(d-1).txt" and "ofilename_Big.txt"(thus, 
you may not want to include '.txt' in the filename). 
The first d files are the MAP partitions for each dimension(represent marginal distribution), while
the last file is the MAP partition for the joint distribution using the copula-transformed data. The first part of "ofilename_Big.txt" is the binary partition result from BSP and the second part is the one after copula transformation.

Examples:

## for non copula case:

>./bsp BSP_NC data/data5normal1000_1.txt result/5normal_MAP_2_Jan.txt

## for copula case:

>./bsp BSP_C data/datamixnormal10000_0.txt result/mixnormal_MAP 200 60




## 2019/12/20 By Yiliang

2000 / 2400 / 3000

./bsp BSP_NC data/low_dim_alpha/dim10/dim10_size100000_bsp.txt result/low_dim_alpha/dim10/dim10_size100000_bsp_result.txt 3000


./bsp density -n data/low_dim_alpha/dim10/dim10_size100000_bsp.txt result/low_dim_alpha/dim10/dim10_size100000_bsp_result.txt_Big.txt
