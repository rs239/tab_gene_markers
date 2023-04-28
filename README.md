# README

![image](https://user-images.githubusercontent.com/6614489/235061911-845f5851-ecd7-4c22-925f-74a924a975ea.png)

## tab_query.py

`tab_query.py` is a command-line tool to query split-intein Gal4 gene-markers for Fly Cell Atlas clusters. 
 
  - Exactly two markers per cluster are returned. This is important for the split-intein Gal4 protocol. 
  - More than one set of results is returned, with the top hit being the strongest. 
  - The script will automatically download the necessary data in the user's directory. By default, two files will be downloaded in a data/ subfolder. You will need approx 4.5G of disk space.

## Usage

To use `tab_query.py`, open a terminal and navigate to the directory containing the `tab_query.py` file. Then, type the following command:

```
python tab_query.py --species [species] --organ [organ] --cluster [cluster] --params [params]
```
or

```
python tab_query.py --describe
```

The arguments for the command are as follows:

- `--species`: The species to query. The only available option right now is 'dm' (Drosophila melanogaster) and the default value is 'dm'.
- `--organ`: The organ to query. These are the organs corresponding to which the bulk RNA-seq data will be used. This needs to match (in your judgement) the single-cell cluster you specify. Use the `--describe` option to see the available choices.
- `--cluster`: The cluster to query. These are currently taken from the Fly Cell Atlas annotation. Use the `--describe` option to see the available choices.
- `--params`: A list of key:value pairs of parameters. These parameters will be used to filter the data before the query. Each pair should be separated by a space. See the the DEFAULT_PARAMS variable in code.
- `--describe`: Provides the list of organs and Fly Cell Atlas currently accepted. 

## Examples

To query the gene markers for enteroblast cells in the fly gut:

```
python tab_query.py  --cluster "enteroblast" --organ Gut
```

To see the list of organs and single-cell clusters:

```
python tab_query.py  --describe
```


## Dependencies

`tab_query.py` requires the following dependencies:

- pandas
- numpy
- scanpy
- anndata

You can install these dependencies using pip. For example:

```
pip install pandas numpy scanpy anndata
```

Note: If you use `tab_query.py` in your research, please consider citing the following preprint and paper: [split-intein Gal4 provides intersectional genetic labeling that is fully repressible by Gal80](https://www.biorxiv.org/content/10.1101/2023.03.24.534001)
