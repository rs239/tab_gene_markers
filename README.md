# README

## tab_query.py

`tab_query.py` is a command-line tool to query split-intein Gal4 gene-markers for Fly Cell Atlas clusters. Exactly two markers per cluster are returned. More than one set of options is returned with teh top hit being the strongest. The script will automatically download the relevant data in the user's directory.

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

- `--species`: The species to query. The only available option right now is 'dm' (Drosophila melanogaster). The default value is 'dm'.
- `--organ`: The organ to query. The available options depend on the chosen species. Use the `--describe` option to see the available choices.
- `--cluster`: The cluster to query. The available options depend on the chosen species and organ. Use the `--describe` option to see the available choices.
- `--params`: A list of key:value pairs of parameters. These parameters will be used to filter the data before the query. Each pair should be separated by a space. See the the DEFAULT_PARAMS variable in code.

## Examples

To query the gene expression data for the 4th age group of male Drosophila melanogaster antennae, type the following command:

```
python tab_query.py  --cluster "enteroblast" --organ Gut
```


## Dependencies

`tab_query.py` requires the following dependencies:

- argparse
- pandas
- numpy
- scanpy
- anndata

You can install these dependencies using pip. For example:

```
pip install argparse pandas numpy scanpy anndata
```

Note: If you use `tab_query.py` in your research, please consider citing the following preprint and paper: [split-intein Gal4 provides intersectional genetic labeling that is fully repressible by Gal80](https://www.biorxiv.org/content/10.1101/2023.03.24.534001)