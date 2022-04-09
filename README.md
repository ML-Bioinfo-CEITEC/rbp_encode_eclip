# RBP Encode eCLIP

This repository contains RBP eCLIP data from Encode in "interval" format (for space eficiency reasons), notebook and a python script:
- [create_csv_datasets.ipynb](create_csv_datasets.ipynb) for reproducing creation of "interval" dataset and
- [add_fasta.py](add_fasta.py) for adding fasta sequences to the dataset.

## Usage

```python
# clone this repository
git clone https://github.com/ML-Bioinfo-CEITEC/rbp_encode_eclip.git
cd rbp_encode_eclip

# create virtual environment
virtualenv venv --python=python3.8
source venv/bin/activate

# install dependencies
pip install genomic_benchmarks

# create final dataset
python add_fasta.py
```