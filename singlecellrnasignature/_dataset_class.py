"""dataset_class.py

This module defines the `DatasetscRNASeqSignature` base class for working with namespace-aware datasets.
"""

import datalair

class DatasetscRNASeqSignature(datalair.Dataset):
    """Datalair Dataset class for all Datasets in this package."""

    def __init__(self):
        """Initialize this dataset class as a datalair.Dataset class with namespace `DatasetscRNASeqSignature`."""
        super().__init__(namespace="DatasetscRNASeqSignature")