from singlecellrnasignature import dataset_class
from datalair import Lair

def test_dataset_class(tmp_path):

    class Dataset(dataset_class.DatasetscRNASeqSignature):

        def derive(self, lair: Lair) -> None:
            output_dir = lair.get_path(self)
            open(output_dir.joinpath("empty.txt"), "w").close()

    lair = Lair(tmp_path.joinpath("lair"))
    lair.create()
    lair.assert_ok_satus()
    ds = Dataset()
    lair.safe_derive(ds)
    assert lair.get_path(ds).joinpath("empty.txt").exists()
    assert lair.get_path(ds).name=="DatasetscRNASeqSignature-Dataset"