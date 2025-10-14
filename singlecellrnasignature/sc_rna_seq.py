import polars as pl
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import re
import gzip
import shutil
import tarfile
import tempfile

from pathlib import Path as _Path
from scipy.io import mmread as _mmread
from scipy.sparse import csr_matrix as _csr_matrix
from datalair import Lair as _Lair
from singlecellrnasignature import raw_sc_rna_seq as _raw
from singlecellrnasignature.dataset_class import DatasetscRNASeqSignature as _Dataset


class AziziSingleCellMapDiverse2018Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.AziziSingleCellMapDiverse2018()
        lair.safe_derive(ds)
        filepaths = lair.get_dataset_filepaths(ds)

        with tempfile.TemporaryDirectory() as tmpdir:
            extract_dir = _Path(tmpdir)

            with tarfile.open(filepaths["GSE114727_RAW.tar"], 'r:*') as tar:  # 'r:*' auto-detects compression
                tar.extractall(path=extract_dir)

            mtx_files = list(filter(lambda x: x.name.endswith("matrix.mtx.gz"), extract_dir.iterdir()))
            counts_files = list(filter(lambda x: x.name.endswith("counts.csv.gz"), extract_dir.iterdir()))

            adatas = []

            for count_file in counts_files:
                print(count_file)
                df = pd.read_csv(count_file, index_col=0).fillna(0)
                adata = ad.AnnData(df)
                metadata = count_file.name.split("_")[:3]
                adata.obs[["geo_id", "patient", "tissue"]] = metadata
                adatas.append(adata)

            for mtx_file in mtx_files:
                data = _mmread(mtx_file).T.tocsr()
                metadata = mtx_file.name.split("_")[:3]
                barcodes = pd.read_csv(mtx_file.joinpath("..", "_".join(metadata[:3] + ["barcodes.tsv.gz"])).resolve(), sep="\t", header=None)[0]
                genes = pd.read_csv(mtx_file.joinpath("..", "_".join(metadata[:3] + ["genes.tsv.gz"])).resolve(), sep="\t", header=None).rename(columns={0: "ensembl_id", 1: "gene_name"}).set_index("ensembl_id")
                adata = ad.AnnData(X=data, var=genes)
                adata.obs["original.barcode"] = list(barcodes)
                adata.obs[["geo_id", "patient", "tissue"]] = metadata
                adatas.append(adata)

        adata = ad.concat(adatas, join="outer", axis=0)
        adata.obs.reset_index(inplace=True)
        adata.write(output_dir.joinpath("adata.h5ad"))


class BeckerSinglecellAnalysesDefine2022Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.BeckerSinglecellAnalysesDefine2022()
        lair.safe_derive(ds)
        filepaths = lair.get_dataset_filepaths(ds)

        with tempfile.TemporaryDirectory() as tmpdir:
            extract_dir = _Path(tmpdir)
            with tarfile.open(filepaths["GSE201349_RAW.tar"], 'r:*') as tar:  # 'r:*' auto-detects compression
                tar.extractall(path=extract_dir)
            files = sorted(list(filter(lambda x: x.name.startswith("GSM6061"), extract_dir.iterdir())))
            prefixes = {"_".join(file.name.split("_")[:-1]) for file in files}
            adatas = []
            for prefix in prefixes:
                data = _mmread(extract_dir.joinpath(prefix+"_matrix.mtx.gz")).T.tocsr()
                barcodes = pd.read_csv(extract_dir.joinpath(prefix+"_barcodes.tsv.gz"), header=None, sep="\t")[0]
                features = (pd.read_csv(extract_dir.joinpath(prefix+"_features.tsv.gz"), header=None, sep="\t")
                            .rename(columns={0: "ensembl_id", 1: "gene_name", 2: "gene_type"})
                            .set_index("ensembl_id"))
                adata = ad.AnnData(data)
                adata.obs_names = barcodes
                adata.obs[["geo_id", "sample"]] = (prefix.split("_")[0], "_".join(prefix.split("_")[1:]))
                adata.var = features
                adatas.append(adata)
        adata = ad.concat(adatas, axis=0, join="outer")
        adata.write(output_dir.joinpath("adata.h5ad"))


class BiTumorImmuneReprogramming2021Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.BiTumorImmuneReprogramming2021()
        lair.safe_derive(ds)
        filepaths = lair.get_dataset_filepaths(ds)

        dir_path = _Path(filepaths["SCP1288"].joinpath("expression", "60c76a18771a5b0ba10ea91b"))
        data = _mmread(dir_path.joinpath("matrix.mtx")).T.tocsr()
        barcodes = pd.read_csv(dir_path.joinpath("barcodes.tsv"), sep="\t", header=None)[0]
        genes = pd.read_csv(dir_path.joinpath("genes.tsv"), sep="\t", header=None)[0]
        adata = ad.AnnData(data)
        adata.obs_names = barcodes
        adata.var_names = genes

        adata.write(output_dir.joinpath("adata.h5ad"))


class BiermannDissectingTreatmentnaiveEcosystem2022Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.BiermannDissectingTreatmentnaiveEcosystem2022()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        data = _mmread(filepaths["GSE200218_sc_sn_counts.mtx.gz"]).T.tocsr()
        metadata = pd.read_csv(filepaths["GSE200218_sc_sn_metadata.csv.gz"], index_col=0)
        genes = pd.read_csv(filepaths["GSE200218_sc_sn_gene_names.csv.gz"], index_col=0)
        mask = (genes.index != "1-Mar") & (genes.index != "2-Mar")
        data = data[:, mask]
        genes = genes.loc[mask]
        adata = ad.AnnData(data)
        adata.obs = metadata
        adata.var = genes

        adata.write(output_dir.joinpath("adata.h5ad"))


class BorcherdingMappingImmuneEnvironment2021Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.BorcherdingMappingImmuneEnvironment2021()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        with tempfile.TemporaryDirectory() as tmpdir:
            extract_dir = _Path(tmpdir)
            with tarfile.open(filepaths["GSE121638_RAW.tar"], 'r:*') as tar:  # 'r:*' auto-detects compression
                tar.extractall(path=extract_dir)

            filepaths = sorted(list(extract_dir.glob("*GSM344084*")))
            prefixes = {"_".join(filepath.name.split("_")[:-1]) for filepath in filepaths}
            adatas = []
            for prefix in prefixes:
                data = _mmread(extract_dir.joinpath("_".join([prefix, "matrix.mtx.gz"]))).T.tocsr()
                genes = (pd.read_csv(extract_dir.joinpath("_".join([prefix, "genes.tsv.gz"])), sep="\t", header=None)
                        .rename(columns={0: "ensembl_id", 1: "gene_name"}))
                barcodes = pd.read_csv(
                    extract_dir.joinpath("_".join([prefix, "barcodes.tsv.gz"])),
                    sep="\t", header=None)
                adata = ad.AnnData(data)
                adata.var = genes
                adata.var.set_index("ensembl_id", inplace=True)
                adata.obs_names = barcodes.index
                adatas.append(adata)
            adata = ad.concat(adatas, axis=0, join="outer")

        adata.write(output_dir.joinpath("adata.h5ad"))


class ChengPancancerSinglecellTranscriptional2021Adata(_Dataset):
    # uuid = datalair.UUID("bf59bdb576656b93")

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _raw_dataset = _raw.ChengPancancerSinglecellTranscriptional2021()
        lair.safe_derive(_raw_dataset)
        filepaths = lair.get_dataset_filepaths(_raw_dataset)
        all_keys = list(filepaths.keys())
        cancer_types = [key.split("_")[1] for key in all_keys]

        adatas = []
        for i, cancer_type in enumerate(cancer_types):
            df = pd.read_csv(filepaths["GSE154763_{}_{}.csv.gz".format(cancer_type, "normalized_expression")],
                             index_col=0)
            metadata = pd.read_csv(filepaths["GSE154763_{}_{}.csv.gz".format(cancer_type, "metadata")], index_col=0)
            assert np.all(df.index == metadata.index)
            adata = ad.AnnData( _csr_matrix(np.array(df.values, dtype=np.float64)), obs=metadata)
            adata.var_names = df.columns
            adata.obs["cancer_type"] = cancer_type
            adatas.append(adata)

        adata = ad.concat(adatas, axis="obs", join="outer")
        adata.write_h5ad(output_dir.joinpath("adata.h5ad"))


class DuranteSinglecellAnalysisReveals2020Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.DuranteSinglecellAnalysisReveals2020()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        with tempfile.TemporaryDirectory() as tmpdir:
            extract_dir = _Path(tmpdir)
            with tarfile.open(filepaths["GSE139829_RAW.tar"], 'r:*') as tar:  # 'r:*' auto-detects compression
                tar.extractall(path=extract_dir)
            filepaths = sorted(list(extract_dir.glob("*")))
            prefixes = {"_".join(filepath.name.split("_")[:-1]) for filepath in filepaths}
            adatas = []
            for prefix in prefixes:
                data = _mmread(extract_dir.joinpath("_".join([prefix, "matrix.mtx.gz"]))).T.tocsr()
                genes = pd.read_csv(extract_dir.joinpath("_".join([prefix, "genes.tsv.gz"])), sep="\t",
                                    header=None).rename(columns={0: "ensembl_id", 1: "gene_name"})
                barcodes = pd.read_csv(extract_dir.joinpath("_".join([prefix, "barcodes.tsv.gz"])), sep="\t",
                                       header=None)
                adata = ad.AnnData(data)
                adata.var = genes
                adata.var.set_index("ensembl_id", inplace=True)
                adata.obs_names = barcodes.index
                adata.obs[["geo_id", "sample"]] = prefix.split("_")
                adatas.append(adata)
        adata = ad.concat(adatas, axis=0, join="outer")
        adata.write_h5ad(output_dir.joinpath("adata.h5ad"))


class JerbyArnonCancerCellProgram2018Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.JerbyArnonCancerCellProgram2018()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        data = pd.read_csv(filepaths["GSE115978_counts.csv.gz"], index_col=0).T
        cell_annotations = pd.read_csv(filepaths["GSE115978_cell.annotations.csv.gz"], index_col=0)
        adata = ad.AnnData(data)
        adata.obs = cell_annotations

        adata.write_h5ad(output_dir.joinpath("adata.h5ad"))


class KhaliqRefiningColorectalCancer2022Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.KhaliqRefiningColorectalCancer2022()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        df = pl.read_csv(filepaths["GSE200997_GEO_processed_CRC_10X__raw_UMI_count_matrix.csv.gz"])
        data =  _csr_matrix(df.select(df.columns[1:]).to_numpy()).T
        genes = np.ndarray.flatten(df.select("").to_numpy())
        cell_annotations = pd.read_csv(filepaths["GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"], index_col=0)
        adata = ad.AnnData(data)
        adata.obs = cell_annotations
        adata.var_names = genes
        adata.write(output_dir.joinpath("adata.h5ad"))


class KimSinglecellRNASequencing2020Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.KimSinglecellRNASequencing2020()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        data = pl.read_csv(filepaths["GSE131907_Lung_Cancer__raw_UMI_matrix.txt.gz"], separator="\t")
        genes = data.select("Index")
        data = data.drop("Index")
        data = data.to_numpy()
        data =  _csr_matrix(data).T
        adata = ad.AnnData(data)
        adata.var_names = genes["Index"]
        adata.obs = pd.read_csv(filepaths["GSE131907_Lung_Cancer_cell_annotation.txt.gz"], sep="\t")

        adata.write(output_dir.joinpath("adata.h5ad"))


class KrishnaSinglecellSequencingLinks2021Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.KrishnaSinglecellSequencingLinks2021()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        assert filepaths["ccRCC_6pat_Seurat.h5ad"].exists(), "You must first convert the RDS to h5ad manually!"
        adata = ad.read(filepaths["ccRCC_6pat_Seurat.h5ad"])
        del adata._raw.var["_index"]
        adata.write(output_dir.joinpath("adata.h5ad"))


class LeaderSinglecellAnalysisHuman2021Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.LeaderSinglecellAnalysisHuman2021()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        with tempfile.TemporaryDirectory() as tmpdir:
            extract_dir = _Path(tmpdir)
            for filename in filepaths.keys():
                if filename == "GSE154826_sample_annots.csv.gz":
                    continue
                sub_extract_dir = extract_dir.joinpath(filename.split(".")[0])
                sub_extract_dir.mkdir()
                with tarfile.open(filepaths[filename], mode="r:*") as tar:
                    tar.extractall(path=extract_dir.joinpath(sub_extract_dir))

            adatas = []
            for sub_dir in extract_dir.iterdir():
                files = sorted(list(sub_dir.iterdir()))
                barcodes = pd.read_csv(list(filter(lambda x: "barcodes" in x.name, files))[0], sep="\t", header=None)
                genes = pd.read_csv(list(filter(lambda x: "features" in x.name, files))[0], sep="\t", header=None)
                data = _mmread(list(filter(lambda x: "matrix" in x.name, files))[0]).tocsr().T
                adata = ad.AnnData(data)
                adata.obs = barcodes
                adata.obs["batch"] = int(sub_dir.name.split("_")[-1])
                adata.var = genes
                adata.var.set_index(0, inplace=True)
                adatas.append(adata)

        adata = ad.concat(adatas)
        adata.var.rename(columns={0: "barcode"}, inplace=True)
        adata.write(output_dir.joinpath("adata.h5ad"))


class LuSinglecellAtlasMulticellular2022Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        ds = _raw.LuSinglecellAtlasMulticellular2022()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir = _Path(tmp_dir)
            count_file = tmp_dir.joinpath("count.txt")
            with gzip.open(filepaths["GSE149614_HCC.scRNAseq.S71915.count.txt.gz"], "rb") as f_in:
                with open(count_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            fixed_count_file = tmp_dir.joinpath("count-fixed.txt")
            with open(count_file, "r") as f_in, open(fixed_count_file, "w") as f_out:
                first_line = f_in.readline()
                f_out.write("gene\t" + first_line)
                shutil.copyfileobj(f_in, f_out)
            data = pl.read_csv(fixed_count_file, separator="\t")
            genes = data.select("gene")
            data =  _csr_matrix(data.drop("gene").to_numpy()).T
            cell_annotations = pd.read_csv(filepaths["GSE149614_HCC.metadata.updated.txt.gz"], sep="\t").set_index("Cell")
            adata = ad.AnnData(data)
            adata.obs = cell_annotations
            adata.var_names = genes["gene"]
        adata.write(output_dir.joinpath("adata.h5ad"))


class PelkaSpatiallyOrganizedMulticellular2021Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.PelkaSpatiallyOrganizedMulticellular2021()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        adata = sc.read_10x_h5(filepaths["GSE178341_crc10x_full_c295v4_submit.h5"])
        metatables = pd.read_csv(filepaths["GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"])
        cluster = pd.read_csv(filepaths["GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"])
        assert np.all(cluster["sampleID"]==metatables["cellID"])
        assert set(metatables.columns).intersection(set(cluster.columns))==set()
        obs = pd.concat([metatables, cluster], axis=1)
        assert np.all(adata.obs_names==metatables["cellID"])
        adata.obs = obs

        adata.write(output_dir.joinpath("adata.h5ad"))


class PuSinglecellTranscriptomicAnalysis2021Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.PuSinglecellTranscriptomicAnalysis2021()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        with tempfile.TemporaryDirectory() as tmp_dir:
            extract_dir = _Path(str(tmp_dir))
            with tarfile.open(filepaths["GSE184362_RAW.tar"], "r") as tar:
                tar.extractall(extract_dir)
            filepaths = sorted(list(extract_dir.glob("*")))
            prefixes = {"_".join(filepath.name.split("_")[:-1]) for filepath in filepaths}
            adatas = []
            for prefix in prefixes:
                data = _mmread(extract_dir.joinpath("_".join([prefix, "matrix.mtx.gz"]))).T.tocsr()
                genes = pd.read_csv(extract_dir.joinpath("_".join([prefix, "features.tsv.gz"])), sep="\t",
                                    header=None).rename(columns={0: "ensembl_id", 1: "gene_name"})
                barcodes = pd.read_csv(extract_dir.joinpath("_".join([prefix, "barcodes.tsv.gz"])), sep="\t",
                                       header=None)
                adata = ad.AnnData(data)
                adata.var = genes
                adata.var.set_index("ensembl_id", inplace=True)
                adata.obs_names = barcodes.index
                adata.obs[["geo_id", "patient", "biopsy_site"]] = prefix.split("_")
                adatas.append(adata)
            adata = ad.concat(adatas, axis=0, join="outer")
            assert np.all(adatas[0].var_names==adata.var_names)
            adata.var = adatas[0].var
            adata.write(output_dir.joinpath("adata.h5ad"))


class QianPancancerBlueprintHeterogeneous2020aAdata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.QianPancancerBlueprintHeterogeneous2020a()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        adatas = []
        sample_metadata = pd.read_csv(filepaths['E-MTAB-8107.sdrf.txt'], sep="\t")
        for name, path in filepaths.items():
            if not name.endswith(".counts.csv"):
                continue
            df = pl.read_csv(path)
            genes = df[:, 0]
            barcodes = df.columns[1:]
            data =  _csr_matrix(df[:, 1:]).T
            adata = ad.AnnData(data)
            adata.var_names = genes
            source_name = name.split(".")[0]
            if re.fullmatch(r"scrEXT[0-9]{3}", source_name):
                source_name = source_name[3:]
            adata.obs = pd.concat(len(barcodes)*[sample_metadata[sample_metadata["Source Name"]==source_name]], ignore_index=True)
            adatas.append(adata)
        adata = ad.concat(adatas, axis=0, join="outer")
        adata.write(output_dir.joinpath("adata.h5ad"))


class SharmaOncofetalReprogrammingEndothelial2020Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.SharmaOncofetalReprogrammingEndothelial2020()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        adata = ad.AnnData(_mmread(filepaths["GSE156625_HCCFmatrix.mtx.gz"]).T.tocsr())
        adata.obs = pd.read_csv(filepaths["GSE156625_HCCFbarcodes.tsv.gz"], sep="\t", header=None).rename(columns={0: "barcode"})
        adata.obs["cancer_type"] = "HCCF"
        adata.var = pd.read_csv(filepaths["GSE156625_HCCFgenes.tsv.gz"], sep="\t", header=None).rename(columns={0: "ensembl_id", 1: "gene_name"})

        bdata = ad.AnnData(_mmread(filepaths["GSE156625_HCCmatrix.mtx.gz"]).T.tocsr())
        bdata.obs = pd.read_csv(filepaths["GSE156625_HCCbarcodes.tsv.gz"], sep="\t", header=None).rename(columns={0: "barcode"})
        adata.obs["cancer_type"] = "HCC"
        bdata.var = pd.read_csv(filepaths["GSE156625_HCCgenes.tsv.gz"], sep="\t", header=None).rename(columns={0: "ensembl_id", 1: "gene_name"})

        assert np.all(adata.var_names==bdata.var_names)

        adata = ad.concat([adata, bdata], axis=1, join="outer")
        adata.var.set_index("ensembl_id", inplace=True)
        adata.write(output_dir.joinpath("adata.h5ad"))


class VazquezOvarianCancerMutational2022Adata(_Dataset):
    """Ovarian cancer scRNAseq dataset """

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.VazquezOvarianCancerMutational2022()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)
        adata = ad.read_h5ad(filepaths["GSE180661_matrix.h5"], backed="r")
        metadata = pd.read_csv(filepaths["GSE180661_GEO_cells.tsv.gz"], sep="\t", index_col=0)
        metadata = metadata.reindex(adata.obs.index)
        metadata = metadata.loc[adata.obs.index]
        assert np.all(adata.obs.index == metadata.index)
        adata.obs = metadata
        adata.write_h5ad(output_dir.joinpath("adata.h5ad"))


class ZhangSinglecellAnalysesReveal2021Adata(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)

        ds = _raw.ZhangSinglecellAnalysesReveal2021()
        lair.safe_derive(ds, overwrite=False)
        filepaths = lair.get_dataset_filepaths(ds)

        adata = ad.AnnData(_mmread(filepaths["GSE169246_TNBC_RNA.counts.mtx.gz"]).T.tocsr())
        adata.obs = pd.read_csv(filepaths["GSE169246_TNBC_RNA.barcode.tsv.gz"], header=None, sep="\t").rename(columns={0: "barcode"})
        adata.var_names = pd.read_csv(filepaths["GSE169246_TNBC_RNA.feature.tsv.gz"], header=None, sep="\t", index_col=0).index.astype(str).rename(None)
        adata.write(output_dir.joinpath("adata.h5ad"))
