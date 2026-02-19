"""raw.py

This module defines the single cell TNA sequencing datasets as downloaded from public repositories
(e.g. gene expression omnibus, or arrayexpress).
"""


import subprocess
from singlecellrnasignature._dataset_class import DatasetscRNASeqSignature as _Dataset
from datalair import (Lair as _Lair,
                      download_supplementary_from_geo as _download_geo,
                      download_files_from_arrayexpress as _download_ae,
                      download_file as _download_file)


class AziziSingleCellMapDiverse2018(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE114727", output_dir)


class BeckerSinglecellAnalysesDefine2022(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE201349", output_dir)


class BiTumorImmuneReprogramming2021(_Dataset):

    shell_command = None

    def derive(self, lair: _Lair) -> None:
        assert self.shell_command is not None,\
            """You must first provide the shell command for download from the website of the Braod institute!
            Go to https://singlecell.broadinstitute.org/single_cell/study/SCP1288/tumor-and-immune-reprogramming-during-immunotherapy-in-advanced-renal-cell-carcinoma#/,
            , click on the "Download" tab, and click on "Bulk download" to get a curl command that is only valid for 30 minutes.
            Also, you need to be signed in! And remember that each link can only be used once!"""
        output_dir = lair.get_path(self)
        result = subprocess.run(self.shell_command, capture_output=True, text=True, shell=True, cwd=output_dir)
        print("=== BEGIN CURL OUTPUT ===/n{}\n=== END CURL OUTPUT ===".format(result))


class BiermannDissectingTreatmentnaiveEcosystem2022(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE200218", output_dir)


class BorcherdingMappingImmuneEnvironment2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE121638", output_dir)


class BraunProgressiveImmuneDysfunction2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("Need to request access (dbGaP: phs002252.v1.p1)!")


class ChanSignaturesPlasticityMetastasis2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("Need to go to HTAN Network!")


class ChengPancancerSinglecellTranscriptional2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE154763", output_dir)


class CheSinglecellAtlasLiver2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE178318", output_dir)


class DuranteSinglecellAnalysisReveals2020(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE139829", output_dir)


class JerbyArnonCancerCellProgram2018(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE115978", output_dir)


class KhaliqRefiningColorectalCancer2022(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE200997", output_dir)


class KimSinglecellRNASequencing2020(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE131907", output_dir)


class KrishnaSinglecellSequencingLinks2021(_Dataset):
    """Data from SRA: SRP308561; BioProject: PRJNA705464
    https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP308561"""

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        urls = [
            "https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/bulk_multiregional_cohort_counts.csv",
            "https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/bulk_multiregional_cohort_metadata.csv",
            "https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/bulk_sortedpops_counts.csv",
            "https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/bulk_sortedpops_metadata.csv",
            "https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/ccRCC_6pat_cell_annotations.txt",
            "https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/ccRCC_6pat_Seurat",
            "https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/ccRCC_regions.txt",
            "https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/ccRCC_TCRs.txt",
            "https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/000190/SRZ190804/File_descriptions.txt"
        ]
        for url in urls:
            _download_file(url, output_dir.joinpath(url.split("/")[-1]))


class LeaderSinglecellAnalysisHuman2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE154826", output_dir)


class LiLiquidBiopsybasedSinglecell2019(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("BioProject accession PRJNA554445")


class LuSinglecellAtlasMulticellular2022(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE149614", output_dir)


class MaynardTherapyInducedEvolutionHuman2020(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("BioProject accession PRJNA591860")


class PelkaSpatiallyOrganizedMulticellular2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE178341", output_dir)


class PomboantunesSinglecellProfilingMyeloid2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE163120", output_dir)


class PuSinglecellTranscriptomicAnalysis2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE184362", output_dir)


class QianPancancerBlueprintHeterogeneous2020a(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_ae("E-MTAB-8107", output_dir)


class SharmaOncofetalReprogrammingEndothelial2020(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE156625", output_dir)


class VazquezOvarianCancerMutational2022(_Dataset):
    """The h5 file is actually a h5ad file"""

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE180661", output_dir)


class WuSinglecellProfilingTumor2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE148071", output_dir)


class WuSinglecellSpatiallyResolved2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE176078", output_dir)


class XuSinglecellRNASequencing2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE180286", output_dir)


class ZhangLongitudinalSinglecellRNAseq2022(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE165897", output_dir)


class ZhangDissectingEsophagealSquamouscell2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE160269", output_dir)


class ZhangSinglecellAnalysesReveal2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE169246", output_dir)


class ZhangSinglecellAnalysisReveals2022(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE215120", output_dir)


class ZhengSinglecellTranscriptomicProfiling2022(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE161277", output_dir)


class ZilionisSingleCellTranscriptomicsHuman2019(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE127465", output_dir)


########

class LiuSinglecellAtlasReveals2025(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE243013", output_dir)

########


class BanchereauMolecularDeterminantsResponse2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("GO to: https://ega-archive.org/studies/EGAS00001004343")


class HugoGenomicTranscriptomicFeatures2016(_Dataset):

    def derive(self, lair: _Lair) -> None:
        output_dir = lair.get_path(self)
        _download_geo("GSE78220", output_dir)


class KimComprehensiveMolecularCharacterization2018(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("Go to: https://www.ebi.ac.uk/ena/browser/view/PRJEB25780")


class LiuIntegrativeMolecularClinical2019(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError(
            """All reasonable requests for raw and analyzed data and materials will be promptly reviewed by the senior
            authors to determine whether the request is subject to any intellectual property or confidentiality
            obligations. Patient-related data not included in the paper may be subject to patient confidentiality.
            Any data and materials that can be shared will be released via a material transfer agreement. All analyzed
            sequencing data are in supplementary tables or data available online. Raw sequencing data are available in
            dbGaP (accession number phs000452.v3.p1)."
            """)


class MariathasanTGFvAttenuatesTumour2018(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("Go to: https://www.ebi.ac.uk/ega/studies/EGAS00001002556")


class McdermottClinicalActivityMolecular2018(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("Go to: https://www.ebi.ac.uk/ega/search/site/EGAS00001002928 and for clinical data to: https://clinicalstudydatarequest.com/")


class MiaoGenomicCorrelatesResponse2018(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("Go to: https://www.science.org/doi/full/10.1126/science.aan5951")


class PatilIntratumoralPlasmaCells2022(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("""Raw and processed transcriptomic data, relevant mutation status, and limited clinical
        data has been deposited at the European Genome-phenome Archive (EGA), which is hosted by the EBI and the CRG,
        under accession number EGA: EGAS00001005013. Additional clinical data is available via request from
        vivli.org.""")


class RiazTumorMicroenvironmentEvolution2017(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError("""Raw and processed transcriptomic data, relevant mutation status, and limited clinical
        data has been deposited at the European Genome-phenome Archive (EGA), which is hosted by the EBI and the CRG,
        under accession number EGA: EGAS00001005013. Additional clinical data is available via request from
        vivli.org.""")


class VanallenGenomicCorrelatesResponse2015(_Dataset):

    def derive(self, lair: _Lair) -> None:
        NotImplementedError(
            """No automatic download possible! Go to https://www.science.org/doi/full/10.1126/science.aan5951 and
            downloa dthe supplementary materials manually""")


class LiMappingSinglecellTranscriptomes2022(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError(
            """The genome sequence data reported in this paper is available at the European
            Genome-Phenome Archive: EGAD00001008029 for whole-exome sequencing data, EGAD00001008030 for
            the single cell RNA sequencing data, and EGAD00001008781 for the spatial transcriptomic data.
            Our single cell RNA sequencing and spatial transcriptomics data are available to download as h5ad objects
            in Mendeley Data: https://doi.org/10.17632/g67bkbnhhg.1.""")


class BassezSinglecellMapIntratumoral2021(_Dataset):

    def derive(self, lair: _Lair) -> None:
        raise NotImplementedError(
            """Raw sequencing reads of all single-cell experiments (scRNA-seq, scTCR-seq and CITE-seq) have been
            deposited in the European Genome-phenome Archive (EGA) under study no. EGAS00001004809 (with a summary of
            the BioKey study and patient characteristics) and with data accession no. EGAD00001006608 (to access the
            data itself under restricted access). Requests for accessing raw sequencing reads will be reviewed by
            the UZLeuven-VIB data access committee. Any data shared will be released via a Data Transfer Agreement
            that will include the necessary conditions to guarantee protection of personal data (according to European
            GDPR law). Alternately, a download of the read count data per individual patient is publicly available at
            http://biokey.lambrechtslab.org. The publicly available gnomAD database (https://gnomad.broadinstitute.org)
            was used to filter tumor exome-seq data for somatic mutations and calculate tumor mutation burden.
            Raw sequencing reads of all exome and low-coverage whole-genome sequencing experiments are also provided
            under EGAS00001004809.""")
