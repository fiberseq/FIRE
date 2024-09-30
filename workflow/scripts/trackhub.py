#!/usr/bin/env python
import os
import defopt
import logging
from pathlib import Path
import numpy as np
from typing import Optional


HUB = """
hub {sample}-FIRE-fiberseq
shortLabel {sample}-FIRE-fiberseq
longLabel {sample}-FIRE-fiberseq
genomesFile genomes.txt
email mvollger.edu
descriptionUrl fire-description.html
"""

GENOMES = """
genome {ref}
trackDb trackDb.txt
"""


# transparentOverlay
PER_ACC_COMP = """
track {sample}-FIRE-percent-accessible
shortLabel {sample}-FIRE-percent-accessible
longLabel  {sample}-FIRE-percent-accessible
graphTypeDefault points
aggregate transparentOverlay
container multiWig
aggregate none 
showSubtrackColorOnUi on
type bigWig 0 1000
alwaysZero on
viewLimits 0:100
autoScale off
maxItems 100000
visibility full
maxHeightPixels 100:50:8
priority 1
yLineOnOff on
yLineMark 100
html fire-description.html
gridDefault on
"""

PER_ACC_TEMPLATE = """
    track {sample}-{hap}-FIRE-percent-accessible
    parent {sample}-FIRE-percent-accessible
    shortLabel {sample}-{hap}-FIRE-percent-accessible
    longLabel  {sample}-{hap}-FIRE-percent-accessible
    bigDataUrl {file}
    type bigWig
    visibility {viz}
    color {color}
"""


MULTI_WIG = """
track {sample}-{hap}-FIRE-coverage
parent {sample}-FIRE-coverage
longLabel {sample}-{hap}-FIRE-coverage
shortLabel {sample}-{hap}-FIRE-coverage
container multiWig
aggregate stacked
showSubtrackColorOnUi on
type bigWig 0 1000
autoScale off
alwaysZero on
viewLimits 0:{upper_coverage}
visibility {viz}
maxHeightPixels 100:50:8
html fire-description.html
priority 90
    
    track {sample}-{hap}-FIRE-accessible
    parent {sample}-{hap}-FIRE-coverage
    bigDataUrl {acc}
    type bigWig
    color 139,0,0
    
    track {sample}-{hap}-FIRE-linker
    parent {sample}-{hap}-FIRE-coverage
    bigDataUrl {link}
    type bigWig
    color 147,112,219
    
    track {sample}-{hap}-FIRE-nucleosome
    parent {sample}-{hap}-FIRE-coverage
    bigDataUrl {nuc}
    type bigWig
    color 169,169,169
    """


FIRE_SCORE_AND_FDR = """
track {sample}-FIRE-FDR
compositeTrack on
shortLabel {sample}-FIRE-FDR
longLabel {sample}-FIRE-FDR
type bigWig 0 1000
maxItems 100000
maxHeightPixels 100:50:1
alwaysZero on
priority 10
html fire-description.html

    track {sample}-log-FIRE-FDR
    parent {sample}-FIRE-FDR
    bigDataUrl {fdr}
    shortLabel {sample} -10log10 FDR
    longLabel {sample} -10log10 FDR
    autoScale on
    visibility hide
    yLineOnOff on
    yLineMark {y_line}
    gridDefault on
"""

TRACK_GROUPS = """
# grouping for fibers 
track {sample}-FIRE-fibers
compositeTrack on
shortLabel {sample}-FIRE-fibers
longLabel {sample}-FIRE-fibers
type bigBed 12 +
maxItems 100000
visibility dense
priority 80
html fire-description.html

# grouping for peaks
track {sample}-FIRE-peaks
compositeTrack on
shortLabel {sample}-FIRE-peaks
longLabel {sample}-FIRE-peaks
type bigBed 12 +
maxItems 100000
visibility dense
priority 30
html fire-description.html

    # track of unreliable regions just above the peak tracks
    track {sample}-unreliable-FIRE-coverage-regions
    parent {sample}-FIRE-peaks
    shortLabel {sample}-unreliable-FIRE-coverage-regions
    longLabel {sample}-unreliable-FIRE-coverage-regions
    type bigBed
    bigDataUrl bb/unreliable-coverage-regions.bb
    visibility dense
    priority 29

# grouping for coverage
track {sample}-FIRE-coverage
superTrack on show
shortLabel {sample}-FIRE-coverage
longLabel {sample}-FIRE-coverage
priority 90
html fire-description.html
"""

DECORATED = """
    track {sample}-{hap}-FIRE-fibers
    parent {sample}-FIRE-fibers
    shortLabel {sample}-{hap}-FIRE-fibers
    longLabel {sample}-{hap}-FIRE-fibers
    visibility dense
    type bigBed 12 +
    itemRgb On
    filterText.keywords {hap}
    bigDataUrl bb/fire-fibers.bb 
    decorator.default.bigDataUrl bb/fire-fiber-decorators.bb 
    decorator.default.filterValues.keywords 5mC,m6A,NUC,LINKER,FIRE
    decorator.default.filterValuesDefault.keywords LINKER,FIRE
"""


# type bigBed 6 + 4
FIRE_TEMPLATE = """
    track {sample}-narrow-FIRE-peaks
    parent {sample}-FIRE-peaks
    type bigNarrowPeak
    bigDataUrl {file}
    shortLabel {sample}-FIRE-peaks
    longLabel {sample}-FIRE-peaks
    visibility dense
    maxHeightPixels 50:50:1
"""

WIDE_TEMPLATE = """
    track {sample}-wide-FIRE-peaks
    parent {sample}-FIRE-peaks
    shortLabel {sample}-wide-FIRE-peaks
    longLabel {sample}-wide-FIRE-peaks
    type bigBed 
    bigDataUrl {file}
    visibility dense
    maxItems 100000
    priority 30
"""

HAP_TEMPLATE = """
    track {sample}-FIRE-hap-differences
    parent {sample}-FIRE-peaks
    type bigBed 9 +
    itemRgb on
    bigDataUrl {file}
    shortLabel {sample}-FIRE-hap-differences
    longLabel {sample}-FIRE-hap-differences
    visibility dense
    maxHeightPixels 25:25:1
"""


def generate_trackhub(
    trackhub_dir,
    ref,
    sample,
    ave_coverage,
):
    if ref == "T2Tv2.0":
        ref = "GCA_009914755.4"
    elif ref == "HG002v1.1":
        ref = "HG002v1.1.PAT"

    upper_coverage = int(ave_coverage + 5 * np.sqrt(ave_coverage))
    os.makedirs(f"{trackhub_dir}/", exist_ok=True)
    open(f"{trackhub_dir}/hub.txt", "w").write(HUB.format(sample=sample))
    open(f"{trackhub_dir}/genomes.txt", "w").write(GENOMES.format(ref=ref))
    trackDb = open(f"{trackhub_dir}/trackDb.txt", "w")
    trackDb.write(TRACK_GROUPS.format(sample=sample))

    for hap in ["all", "hap1", "hap2", "unk"]:
        if hap == "all":
            viz = "full"
        else:
            viz = "hide"
        # add coverage tracks
        if hap == "all":
            acc = f"bw/{hap}.fire.coverage.bw"
            nuc = f"bw/{hap}.nucleosome.coverage.bw"
            link = f"bw/{hap}.linker.coverage.bw"
            trackDb.write(
                MULTI_WIG.format(
                    acc=acc,
                    link=link,
                    nuc=nuc,
                    sample=sample,
                    hap=hap,
                    viz=viz,
                    upper_coverage=upper_coverage,
                )
            )

        if hap == "all":
            file = "bb/fire-peaks.bb"
            trackDb.write(FIRE_TEMPLATE.format(file=file, sample=sample))
            # add hap tracks
            file = "bb/hap_differences.bb"
            trackDb.write(HAP_TEMPLATE.format(file=file, sample=sample))
            file = "bb/fire-wide-peaks.bb"
            trackDb.write(WIDE_TEMPLATE.format(file=file, sample=sample))

        # add percent accessible tracks
        file = f"bw/{hap}.percent.accessible.bw"
        if hap == "all":
            color = "0,0,0"
            trackDb.write(PER_ACC_COMP.format(sample=sample))
        elif hap == "hap1":
            color = "0,0,255"
        elif hap == "hap2":
            color = "255,0,0"

        if hap != "unk":
            viz = "full" if hap != "all" else "hide"
            trackDb.write(
                PER_ACC_TEMPLATE.format(
                    sample=sample, hap=hap, file=file, color=color, viz=viz
                )
            )

        # new bin files
        if hap == "all":
            for z in ["H1", "H2", "UNK"]:
                trackDb.write(DECORATED.format(sample=sample, hap=z))

    # FDR scores
    trackDb.write(
        FIRE_SCORE_AND_FDR.format(
            sample=sample,
            fdr="bw/log_FDR.bw",
            score="bw/score.bw",
            y_line=-10 * np.log10(0.05),
        )
    )

    # done with track db
    trackDb.close()


def main(
    *,
    trackhub_dir: Optional[Path] = None,
    reference: Optional[str] = None,
    sample: Optional[str] = None,
    average_coverage: Optional[int] = 60,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger
    """
    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)
    generate_trackhub(trackhub_dir, reference, sample, average_coverage)
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
