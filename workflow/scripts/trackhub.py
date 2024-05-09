#!/usr/bin/env python
import os
import defopt
import logging
from pathlib import Path
import numpy as np
from typing import Optional


HUB = """
hub {sample}-fiberseq
shortLabel {sample}-fiberseq
longLabel {sample}-fiberseq
genomesFile genomes.txt
email mvollger.edu
"""

GENOMES = """
genome {ref}
trackDb trackDb.txt
"""


BW_COMP = """
track {sample}-{hap}-FDR
compositeTrack on
shortLabel {hap} FDR tracks 
longLabel {hap} FDR tracks
type bigWig 0 1000
autoScale off
viewLimits {FDR_min}:{FDR_max}
maxItems 100000
maxHeightPixels 50:50:1
"""

BW_TEMPLATE = """
    track FDR.{sample}.{hap}.{nm}
    parent {sample}-{hap}-FDR
    bigDataUrl {file}
    shortLabel FDR.{sample}.{hap}.{nm}
    longLabel FDR.{sample}.{hap}.{nm}
    type bigWig
    visibility {viz}
    priority {i}
    maxHeightPixels 50:50:1
"""

# transparentOverlay
PER_ACC_COMP = """
track {sample}-percent-accessible
shortLabel {sample}-percent-accessible
longLabel  {sample}-percent-accessible
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
"""

PER_ACC_TEMPLATE = """
    track {sample}-{hap}-percent-accessible
    parent {sample}-percent-accessible
    shortLabel {sample}-{hap}-percent-accessible
    longLabel  {sample}-{hap}-percent-accessible
    bigDataUrl {file}
    type bigWig
    visibility {viz}
    color {color}
"""

FIRE_SCORE_COMP = """
track {sample}-FIRE-score
shortLabel {sample}-FIRE-score
longLabel  {sample}-FIRE-score
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
priority 100
"""

FIRE_SCORE = """
    track {sample}-{hap}-FIRE-score
    parent {sample}-FIRE-score
    shortLabel {sample}-{hap}-FIRE-score
    longLabel  {sample}-{hap}-FIRE-score
    bigDataUrl {file}
    type bigWig
    visibility {viz}
    color {color}
"""


MULTI_WIG = """
track {sample}-{hap}-coverage
parent {sample}-coverage
longLabel {sample}-{hap}-coverage
shortLabel {sample}-{hap}-coverage
container multiWig
aggregate stacked
showSubtrackColorOnUi on
type bigWig 0 1000
autoScale off
alwaysZero on
viewLimits 0:{upper_coverage}
visibility full
maxHeightPixels 100:50:8
priority 90
    
    track {sample}-{hap}-accessible
    parent {sample}-{hap}-coverage
    bigDataUrl {acc}
    type bigWig
    color 139,0,0
    
    track {sample}-{hap}-linker
    parent {sample}-{hap}-coverage
    bigDataUrl {link}
    type bigWig
    color 147,112,219
    
    track {sample}-{hap}-nucleosome
    parent {sample}-{hap}-coverage
    bigDataUrl {nuc}
    type bigWig
    color 169,169,169
    """


FIRE_SCORE_AND_FDR = """
track {sample}-FIRE-FDR
compositeTrack on
shortLabel {sample}-FIRE-FDR
longLabel {sample}-FIRE-FDR
visibility full
type bigWig 0 1000
maxItems 100000
maxHeightPixels 100:50:1
alwaysZero on
priority 10

    track {sample}-log-fdr
    parent {sample}-FIRE-FDR
    bigDataUrl {fdr}
    shortLabel {sample} -10log10 FDR
    longLabel {sample} -10log10 FDR
    autoScale on
    visibility full
    yLineOnOff on
    yLineMark {y_line}
    gridDefault on
"""

TRACK_GROUPS = """
# grouping for fibers 
track {sample}-fibers
compositeTrack on
shortLabel {sample}-fibers
longLabel {sample}-fibers
type bigBed 12 +
maxItems 100000
visibility dense
priority 80

# grouping for peaks
track {sample}-peaks
compositeTrack on
shortLabel {sample}-peaks
longLabel {sample}-peaks
type bigBed 12 +
maxItems 100000
visibility dense
priority 30

    # track of unreliable regions just above the peak tracks
    track {sample}-unreliable-coverage-regions
    parent {sample}-peaks
    shortLabel {sample}-unreliable-coverage-regions
    longLabel {sample}-unreliable-coverage-regions
    type bigBed
    bigDataUrl bb/unreliable-coverage-regions.bb
    visibility dense
    priority 29

# grouping for coverage
track {sample}-coverage
superTrack on show
shortLabel {sample}-coverage
longLabel {sample}-coverage
priority 90
"""

DECORATED = """
    track {sample}-{hap}-fibers
    parent {sample}-fibers
    shortLabel {sample}-{hap}-fibers
    longLabel {sample}-{hap}-fibers
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
    track {sample}-FIRE-peaks
    parent {sample}-peaks
    type bigNarrowPeak
    bigDataUrl {file}
    shortLabel {sample}-FIRE-peaks
    longLabel {sample}-FIRE-peaks
    visibility dense
    maxHeightPixels 50:50:1
"""

WIDE_TEMPLATE = """
    track {sample}-{name}
    parent {sample}-peaks
    shortLabel {sample}-{name}
    longLabel {sample}-{name}
    type bigBed 
    bigDataUrl {file}
    visibility dense
    maxItems 100000
    priority 30
"""

HAP_TEMPLATE = """
    track {sample}-hap-differences
    parent {sample}-peaks
    type bigBed 9 +
    itemRgb on
    bigDataUrl {file}
    shortLabel {sample}-hap-differences
    longLabel {sample}-hap-differences
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

    upper_coverage = int(ave_coverage + 5 * np.sqrt(ave_coverage))
    os.makedirs(f"{trackhub_dir}/", exist_ok=True)
    open(f"{trackhub_dir}/hub.txt", "w").write(HUB.format(sample=sample))
    open(f"{trackhub_dir}/genomes.txt", "w").write(GENOMES.format(ref=ref))
    trackDb = open(f"{trackhub_dir}/trackDb.txt", "w")
    trackDb.write(TRACK_GROUPS.format(sample=sample))

    for hap in ["all", "hap1", "hap2", "unk"]:
        # add coverage tracks
        if hap != "unk":
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
                    upper_coverage=upper_coverage,
                )
            )

        if hap == "all":
            file = f"bb/FDR-FIRE-peaks.bb"
            trackDb.write(FIRE_TEMPLATE.format(file=file, sample=sample))
            # add hap tracks
            file = f"bb/hap_differences.bb"
            trackDb.write(HAP_TEMPLATE.format(file=file, sample=sample))
            file = "bb/FDR-wide-peaks.bb"
            trackDb.write(
                WIDE_TEMPLATE.format(file=file, name="FDR-wide-peaks", sample=sample)
            )

        # add percent accessible tracks
        file = f"bw/{hap}.percent.accessible.bw"
        if hap == "all":
            color = "0,0,0"
            trackDb.write(PER_ACC_COMP.format(sample=sample))
            # trackDb.write(FIRE_SCORE_COMP.format(sample=sample, file=f"bw/score.bw"))
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
            # zhap = "" if hap == "all" else f"_{hap}".replace("hap", "H")
            # trackDb.write(
            #    FIRE_SCORE.format(
            #        sample=sample,
            #        hap=hap,
            #        file=f"bw/score{zhap}.bw",
            #        viz=viz,
            #        color=color,
            #    )
            # )

        # new bin files
        if hap == "all":
            for z in ["H1", "H2", "UNK"]:
                trackDb.write(DECORATED.format(sample=sample, hap=z))

    # FDR scores
    trackDb.write(
        FIRE_SCORE_AND_FDR.format(
            sample=sample,
            fdr=f"bw/log_FDR.bw",
            score=f"bw/score.bw",
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
